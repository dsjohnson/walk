library(tidyverse)
library(terra)
library(viridis)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(sf)
library(Rcpp)
library(fields)
library(future)
library(furrr)
library(langevinSSM) # remotes::install_github('bmcclintock/langevinSSM)
library(here)
library(walk)  # remotes::install_github('bmcclintock/walk@speed-up-slim', subdir="walk")

#' -----------------------------------------------------------------------------
#' Simulation parameters
#' -----------------------------------------------------------------------------
nsims <- 100
nbAnimals <- 5 # number of tracks
obsPerAnimal <- 100000 # number of simulated locations per track
timeStep <- 0.01 # time scale of simulation (should be small to help prevent discretization error)

models <- c("underdamped","overdamped")

betaList <- list(             # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
  c(-4, 6, 5, -0.1),    # High Gradient - Fast Diffusion
  c(-1, 2, 0, -0.1))    # Low Gradient - Slow Diffusion

sigmaList <- list(      # diffusion (or speed) parameter
  2.5, # Fast Diffusion
  1)   # Slow Diffusion

gamma <- 0.5 # friction parameter (smaller value -> more directional persistence); ignored unless model=="underdamped"
psi <- 1 # error ellipse scaling parameter

includeBarrier <- TRUE
if(includeBarrier) betaList <- lapply(betaList,function(x) c(x,-0.1)) # add d2coast covariate

## sampling rate, missing data, and measurement error
samplingRate <- 1 # for subsampling observations from true continuous-time model (e.g. if samplingRate = 2 then data are roughly thinned by 2); must be >= 1; note bias increases with samplingRate, but relationships largely preserved
propMissingVec <- c(0.99,0.995,0.998)

smaj.sdVec <- c(0.5,2.5,5,10) # SD for semi-major error ellipse axis; smaj ~ abs(Normal(0,smaj.sd))
eor.lim <- c(0,180) # range for error ellipse orientation (in degrees from north); eor ~ Uniform(eor.lim[1],eor.lim[2])

## specify scale and spatial autocorrelation for covariates
sca <- 200 # bounding box scale for Langevin diffusion
obsBuff <- 10 # buffer around observed locations for defining raster extent for {walk}
covRange <- c(0.1,0.5) # lower and upper bounds for covariate spatial range parameter (lower has less spatial autocorrelation)

langRes <- 1 # cell resolution for Langevin diffusion simulation and fitting
walkResVec <- c(1,2,5,8) # cell resolution for walk model fitting (should be >= langRes)

nthreads <- nbAnimals

#' -----------------------------------------------------------------------------
#' walk model form (requires walk with walk_logistic_form.patch applied)
#' -----------------------------------------------------------------------------
#' "mult"     : q_ij = exp(b0) * exp((x_j - x_i)'b)            UD ~ exp(2 x'b)        => beta = 2 b
#' "sde"      : q_ij = exp(b0)/(k h^2) + f(mu_ij)/((k/2) h),    mu_ij = b'(x_j - x_i)/h (Langevin drift along the edge)
#'                                                              UD ~ exp(2 x'b/exp(b0)) => beta = 2 b / exp(b0)
#' "logistic" : q_ij = lambda * logistic((x_j - x_i)'b)         UD ~ exp(x'b)/lambda    => beta = b
#'              lambda = U*logistic(a0) (link "logit", bounded by U) or exp(a0) (link "log")
#' 'beta' columns of parMat are on the Langevin (UD) scale for every form, so they are directly
#' comparable with the true beta and with parMatLang; the raw walk estimates are kept in parMatRaw.
walkForm <- "mult"          # one of "mult", "sde", "logistic"
walkForm <- match.arg(walkForm, c("mult", "sde", "logistic"))
sdeA <- Inf                 # "sde" only: drift transform; Inf = hard plus max(0,x) (default), finite = soft plus (see ?ctmc_control)
logisticLink <- "logit"     # "logistic" only: residency link, "logit" (rate bounded by U) or "log" (unbounded)
logisticUmult <- 10         # "logistic" + "logit" only: U = logisticUmult * (2 D / h^2), where 2 D / h^2 is the value of
# lambda implied by the true diffusion D (see below). U is fixed, not estimated; with real data
# set it from an upper bound on the diffusion rate instead.
walkClip <- 0               # optional cap on the total leaving rate of a cell (0 = none), see ?ctmc_control

## coastline placement (see simLand.R): the coast passes coastDist*udSD from the UD centre (origin)
source("simulations/simLand.R")
coastDist <- c(0.3, 1)   # in units of the d2c UD scale sqrt(sca / (2 |beta_d2c|))
landBuff <- 3            # walk mask keeps land cells within landBuff (coordinate units) of the coast

rasterMetrics <- function(r1, r2) {
  
  # Ensure NAs match perfectly before calculating
  overlap_mask <- terra::mask(r1, r2)
  r2 <- terra::mask(r2, overlap_mask)
  
  # Extract values as vectors, dropping NAs
  v1 <- terra::values(overlap_mask, mat = FALSE, na.rm = TRUE)
  v2 <- terra::values(r2, mat = FALSE, na.rm = TRUE)
  
  if (length(v1) == 0 || length(v2) == 0) return(c(BA = NA, SchoenerD = NA, RMSE_log = NA))
  
  # Normalize to true probability distributions (sum to 1)
  p <- v1 / sum(v1)
  q <- v2 / sum(v2)
  
  # Bhattacharyya's Affinity (Pattern / Colocation)
  BA <- sum(sqrt(p * q))
  
  # Schoener's D (Absolute Probability Intersection)
  SchoenerD <- 1 - 0.5 * sum(abs(p - q))
  
  # RMSE of log(UD) (Magnitude / Intensity deviation)
  # Prevent log(0) by filtering out 0s if any exist
  valid_idx <- p > 0 & q > 0
  RMSE_log <- sqrt(mean((log(p[valid_idx]) - log(q[valid_idx]))^2))
  
  return(c(BA = BA, SchoenerD = SchoenerD, RMSE_log = RMSE_log))
}

for(ind in 1:length(sigmaList)){
  
  sigma <- sigmaList[[ind]]
  beta <- betaList[[ind]]
  
  ncov <- length(beta) - 1 - includeBarrier # number of spatial covariates to be generated using simCov
  
  for(model in models){
    
    #' -----------------------------------------------------------------------------
    #' Discretization & Stability Checks
    #' -----------------------------------------------------------------------------
    if (model == "underdamped") {
      # Expected distance traveled in one timeStep
      exp_step <- (sigma / sqrt(2 * gamma)) * timeStep
      
      # Temporal stability check: gamma * dt should be << 1
      if ((gamma * timeStep) > 0.1) {
        warning(sprintf("Temporal instability risk: gamma * timeStep (%.3f) is large. Euler integration may be unstable. Consider decreasing timeStep.", 
                        gamma * timeStep),immediate. = TRUE)
      }
    } else if (model == "overdamped") {
      # Expected distance traveled in one timeStep (2D Brownian motion)
      exp_step <- sigma * sqrt((pi / 2) * timeStep)
    }
    
    # Spatial sampling check: step length should be << langRes
    step_ratio <- exp_step / langRes
    
    if (step_ratio > 0.25) {
      warning(sprintf(
        "Expected step length (%.4f) is large relative to langRes (%g). Consider decreasing timeStep or increasing langRes.", 
        exp_step, langRes
      ), immediate. = TRUE)
    }
    
    for(propMissing in propMissingVec){
      
      for(smaj.sd in smaj.sdVec){
        
        smin.sd <- smaj.sd/2 # SD for semi-minor error ellipse axis; smin ~ abs(Normal(0,smin.sd))
        
        measurementError <- list(smaj.sd=smaj.sd,smin.sd=smin.sd,eor.lim=eor.lim) 
        
        if(smaj.sd==2.5) resInd <- walkResVec
        else resInd <- 2
        
        for(walkRes in resInd){
          
          simLang <- subDatWalk <- spatialCovs <- spatialCovsWalk <- langFit <- walkFit <- vector("list",nsims)
          
          parMat <- matrix(NA,nrow=nsims,ncov+5+includeBarrier)
          colnames(parMat) <- c(paste0("beta",1:(ncov+1+includeBarrier)),"sigma_eq","BA","SchoenerD","RMSE_log")
          # raw walk estimates (beta_q_r, beta_q_m) on the scale of the walk model
          parMatRaw <- matrix(NA,nrow=nsims,ncov+2+includeBarrier)
          colnames(parMatRaw) <- c("beta_q_r",paste0("beta_q_m",1:(ncov+1+includeBarrier)))
          
          parMatLang <- matrix(NA,nrow=nsims,ncov+5+1+includeBarrier)
          colnames(parMatLang) <- c(paste0("beta",1:(ncov+1+includeBarrier)),"sigma","gamma","BA","SchoenerD","RMSE_log")
          
          cur_params <- list(
            model = model,
            nsims = nsims,
            nbAnimals = nbAnimals,
            obsPerAnimal = obsPerAnimal,
            timeStep = timeStep,
            sigma = sigma,
            gamma = gamma,
            psi = psi,
            includeBarrier = includeBarrier,
            samplingRate = samplingRate,
            propMissing = propMissing,
            smaj.sd = smaj.sd,
            sca = sca,
            langRes = langRes,
            walkRes = walkRes,
            walkForm = walkForm
          )
          if(walkForm == "sde" && is.finite(sdeA)) cur_params$sdeA <- sdeA
          if(walkForm == "logistic"){
            cur_params$logisticLink <- logisticLink
            if(logisticLink == "logit") cur_params$logisticUmult <- logisticUmult
          }
          if(walkClip > 0) cur_params$walkClip <- walkClip
          
          # Collapse all parameters into key=value pairs separated by underscores
          param_tags <- paste(names(cur_params), cur_params, sep = "=", collapse = "_")
          
          # Construct the fileName
          fileName <- paste0("simulations/results/sim_results_", param_tags, ".RData")
          
          print(cur_params)
          
          if(!file.exists(fileName)){
            
            #' -----------------------------------------------------------------------------
            #' Begin simulation here
            #' -----------------------------------------------------------------------------
            set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
            for(isim in 1:nsims){
              
              message("Simulation ",isim)
              
              if(walkRes < langRes) stop("walkRes must be >= langRes)")
              
              # Simulate path
              if(model=="underdamped"){
                par <- list(beta=beta,sigma=sigma,gamma=gamma,psi=psi)
              } else {
                par <- list(beta=beta,sigma=sigma,psi=psi)
              }
              
              # simulate "high resolution" tracks
              simLang[[isim]] <- walkFit[[isim]] <- tryCatch(stop(),error=function(e) e)
              while(inherits(simLang[[isim]],"error") | inherits(walkFit[[isim]],"error")){
                simLang[[isim]] <- walkFit[[isim]] <- tryCatch(stop(),error=function(e) e)
                
                # Simulate spatial habitat
                covNames <- c("cov1","cov2","cov3","d2c")
                spatialCovs[[isim]] <- list()
                
                for(i in 1:ncov) {
                  spatialCovs[[isim]][[i]] <- list()
                  irange <- runif(1,covRange[1],covRange[2])
                  
                  # Generate the covariate natively at langRes by passing it to simCov
                  spatialCovs[[isim]][[i]] <- simCov(sca = sca, res = langRes, irange=irange, sigma2 = 0.1, kappa = 1.5)
                  
                  names(spatialCovs[[isim]][[i]]) <- paste0("cov",i)
                  # terra::crs(spatialCovs[[isim]][[i]]) <- "epsg:3416"
                }
                
                coords <- terra::crds(spatialCovs[[isim]][[1]])
                dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca
                
                spatialCovs[[isim]][[4]] <- terra::setValues(spatialCovs[[isim]][[1]][[1]], dist2)
                names(spatialCovs[[isim]][[4]]) <- covNames[4]
                names(spatialCovs[[isim]]) <- covNames
                
                while(inherits(simLang[[isim]],"error")){
                  
                  if(includeBarrier){
                    
                    # simulate coastline (mainland + islands) through the animal's range
                    udSD <- sqrt(sca / (2 * abs(beta[ncov + 1])))   # d2c term: UD ~ N(0, udSD^2) per axis
                    water_mask <- simLand(spatialCovs[[isim]]$cov1, udSD = udSD, coastDist = coastDist)
                    
                    barrier <- suppressMessages(prepBarrier(water_mask))
                    water_mask <- ifel(water_mask>0,1,NA)
                    
                    spatialCovs[[isim]]$d2coast <- barrier
                    spatialCovs[[isim]]$coast_barrier = barrier
                    
                    barrier <- "coast_barrier"
                    
                  } else {
                    barrier <- NULL
                    water_mask <- NULL
                  }
                  
                  simLang[[isim]] <- tryCatch({
                    out_data <- suppressMessages(simLangevin(
                      model = model,
                      nbAnimals = nbAnimals,
                      obsPerAnimal = obsPerAnimal,
                      timeStep = timeStep,
                      measurementError = measurementError,
                      par=par,
                      barrier = barrier,
                      spatialCovs = spatialCovs[[isim]],
                      subSample = list(samplingRate = samplingRate, propMissing = propMissing)
                    ))
                    
                    if(includeBarrier){  
                      pts <- cbind(out_data$x, out_data$y)
                      dist_vals <- terra::extract(spatialCovs[[isim]]$coast_barrier, pts)[, 1]
                      land_idx <- which(dist_vals <= 0)
                      
                      if (length(land_idx) == 0) stop("No observed locations on land. Resimulating...")
                    }
                    
                    out_data
                  },error=function(e) e)
                  if(inherits(simLang[[isim]],"error")){
                    message("    Retrying Simulation ",isim,": ",simLang[[isim]]$message)
                  }
                }
                
                subDatWalk[[isim]] <- filter(simLang[[isim]], !is.na(x))
                
                # Format for {walk}
                obsXY <- cbind(c(simLang[[isim]]$x,simLang[[isim]]$mu.x),c(simLang[[isim]]$y,simLang[[isim]]$mu.y))
                obsXY <- obsXY[stats::complete.cases(obsXY), , drop=FALSE]
                obsExt <- terra::ext(min(obsXY[,1]), max(obsXY[,1]), min(obsXY[,2]), max(obsXY[,2]))
                subDatWalk[[isim]]$dt <- c(0, diff(subDatWalk[[isim]]$date))
                cropCovs <- lapply(spatialCovs[[isim]],function(x) crop(x,obsExt+obsBuff))
                
                if(includeBarrier) {
                  water_mask <- crop(water_mask,obsExt+obsBuff)
                }
                
                langFit[[isim]] <- tryCatch(suppressMessages(fitLangevin(subDatWalk[[isim]],model=model,spatialCovs=spatialCovs[[isim]],barrier=barrier)),error=function(e) e)
                if(inherits(langFit[[isim]],"error")) next
                if(langFit[[isim]]$convergence!=0){
                  tmp <- tryCatch(suppressMessages(fitLangevin(subDatWalk[[isim]],model=model,par=par,map=list(psi=factor(NA)),spatialCovs=spatialCovs[[isim]],barrier=barrier)),error=function(e) e)
                  if(inherits(tmp,"error")) next
                  if(tmp$convergence!=0){
                    message("    Retrying Simulation ",isim,": LangFit did not converge (convergence code: ",langFit[[isim]]$convergence,")")
                    next
                  } else langFit[[isim]] <- tmp
                }
                langPar <- getPar(langFit[[isim]])
                parMatLang[isim,1:(ncov+1+includeBarrier)] <- getPar(langFit[[isim]])$beta
                parMatLang[isim,"sigma"] <- getPar(langFit[[isim]])$sigma
                if(model=="underdamped") parMatLang[isim,"gamma"] <- getPar(langFit[[isim]])$gamma
                
                trueUD <- getUD(cropCovs,beta=par$beta,barrier=barrier,lambda=attr(simLang[[isim]],"lambda"),plot=FALSE,maskRast=water_mask,log=TRUE)
                langUD <- suppressMessages(getUD(cropCovs,beta=getPar(langFit[[isim]])$beta,barrier=barrier,lambda=attr(simLang[[isim]],"lambda"),maskRast=water_mask,log=TRUE,plot=FALSE))
                
                langMetrics <- rasterMetrics(exp(langUD$log_UD),exp(trueUD))
                parMatLang[isim,"BA"] <- langMetrics["BA"] 
                parMatLang[isim,"SchoenerD"] <- langMetrics["SchoenerD"] 
                parMatLang[isim,"RMSE_log"] <- langMetrics["RMSE_log"] 
                
                cropCovs$coast_barrier <- NULL
                
                # Aggregation step: coarsen the rasters for the walk model relative to langRes
                aggFactWalk <- walkRes / langRes
                if(aggFactWalk > 1) {
                  cropCovsWalk <- lapply(cropCovs, function(x) terra::aggregate(x, fact = aggFactWalk, fun = mean, na.rm = TRUE))
                } else if (aggFactWalk < 1) {
                  cropCovsWalk <- lapply(cropCovs, function(x) terra::disagg(x, fact = 1/aggFactWalk))
                } else {
                  cropCovsWalk <- cropCovs
                }
                
                maskRast <- NULL
                if(includeBarrier){
                  # Allowed cells for walk: water plus land within landBuff of the coast. The barrier is
                  # soft, so true locations reach ~1.5 cells inland; masking all land then gives
                  # observations with zero likelihood in every cell.
                  maskRast <- crop(ifel(spatialCovs[[isim]]$coast_barrier > -landBuff, 1, NA), obsExt+obsBuff)
                  # Coarsen the barrier mask preserving the barrier (use max)
                  if(aggFactWalk > 1) {
                    maskRast <- terra::aggregate(maskRast, fact = aggFactWalk, fun = max, na.rm = TRUE)
                  } else if (aggFactWalk < 1) {
                    maskRast <- terra::disagg(maskRast, fact = 1/aggFactWalk)
                  }
                }
                spatialCovsWalk[[isim]] <- rast(cropCovsWalk)
                pt_data <- tryCatch(proc_Langevin(subDatWalk[[isim]], spatialCovsWalk[[isim]],neighborhood="buffer"),error=function(e) e)
                
                if(inherits(pt_data,"error")){
                  message("    Retrying Simulation ",isim,": ",pt_data$message)
                  next
                }
                
                walk_data <- tryCatch({
                  withCallingHandlers(
                    make_walk_data(pt_data, spatialCovsWalk[[isim]], rast_mask=maskRast),
                    warning = function(w) {
                      if (grepl("likelihood values of 0 in all ras cells", w$message)) {
                        stop(w$message) 
                      }
                    }
                  )
                }, error = function(e) e)
                
                # Standard retry logic to catch the new error and trigger the while loop to resimulate
                if (inherits(walk_data, "error")) {
                  message("    Retrying Simulation ", isim, ": ", walk_data$message)
                  next
                }
                
                # Covariates of the potential (beta order: cov1..cov3, d2c[, d2coast])
                potCovs <- c("cov1","cov2","cov3","d2c")
                if(includeBarrier) potCovs <- c(potCovs, "d2coast")
                
                # Long-run diffusion coefficient D of the simulated Langevin process (per axis, position units^2 / time):
                # overdamped D = sigma^2/2; underdamped D = sigma^2/gamma. Used only for starting values
                # (and for U of the bounded logistic form).
                D_true <- if(model == "underdamped") sigma^2/gamma else sigma^2/2
                h_walk <- walkRes
                
                if(walkForm == "mult"){
                  # differences x_j - x_i on each edge
                  for(v in potCovs) walk_data$q_m[[paste0("d_", v)]] <- walk_data$q_m[[v]] - walk_data$q_m[[paste0("from_", v)]]
                  qmForm <- reformulate(paste0("d_", potCovs))
                  model_parameters <- ctmc_control(
                    q_r = ctmc_model(form = ~ 1, link="log"), 
                    q_m = ctmc_model(form= qmForm, link="log"),
                    norm=FALSE, delta="uniform", form="mult", clip=walkClip)
                  startWalk <- list(beta_q_r=sigma/2)
                } else if(walkForm == "sde"){
                  # directional gradient of each covariate along the edge, (x_j - x_i)/h
                  # (q_m formula must have no intercept; q_r must be ~1)
                  for(v in potCovs) walk_data$q_m[[paste0("g_", v)]] <- (walk_data$q_m[[v]] - walk_data$q_m[[paste0("from_", v)]])/walk_data$q_m$h
                  qmForm <- reformulate(c(paste0("g_", potCovs), "-1"))
                  model_parameters <- ctmc_control(
                    q_r = ctmc_model(form = ~ 1, link="log"),        # exp(beta_q_r) = 2 D (= sigma^2 overdamped)
                    q_m = ctmc_model(form= qmForm),                  # beta_q_m = D * beta
                    delta="uniform", form="sde", sde_a=sdeA, clip=walkClip)
                  startWalk <- list(beta_q_r=log(2*D_true))
                } else {
                  # cell covariates (columns of walk_data$q_r); the edge design x_j - x_i is built by fit_ctmc
                  qmForm <- reformulate(potCovs)
                  lam0 <- 2*D_true/h_walk^2    # lambda giving diffusion D: each of the 4 rook rates is lambda/2 at b = 0
                  if(logisticLink == "logit"){
                    U <- logisticUmult*lam0
                    q_r_mod <- ctmc_model(form = ~ 1, link="logit", U=U)
                    startWalk <- list(beta_q_r=qlogis(lam0/U))
                  } else {
                    q_r_mod <- ctmc_model(form = ~ 1, link="log")
                    startWalk <- list(beta_q_r=log(lam0))
                  }
                  model_parameters <- ctmc_control(
                    q_r = q_r_mod,
                    q_m = ctmc_model(form= qmForm),                  # link is set to logistic by form="logistic"
                    delta="uniform", form="logistic", clip=walkClip)
                }
                
                walkFit[[isim]] <- tryCatch(suppressMessages(fit_ctmc(walk_data, model_parameters = model_parameters, nthreads = nthreads,
                                                                      start=startWalk#,control=list(trace=1)
                )),error=function(e) e)
                # beta <- c(2*walkFit[[isim]]$results$q_m$est)
                
                if(!inherits(walkFit[[isim]],"error")){
                  # Get UD
                  ud <- tryCatch(get_lim_ud(walkFit[[isim]]),error=function(e) e)
                  
                  if(!inherits(ud,"error")){
                    # Get parameter estimates
                    beta_q_r <- walkFit[[isim]]$results$q_r$est
                    beta_q_m <- walkFit[[isim]]$results$q_m$est
                    # convert to the Langevin scale: UD ~ exp(x'beta); sigma_eq = sqrt(2 D_hat) is the
                    # overdamped-equivalent sigma implied by the fitted diffusion rate
                    if(walkForm == "mult"){
                      beta_walk <- 2*beta_q_m
                      sigma_eq <- h_walk*sqrt(2*exp(beta_q_r))
                    } else if(walkForm == "sde"){
                      beta_walk <- 2*beta_q_m/exp(beta_q_r)
                      sigma_eq <- sqrt(exp(beta_q_r))
                    } else {
                      lam_hat <- if(logisticLink == "logit") U*plogis(beta_q_r) else exp(beta_q_r)
                      beta_walk <- beta_q_m
                      sigma_eq <- h_walk*sqrt(lam_hat)
                    }
                    # Calculate preference function from parameters
                    #ud$pref <- beta_q_m[1]*walk_data$q_r$cov1 + beta_q_m[2]*walk_data$q_r$cov2 + beta_q_m[3]*walk_data$q_r$cov3 +  beta_q_m[4]*walk_data$q_r$d2c
                    #if(includeBarrier) ud$pref <- ud$pref + beta_q_m[4]*walk_data$q_r$d2c
                    # Get movement-rate (generator) matrix
                    Q <- get_Q(walkFit[[isim]])
                    # Get expected residence time from diagonal of movement-rate matrix
                    #ud$res <- -1/diag(Q)
                    
                    cells <- rast(spatialCovsWalk[[isim]]$d2c)
                    cells[['UD']] <- 0
                    cells[ud$cell] <- ud$ud
                    exptrueUD <- exp(trueUD$log_UD)
                    
                    # convert coarse NAs to 0
                    coarse_ud <- cells$UD
                    coarse_ud[is.na(coarse_ud)] <- 0
                    
                    # expand the coarse bounding box to match the True UD before resampling
                    coarse_ud <- terra::extend(coarse_ud, exptrueUD, fill = 0)
                    
                    # resample 
                    est_ud_resampled <- terra::resample(coarse_ud, exptrueUD, method = "bilinear")
                    
                    if(includeBarrier) {
                      est_ud_resampled <- terra::mask(est_ud_resampled, water_mask) 
                    }
                    
                    # renormalize 
                    est_ud_sum <- terra::global(est_ud_resampled, "sum", na.rm = TRUE)[[1]]
                    est_ud_resampled <- est_ud_resampled / est_ud_sum
                    
                    # Calculate metrics safely
                    walkMetrics <- tryCatch(rasterMetrics(est_ud_resampled, exptrueUD), error=function(e) e)
                    
                    if(!inherits(walkMetrics,"error")){
                      
                      if(is.na(walkMetrics["BA"])) {
                        walkFit[[isim]] <- tryCatch(stop(),error=function(e) e)
                        message("    Retrying Simulation ",isim,": Raster metrics evaluated to NA.")
                        next
                      }
                      
                      plot_ud <- est_ud_resampled
                      plot_ud[plot_ud == 0] <- min(values(exptrueUD), na.rm = TRUE)
                      
                      shared_limits <- range(c(log(values(plot_ud)), values(trueUD$log_UD)), na.rm = TRUE)
                      
                      p1 <- ggplot2::ggplot() +
                        ggplot2::geom_raster(data = log(plot_ud), ggplot2::aes(x = x, y = y, fill = UD)) +
                        ggplot2::scale_fill_viridis_c(
                          name = "log(UD)", 
                          option = "viridis", 
                          na.value = "transparent",
                          limits = shared_limits 
                        ) + ggplot2::coord_equal() +
                        labs(subtitle=paste0("walk UD (", walkForm, ")"), title=paste("Simulation",isim)) + 
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                          axis.title.x = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank() 
                        ) + geom_path(aes(x=x,y=y, group = id),data=subDatWalk[[isim]],col="orange")
                      
                      p2 <- ggplot2::ggplot() +
                        ggplot2::geom_raster(data = langUD, ggplot2::aes(x = x, y = y, fill = log_UD)) +
                        ggplot2::scale_fill_viridis_c(
                          name = "log(UD)", 
                          option = "viridis", 
                          na.value = "transparent",
                          limits = shared_limits 
                        ) + ggplot2::coord_equal() +
                        labs(subtitle="langevinSSM UD") + 
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                          axis.title.x = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank() 
                        ) + geom_path(aes(x=x,y=y, group = id),data=subDatWalk[[isim]],col="orange")
                      
                      
                      p3 <- ggplot2::ggplot() +
                        ggplot2::geom_raster(data = trueUD$log_UD, ggplot2::aes(x = x, y = y, fill = log_UD)) +
                        ggplot2::scale_fill_viridis_c(
                          name = "log(UD)", 
                          option = "viridis", 
                          na.value = "transparent",
                          limits = shared_limits 
                        ) + ggplot2::coord_equal() +
                        labs(subtitle="True UD") + 
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                          axis.title.x = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank() 
                        ) + geom_path(aes(x=mu.x,y=mu.y, group=id),data=simLang[[isim]],col="orange")
                      
                      
                      print(p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, nrow = 1, guides = "collect"))
                      
                      parMat[isim,1:(ncov+1+includeBarrier)] <- beta_walk
                      parMat[isim,"sigma_eq"] <- sigma_eq
                      parMatRaw[isim,] <- c(beta_q_r, beta_q_m)
                      parMat[isim,"BA"] <- walkMetrics["BA"] 
                      parMat[isim,"SchoenerD"] <- walkMetrics["SchoenerD"] 
                      parMat[isim,"RMSE_log"] <- walkMetrics["RMSE_log"] 
                      
                      cat("current:\n")
                      print(parMat[isim,])
                      cat("\n overall:\n")
                      print(apply(parMat,2,mean,na.rm=TRUE))
                      
                      cat("current:\n")
                      print(parMatLang[isim,])
                      cat("\n overall:\n")
                      print(apply(parMatLang,2,mean,na.rm=TRUE))
                      
                    } else {
                      walkFit[[isim]] <- tryCatch(stop(),error=function(e) e)
                      message("    Retrying Simulation ",isim,": ",walkMetrics$message)
                      next
                    }
                  } else {
                    walkFit[[isim]] <- tryCatch(stop(),error=function(e) e)
                    message("    Retrying Simulation ",isim,": ",ud$message)
                    next
                  }
                } else {
                  message("    Retrying Simulation ",isim,": ",walkFit[[isim]]$message)
                }
              }
            }
            
            save(parMat, parMatRaw, parMatLang, cur_params, file = fileName)
          } else {
            load(fileName)
            print(apply(parMat,2,mean,na.rm=TRUE))
            print(apply(parMatLang,2,mean,na.rm=TRUE))
          }
        }
      }
    }
  }
}
