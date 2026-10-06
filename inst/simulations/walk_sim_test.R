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
library(walk)  # remotes::install_github('dsjohnson/walk)

#' -----------------------------------------------------------------------------
#' Simulation parameters
#' -----------------------------------------------------------------------------
nsims <- 100
nbAnimals <- 1 # number of tracks
obsPerAnimal <- 100000 # number of simulated locations per track
timeStep <- 0.01 # time scale of simulation (should be small to help prevent discretization error)

models <- c("underdamped","overdamped")

betaList <- list(             # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
        c(-4, 6, 5, -0.1),    # High Gradient - Fast Diffusion
        c(-1, 2, 0, -0.1))    # Low Gradient - Slow Diffusion
ncov <- length(beta) - 1 # number of spatial covariates to be generated using simCov
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
          
          parMat <- matrix(NA,nrow=nsims,ncov+4+includeBarrier)
          colnames(parMat) <- c(paste0("beta",1:(ncov+1+includeBarrier)),"BA","SchoenerD","RMSE_log")
          
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
            walkRes = walkRes
          )
          
          # Collapse all parameters into key=value pairs separated by underscores
          param_tags <- paste(names(cur_params), cur_params, sep = "=", collapse = "_")
          
          # Construct the fileName
          fileName <- paste0("inst/simulations/results/sim_results_", param_tags, ".RData")
          
          print(cur_params)
          
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
                spatialCovs[[isim]][[i]] <- simCov(sca = sca, res = langRes, irange=irange, sigma2 = 0.1, kappa = 0.5)
                
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
                  
                  # simulate complex coastline and islands
                  r <- terra::rast(nrows = nrow(spatialCovs[[isim]]$cov1), ncols = ncol(spatialCovs[[isim]]$cov1), ext = ext(spatialCovs[[isim]]$cov1))
                  crs(r) <- ""
                  terra::values(r) <- runif(terra::ncell(r))
                  w_size <- sample(seq(11, 21, by = 2), 1)
                  w <- matrix(1, nrow = w_size, ncol = w_size)
                  noise <- r
                  for (i in 1:4) {
                    noise <- terra::focal(noise, w = w, fun = mean, na.policy = "omit", expand=TRUE)
                  }
                  n_min <- terra::global(noise, "min", na.rm = TRUE)[[1]]
                  n_max <- terra::global(noise, "max", na.rm = TRUE)[[1]]
                  noise <- (noise - n_min) / (n_max - n_min)
                  x_coords <- terra::init(noise, "x")
                  base_width <- runif(1, 0.10, 0.25)
                  wiggle_amp <- runif(1, 0.15, 0.35)
                  x_min <- terra::xmin(r)
                  x_range <- terra::xmax(r) - x_min
                  mainland <- x_coords < (x_min + x_range * base_width + noise * x_range * wiggle_amp)
                  isl_thresh <- runif(1, 0.65, 0.85)
                  islands <- noise > isl_thresh
                  land_mask <- mainland | islands
                  water_mask <- terra::ifel(land_mask, 0, 1)
                  names(water_mask) <- "coast_barrier"
                  
                  barrier <- suppressMessages(prepBarrier(water_mask))
                  water_mask <- ifel(water_mask>0,1,NA)
                  
                  spatialCovs[[isim]]$d2coast <- barrier/100
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
              obsExt <- terra::ext(cbind(c(simLang[[isim]]$x,simLang[[isim]]$mu.x),c(simLang[[isim]]$y,simLang[[isim]]$mu.y)))
              subDatWalk[[isim]]$dt <- c(0, diff(subDatWalk[[isim]]$date))
              cropCovs <- lapply(spatialCovs[[isim]],function(x) crop(x,obsExt+obsBuff))
              
              if(includeBarrier) {
                water_mask <- crop(water_mask,obsExt+obsBuff)
              }
              
              langFit[[isim]] <- tryCatch(suppressMessages(fitLangevin(subDatWalk[[isim]],model=model,spatialCovs=spatialCovs[[isim]],barrier=barrier)),error=function(e) e)
              if(inherits(langFit[[isim]],"error")) next
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
                maskRast <- water_mask # Already cropped above
                # Coarsen the barrier mask preserving the barrier (use max)
                if(aggFactWalk > 1) {
                  maskRast <- terra::aggregate(maskRast, fact = aggFactWalk, fun = max, na.rm = TRUE)
                } else if (aggFactWalk < 1) {
                  maskRast <- terra::disagg(maskRast, fact = 1/aggFactWalk)
                }
              }
              spatialCovsWalk[[isim]] <- rast(cropCovsWalk)
              pt_data <- tryCatch(proc_simLangevin(subDatWalk[[isim]], spatialCovsWalk[[isim]]),error=function(e) e)
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
              
              # Create differences
              walk_data$q_m$d_cov1 <-  walk_data$q_m$cov1 - walk_data$q_m$from_cov1
              walk_data$q_m$d_cov2 <-  walk_data$q_m$cov2 - walk_data$q_m$from_cov2
              walk_data$q_m$d_cov3 <-  walk_data$q_m$cov3 - walk_data$q_m$from_cov3
              walk_data$q_m$d_d2c <-  walk_data$q_m$d2c - walk_data$q_m$from_d2c
              
              form <- ~d_cov1 + d_cov2 + d_cov3 + d_d2c
              
              if(includeBarrier){
                walk_data$q_m$d_d2coast <-  walk_data$q_m$d2coast - walk_data$q_m$from_d2coast
                form <- ~d_cov1 + d_cov2 + d_cov3 + d_d2c + d_d2coast
              }
              
              model_parameters <- ctmc_control(
                q_r = ctmc_model(form = ~ 1, link="log"), 
                q_m = ctmc_model(form= form, link="log"),
                # q_m = ctmc_model(form= ~1, link="log"),
                norm=FALSE, delta="uniform")
              
              walkFit[[isim]] <- tryCatch(suppressMessages(fit_ctmc(walk_data, model_parameters = model_parameters, 
                                                                    start=list(beta_q_r=sigma/2)#,control=list(trace=1)
              )),error=function(e) e)
              # beta <- c(2*walkFit[[isim]]$results$q_m$est)
              
              if(!inherits(walkFit[[isim]],"error")){
                # Get UD
                ud <- tryCatch(get_lim_ud(walkFit[[isim]]),error=function(e) e)
                
                if(!inherits(ud,"error")){
                  # Get parameter estimates
                  beta_q_m <- walkFit[[isim]]$results$q_m$est
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
                    
                    # 5. VISUAL FIX: Clamp 0s so log(0) doesn't plot as transparent white space
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
                      labs(subtitle="Estimated UD", title=paste("Simulation",isim)) + 
                      ggplot2::theme_minimal() +
                      ggplot2::theme(
                        axis.title.x = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank() 
                      ) + geom_point(aes(x=x,y=y),data=subDatWalk[[isim]],col="orange")
                    
                    
                    p2 <- ggplot2::ggplot() +
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
                      ) + geom_path(aes(x=mu.x,y=mu.y),data=simLang[[isim]],col="orange")
                    
                    print(p1 + p2 + patchwork::plot_layout(ncol = 1, nrow = 2, guides = "collect"))
                    
                    parMat[isim,1:(ncov+1+includeBarrier)] <- beta_q_m
                    parMat[isim,"BA"] <- walkMetrics["BA"] 
                    parMat[isim,"SchoenerD"] <- walkMetrics["SchoenerD"] 
                    parMat[isim,"RMSE_log"] <- walkMetrics["RMSE_log"] 
                    
                    cat("current:\n")
                    print(parMat[isim,])
                    cat("\n overall:\n")
                    print(apply(parMat,2,mean,na.rm=TRUE))
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
          
          save(parMat, parMatLang, file = fileName)
          
          print(apply(parMatLang,2,mean))
        }
      }
    }
  }
}
