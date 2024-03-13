library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(Matrix)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(lubridate)
library(animation)
#library(walk)
devtools::load_all("~/research/projects/r_packages/walk")

setwd("~/research/projects/r_packages/walk/testing/analysis/cod")

# load data
bathy_terra <- readRDS("ai_bathy_3km.Rds")
likelihood_terra <- readRDS("likelihood_3km.Rds")

# Get land layer
ak <- ne_states( "united states of america", return="sf")
ak <- ak[pmatch("Alas", ak$name_en),]
ak <- st_transform( ak, crs=st_crs("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") )
ak <- st_geometry(ak)

# change resolution
fact <- 2
bathy_terra <- terra::aggregate( bathy_terra, fact = fact )
likelihood_terra <- as.array( terra::aggregate(likelihood_terra, fact=fact) )

# Eliminate land cells
bathy_terra[bathy_terra<=1] <- NA
bathy_terra$bathy_km <- bathy_terra/1000
bathy_terra$bathy_km2 <- bathy_terra$bathy_km^2
bathy_terra$bathy_t <- terra::ifel(bathy_terra$bathy_km>3, 3, bathy_terra$bathy_km)

ggplot() +
  geom_spatraster(data=bathy_terra$bathy_km) + 
  scale_fill_gradientn(colors=rev(sf.colors(10)), na.value = "transparent", name="Bathymetry (km)") + 
  # scale_fill_distiller(palette="YlOrRd", direction = 1, na.value = "transparent", name="UD") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) +
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()




# Format for observed data for {walk}
# list(
# `L` = sparse matrix with rows as observed likelihood values
# `times` = data.frame with columns (1) "obs" = observation number (2) "timestamp" = POSIX
# time of observation, and (3) "dt"  = diff(timestamp) with a 0 in the first element.
#)
# If you want to predict locations then you need to add the units of "dt" as an attribute.

L <- t(apply(likelihood_terra, MARGIN=3, FUN=function(mat){as.vector(t(mat))}))
L <- Matrix(L) |> as("CsparseMatrix")
pdata <- list(L=L)
pdata$times <- data.frame(
  obs=1:nrow(L), 
  timestamp = seq(mdy_hms("2/21/2019 12:00:00"), mdy_hms("5/23/2019 12:00:00"), "1 day")
  )
pdata$times$dt = c(0, diff(pdata$times$timestamp)) 
units(pdata$times$dt) <- "days"
attr(pdata$times,"time_unit") <- "days"


# Format habitat for {walk}
walk_data <- make_walk_data(pdata, bathy_terra, grad=c("bathy_km","bathy_km2"), rast_mask = bathy_terra[[1]]) 


### Bathymetry differences for preference function.
walk_data$q_m$d_bathy <- (walk_data$q_m$bathy_km-walk_data$q_m$from_bathy_km) / (res(bathy_terra)[1]/1000)
walk_data$q_m$d_bathy2 <- (walk_data$q_m$bathy_km^2-walk_data$q_m$from_bathy_km^2) / (res(bathy_terra)[1]/1000)


# Diffusion depends on average d_bathy more like traditional CTMC
fit <- fit_ctmc(walk_data,
                 model_parameters = list(
                   q_r = list(form=~ bathy_km), # residency model
                   q_m = list(form=~ d_bathy + d_bathy2) # preference model
                   ), 
                reals = FALSE, fit=TRUE, control=list(trace=1)
                )
ud <- get_lim_ud(fit)
beta_q_m <- fit$results$beta$q_m$est
ud$pref <- beta_q_m[1]*walk_data$q_r$bathy_km + beta_q_m[2]*walk_data$q_r$bathy_km^2

Q <- get_Q(fit)

ud$res <- -1/diag(Q)
plot(walk_data$q_r$bathy_km, ud$pref)

cells <- rast(bathy_terra)
cells[['ud']] <- 0
cells[ud$cell] <- ud$ud
cells[['pref']] <- 0
cells[['pref']][ud$cell] <- ud$pref
cells[['res']] <- 0
cells[['res']][ud$cell] <- ud$res

ggplot() +
  geom_spatraster(data=cells$ud) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="UD") + 
  # scale_fill_distiller(palette="YlOrRd", direction = 1, na.value = "transparent", name="UD") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) +
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()

ggplot() +
  geom_spatraster(data=cells$pref) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Preference") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) + 
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()

ggplot() +
  geom_spatraster(data=cells$res) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="E[residence] (d)") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) + 
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()





# Actual constant diffusion
fit_sde <- fit_ctmc(walk_data,
                 model_parameters = list(
                   q_r = list(form=~ 1), # residency model
                   q_m = list(form=~bathy_km_grad), # preference model
                   form="sde"
                 ), reals = FALSE, control=list(trace=10)
                 )
# AIC = 333


ud_sde <- get_lim_ud(fit_sde)
beta_q_m <- fit_sde$results$beta$q_m$est
sig <- as.vector(soft_plus(fit_sde$results$beta$q_r$est))
b <- 2*(beta_q_m/sig^2)

# Check if cell size is small enough:
sig^2/max(abs(b*walk_data$q_m$bathy_km_grad))

ud_sde$pi <- exp(-b[1]*walk_data$q_r$bathy_km)# - b[2]*walk_data$q_r$bathy_km^2) #+ beta_q_m[2]*walk_data$q_r$bathy_km2
ud_sde$pi <- ud_sde$pi/sum(ud_sde$pi)

# Limiting distribution from:
# Brillinger, D. R., Preisler, H. K., Ager, A. A., 
# Kie, J. G., & Stewart, B. S. (2002). Employing stochastic differential 
# equations to model wildlife motion. Bulletin of the Brazilian Mathematical 
# Society, 33, 385-408.


cells[['ud_sde']] <- 0
cells[['ud_sde']][ud_sde$cell] <- ud_sde$ud
cells[['pi_sde']] <- 0
cells[["pi_sde"]][ud_sde$cell] <- ud_sde$pi


ggplot() +
  geom_spatraster(data=cells$ud_sde) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="UD") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) +
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()

ggplot() +
  geom_spatraster(data=cells$pi_sde) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Pi") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) + 
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()


######################
### Reconstruct path
######################

aux_time <- seq(head(walk_data$times$timestamp,1), tail(walk_data$times$timestamp,1), "1 hours")
P <- predict_ctmc(fit, walk_data, aux_timestamp=aux_time)

Punif <- P$local_state_prob[P$times$type=="p",]
pred_loc <- rast(cells, nlyr=nrow(Punif))
for(i in 1:nrow(Punif)){
  pred_loc[[i]][walk_data$q_r$cell] <-  Punif[i,]
}
names(pred_loc) <- filter(P$times, type=="p") |> pull(timestamp) |> as.character()

pred_loc[pred_loc<1.0e-4] <- NA

cells[["obs_use"]] <- NA
cells[["obs_use"]][walk_data$q_r$cell] <- colSums(Punif) |> as.vector()

ggplot() +
  geom_spatraster(data=cells$obs_use) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Observed\nuse (h)") + 
  annotation_spatial(ak, fill=gray(0.8), color=1) +
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()


plt_foo <- function(n){
  for(i in 1:n){
  p <- ggplot() +
    geom_spatraster(data=pred_loc[[i]]) + 
    scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Pred. loc") + 
    annotation_spatial(ak, fill=gray(0.8), color=1) +
    scale_y_continuous(breaks=seq(-180,180,1)) +
    scale_x_continuous(breaks=seq(-180,180,5)) +
    theme_bw(base_size=24) + theme(legend.position = "none")
  print(p)
  }
}

saveVideo({
  plt_foo(nlyr(pred_loc))
}, video.name = "cod_pred_loc.mp4", 
interval=0.025, ani.res=72, ani.width=1080, ani.height=600)



