library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(Matrix)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(lubridate)
#library(walk)
devtools::load_all("~/research/projects/r_packages/walk")

# setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_9")
setwd("~/research/projects/r_packages/walk/testing/analysis/cod")
source( "Shared_functions/add_legend.R" )

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
bathy_terra$bathy_km_s <- scale(bathy_terra$bathy_km)
bathy_terra$bathy_t <- terra::ifel(bathy_terra$bathy_km>4, 4, bathy_terra$bathy_km)

ggplot() +
  geom_spatraster(data=log(bathy_terra$bathy_t)) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="UD") + 
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
walk_data <- make_walk_data(pdata, bathy_terra, grad=c("bathy_km","bathy_km_s"), rast_mask = bathy_terra[[1]])

walk_data$q_m$d_bathy <- (walk_data$q_m$bathy_t-walk_data$q_m$from_bathy_t) 

### Get average "d_bathy" for use in residency model
walk_data$q_r$avg_d_bathy <- walk_data$q_m |> group_by(from_cellx) |> 
  summarize(avg_d_bathy = mean(d_bathy)) |> arrange(from_cellx) |> pull(avg_d_bathy)
walk_data$q_r$avg_d_bathy_2 <- walk_data$q_m |> group_by(from_cellx) |> 
  summarize(avg_d_bathy_2 = mean(d_bathy^2)) |> arrange(from_cellx) |> pull(avg_d_bathy_2)


# Diffusion depends on average d_bathy more like traditional CTMC
fit <- fit_ctmc(walk_data,
                 model_parameters = list(
                   q_r = ~ log(bathy_t), # residency model
                   q_m = ~d_bathy + I(d_bathy^2) # preference model
                   ), reals = TRUE, # get est. of residency expectation and multinomial movement probability
                 control=list(trace=1))
# AIC = 298

# Actual constant diffusion
fit2 <- fit_ctmc(walk_data,
                 model_parameters = list(
                   q_r = ~ 1, # residency model
                   q_m = ~d_bathy + I(d_bathy^2), # preference model
                   link="log",
                   norm=FALSE
                 ), reals = TRUE, control=list(trace=10), fit=FALSE)
# AIC = 333

ud <- get_lim_ud(fit)
beta_q_m <- fit2$results$beta$q_m$est
stat2$pref <- beta_q_m[1]*walk_data$q_r$bathy_km_trunc + beta_q_m[2]*walk_data$q_r$bathy_km_trunc^2



cells[['ud']] <- 0
cells[['ud']][stat2$cell] <- stat2$ud

cells[['pref']] <- 0
cells[["pref"]][stat2$cell] <- stat2$pref

cells[['res']] <- 0
cells[["res"]][stat2$cell] <- fit2$results$real$residency$real
# cells[["pref"]][cells[["pref"]]<log(1.0e-8)] <- NA

pref_df <- data.frame(
  x = walk_data$q_r$bathy_km_trunc,
  x_us = walk_data$q_r$ai_bathy_fill
) %>% mutate(
  y = beta_q_m[1]*x + beta_q_m[2]*x^2
) %>% arrange(x)

plot(y~x_us, type="l", data=pref_df[,])

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

# ggplot() +
#   geom_spatraster(data=cells$res) + 
#   scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Avg. residency") + 
#   annotation_spatial(ak, fill=gray(0.8), color=1) + 
#   scale_y_continuous(breaks=seq(-180,180,1)) +
#   scale_x_continuous(breaks=seq(-180,180,5)) +
#   theme_bw()


######################
### Reconstruct path
######################

aux_time <- seq(head(walk_data$times$timestamp,1), tail(walk_data$times$timestamp,1), "1 hours")
P <- predict_ctmc(fit2, walk_data, aux_timestamp=aux_time, debug=0)

Punif <- P$local_state_prob[P$times$type=="p",]

pred_loc <- rast(cells, nlyr=nrow(Punif))
for(i in 1:nrow(Punif)){
  pred_loc[[i]][walk_data$q_r$cell] <-  Punif[i,]
}
names(pred_loc) <- filter(P$times, type=="p") |> pull(timestamp) |> as.character()

pred_loc[pred_loc<1.0e-4] <- NA

cells[["obs_use"]] <- NA
# cells[["obs_use"]][apply(P$local_state_prob, 1, which.max)] <- 1
cells[["obs_use"]][walk_data$q_r$cell] <- colSums(Punif) |> as.vector()
cells[["obs_use"]][walk_data$q_r$cell] <- Punif[1400,]

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
  plt_foo(300)
}, video.name = "cod_pred_loc.mp4", 
interval=0.025, ani.res=72, ani.width=1080, ani.height=600)



