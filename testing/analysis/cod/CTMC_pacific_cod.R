library(terra)
library(sf)
# library(TMB)
# library(Matrix)
library(rnaturalearth)
library(ggplot2)
library(marginaleffects)
library(Matrix)
library(ggplot2)
library(tidyterra)
library(ggspatial)
#library(walk)
devtools::load_all("~/research/projects/r_packages/walk")

# setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_9")
setwd("~/research/projects/r_packages/walk/testing/analysis/cod")
source( "Shared_functions/add_legend.R" )

# load data
bathy_terra <- readRDS("ai_bathy_3km.Rds")
likelihood_terra <- readRDS("likelihood_3km.Rds")

# Get land layer
sf_states <- ne_states( "united states of america", return="sf")
sf_states <- sf_states[pmatch("Alas", sf_states$name_en),]
sf_states <- st_transform( sf_states, crs=st_crs("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") )
sf_states <- st_geometry(sf_states)

# change resolution
fact <- 2
bathy_terra <- terra::aggregate( bathy_terra, fact = fact )
likelihood_terra <- as.array( terra::aggregate(likelihood_terra, fact=fact) )

# Eliminate land cells
bathy_terra[bathy_terra<=1] <- NA
bathy_terra$bathy_km <- bathy_terra/1000
bathy_terra$bathy_km_s <- scale(bathy_terra$bathy_km)

# Format likelihood as sparse matrix
L <- t(apply(likelihood_terra, MARGIN=3, FUN=function(mat){as.vector(t(mat))}))
L <- Matrix(L) |> as("CsparseMatrix")
# Format for observed data for {walk}
pdata <- list(L=L)
pdata$times <- data.frame(obs=1:nrow(L), dt=c(0, rep(1,nrow(L)-1)))


# Format habitat for {walk}
walk_data <- make_walk_data(pdata, bathy_terra, grad=c("bathy_km","bathy_km_s"), rast_mask = bathy_terra[[1]])

#####################
# Diffusion-only model
#####################

walk_data$q_r <- mutate(walk_data$q_r,
                        bathy_shallow = ifelse(bathy_km<=1, 1, 0),
                        bathy_deep = ifelse(bathy_km>1, 1, 0)
)

fit1 <- fit_ctmc(walk_data,
                 model_parameters = list(q_r = ~1, q_m = ~1, p = FALSE, delta=NULL, link="log"), 
                 pen_fun = NULL, hessian=TRUE, get_reals=TRUE, 
                 # start=list(beta_q_r = -1.3, beta_q_m=c(0.17,-0.19), logit_p=-3.9),
                 method="nlminb", fit=TRUE, eq_prec = 1.0e-8,
                 control=list(trace=1))
stat1 <- get_lim_ud(fit1)

#####################
# Preference model
#####################
deltaD <- mean(res(bathy_terra)) / 1000
walk_data$q_m$d_bathy <- (walk_data$q_m$bathy_km-walk_data$q_m$from_bathy_km) 



fit2 <- fit_ctmc(walk_data,
                 model_parameters = list(q_r = ~ log(bathy_km), q_m = ~poly(d_bathy, 2, raw=T) , p = FALSE, delta=NULL), 
                 pen_fun = NULL, hessian=TRUE, get_reals=TRUE, 
                 # start=list(beta_q_r=fit1$par),
                 method="nlminb", fit=T, eq_prec = 1.0e-8,
                 control=list(#rel.tol=1.0e-6, 
                   trace=1) )
stat2 <- get_lim_ud(fit2, hpd=c(0.5, 0.9, 0.95))
beta_q_m <- fit2$results$beta$q_m$est
stat2$pref <- beta_q_m[1]*walk_data$q_r$bathy_km + beta_q_m[2]*walk_data$q_r$bathy_km^2

reals <- get_reals(fit2$par, fit2$vcov, walk_data, fit2$model_parameters)

cells <- rast(bathy_terra)

cells[['ud']] <- 0
cells[stat2$cell] <- stat2$ud

cells[['pref']] <- 0
cells[["pref"]][stat2$cell] <- stat2$pref
# cells[["pref"]][cells[["pref"]]<log(1.0e-8)] <- NA

pref_df <- data.frame(
  x = walk_data$q_r$bathy_km,
  x_us = walk_data$q_r$ai_bathy_fill
) %>% mutate(
  y = beta_q_m[1]*x + beta_q_m[2]*x^2
) %>% arrange(x)

plot(y~x_us, type="l", data=pref_df[,])

ggplot() +
  geom_spatraster(data=cells$ud) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="UD") + 
  # scale_fill_distiller(palette="YlOrRd", direction = 1, na.value = "transparent", name="UD") + 
  # annotation_spatial(sf_states, fill=gray(0.8), color=1) + 
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()

ggplot() +
  geom_spatraster(data=cells$pref) + 
  scale_fill_gradientn(colors=sf.colors(10), na.value = "transparent", name="Preference") + 
  annotation_spatial(sf_states, fill=gray(0.8), color=1) + 
  scale_y_continuous(breaks=seq(-180,180,1)) +
  scale_x_continuous(breaks=seq(-180,180,5)) +
  theme_bw()





