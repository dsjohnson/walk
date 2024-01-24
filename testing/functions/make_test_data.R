library(ctmm)
# library(walk)
library(terra)
devtools::load_all(".") #load walk functions




###

data(pelican)
data <- ctmm::tbind(pelican)
data$COV.x.x <- ifelse(data$COV.x.x==Inf, 20^2, data$COV.x.x)
data$COV.y.y <- ifelse(data$COV.y.y==Inf, 20^2, data$COV.y.y)
xy <- as.matrix(data[,c("x","y")])
sa <- terra::rast(
  xmin=min(xy[,1])-20000, xmax=max(xy[,1])+20000,
  ymin=min(xy[,2])-20000, ymax=max(xy[,2])+20000,
  crs=projection(data), res=10000
)
telem_pts <- terra::vect(cbind(data$longitude, data$latitude), crs="epsg:4326") |>
  terra::project(sa)
t_buf <- terra::buffer(telem_pts, 300000)
ras <- terra::rasterize(t_buf, sa)
# sp_data <- telem_to_ras(data, ras)
