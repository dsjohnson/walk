#' Simulate a coastline (mainland + islands) placed relative to where the animal lives.
#'
#' The old generator put the mainland against the western edge of the raster and
#' scattered islands anywhere, while the 'd2c' term keeps the animal near the origin
#' (UD sd = sqrt(sca / (2 |beta_d2c|)) ~ 32 cells). The coast was therefore usually
#' 40-130 cells from the track, so observations could only reach land in the few
#' replicates where the coast happened to cut into the core. Here the mainland is a
#' half-plane in a random direction whose (wiggly) coast passes a chosen distance from
#' the UD center, and islands are confined to a disc around the center.
#'
#' @param template SpatRaster defining the grid (e.g. spatialCovs$cov1).
#' @param udSD     Scale of the animal's range (cells); distances below are multiples of it.
#' @param coastDist Range (in udSD units) for the distance from the UD center to the mean coastline.
#' @param wiggle   Range (in udSD units) for the amplitude of the coastline wiggle.
#' @param islandRadius Islands are only allowed within this distance (udSD units) of the center.
#' @param islandProb  Probability that a replicate has islands at all.
#' @param islandFrac  Range for the fraction of the island disc covered by islands.
#' @param clear    Minimum water distance (cells) kept around the center (start/UD mode stays at sea).
#' @param center   UD center (default origin, where d2c = 0).
#' @return SpatRaster named "coast_barrier": 1 = water (allowed), 0 = land (restricted).
simLand <- function(template, udSD, coastDist = c(0.3, 1), wiggle = c(0.1, 0.4),
                    islandRadius = 1.5, islandProb = 0.5, islandFrac = c(0.02, 0.08),
                    clear = 5, center = c(0, 0), maxTries = 20) {
  r <- terra::rast(template[[1]]); terra::crs(r) <- ""
  X <- terra::init(r, "x") - center[1]; Y <- terra::init(r, "y") - center[2]
  R <- sqrt(X^2 + Y^2)
  smooth_noise <- function() {                      # smooth field scaled to [-1, 1]
    n <- r; terra::values(n) <- stats::runif(terra::ncell(n))
    ws <- sample(seq(11, 21, by = 2), 1); w <- matrix(1, ws, ws)
    for (i in 1:4) n <- terra::focal(n, w = w, fun = mean, na.policy = "omit", expand = TRUE)
    mm <- terra::minmax(n); 2 * (n - mm[1]) / (mm[2] - mm[1]) - 1
  }
  for (k in seq_len(maxTries)) {
    ## mainland: half-plane in a random direction with a noisy coastline
    th <- stats::runif(1, 0, 2 * pi)
    d0 <- stats::runif(1, coastDist[1], coastDist[2]) * udSD
    amp <- stats::runif(1, wiggle[1], wiggle[2]) * udSD
    proj <- X * cos(th) + Y * sin(th)
    land <- proj > (d0 + amp * smooth_noise())
    ## islands: high values of a second smooth field, only inside the island disc
    if (stats::runif(1) < islandProb) {
      n2 <- smooth_noise()
      inDisc <- R < islandRadius * udSD
      v <- terra::values(n2)[terra::values(inDisc) == 1]
      q <- stats::quantile(v, 1 - stats::runif(1, islandFrac[1], islandFrac[2]), na.rm = TRUE)
      land <- land | ((n2 > q) & inDisc)
    }
    land <- land & (R > clear)                         # keep the center at sea
    wm <- terra::ifel(land, 0, 1); names(wm) <- "coast_barrier"
    fl <- mean(terra::values(wm) == 0)
    if (fl > 0.02 && fl < 0.7) return(wm)              # must contain both land and water
  }
  stop("simLand: could not generate a valid land mask")
}
