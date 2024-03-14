
setwd( R'(C:\Users\James.Thorson\Desktop\Git\walk\data-raw)' )

bathymetry = readRDS( "ai_bathy_3km.Rds" )
data_likelihood = readRDS( "likelihood_3km.Rds" )

pcod_archival_tag = list(
  bathymetry = bathymetry,
  data_likelihood = data_likelihood                       
)

usethis::use_data( pcod_archival_tag, overwrite=TRUE )
