#' @title read the web-of-life data files of Bascompte et al.
#' @references http://www.web-of-life.es/map.php
#' 
#require(devtools)
weboflife = numeric()
for (i in 1:59) {
  networkname = paste('M_PL_', sprintf('%03d', i), sep='')
  weboflife[i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_PL_", sprintf("%03d", i), ".csv", sep=""),  header=T, row.names = 1))
  #devtools::use_data(get(networkname))
}
for (i in 1:30) {
  networkname = paste('M_SD_', sprintf('%03d', i), sep='')
  weboflife[59 + i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_SD_", sprintf("%03d", i), ".csv", sep=""),  header=T, row.names = 1))  
  #devtools::use_data(get(networkname))
}


