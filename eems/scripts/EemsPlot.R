if (file.exists("./rEEMSplots")) {
  install.packages("rEEMSplots",repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

coord=read.table("oneperrad_trachylepis.coord")
mcmcpath =c( "oneperrad_trachylepis-EEMS-nDemes200-1","oneperrad_trachylepis-EEMS-nDemes200-2","oneperrad_trachylepis-EEMS-nDemes200-3","oneperrad_trachylepis-EEMS-nDemes200-4","oneperrad_trachylepis-EEMS-nDemes200-5")
plotpath = "plots_2019.10.29wmapandpoints"

eems.plots(mcmcpath, plotpath, longlat = TRUE,projection.in='+proj=longlat +datum=WGS84',projection.out='+proj=longlat +datum=WGS84',add.map=TRUE,m.plot.xy={ points(coord, col = "black") },q.plot.xy={ points(coord, col = "black") })

