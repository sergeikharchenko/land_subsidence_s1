library(rgdal)
library(raster)

### SETTING WORK DIRECTORY
pf_folder <- "D:/Baltics/"
dir.create(pf_folder)
setwd(pf_folder)

##### 
### CREATE CONVERTING SCRIPT (from *.DIM to *.TIFF)
writeLines(text = '
           <graph id="Graph">
  <version>1.0</version>
  <node id="Read">
    <operator>Read</operator>
    <sources/>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>disp_tc.dim</file>
    </parameters>
  </node>
  <node id="Write">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="Read"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>disp_tc.tif</file>
      <formatName>GeoTIFF</formatName>
    </parameters>
  </node>
</graph>
           ', con = "write.xml")

#####

### CHOOSE SPATAIL AND TEMPORARY BORDERS
bbox_sen1 <- c(20, 54.7, 20.5, 55.2)
start_day <- "2022-01-01"
end_day <- "2022-05-01"
sequent <- T


earthdatalogin <- "eal" # your earthdata login
earthdatapw <- "edpw" # your earthdata password


link_api <- paste0("https://api.daac.asf.alaska.edu/services/search/param?platform=Sentinel-1&beamMode=IW&processingLevel=SLC&bbox=",paste(bbox_sen1, collapse = ","),"&start=",start_day,"T00:00:00Z&end=",end_day,"T23:59:59Z&output=CSV")

query1 <- read.csv(link_api)
write.csv(query1, "query1.csv")


(uniq_path_frame <- cbind(query1$Path.Number, query1$Frame.Number)[!duplicated(cbind(query1$Path.Number, query1$Frame.Number)),])
w <- 3 # NULL
query1 <- query1[query1$Path.Number == uniq_path_frame[w,1] & query1$Frame.Number == uniq_path_frame[w,2],]

#if (is.null(w)) w <- which.min(as.matrix(dist(rbind(c(mean(bbox_sen1[c(1,3)]), mean(bbox_sen1[c(2,4)])), query1[,c("Center.Lon", "Center.Lat")])))[-1,1])

#if (is.null(dim(uniq_path_frame))) {
#  query1 <- query1[query1$Path.Number == uniq_path_frame[1] & query1$Frame.Number == uniq_path_frame[2],]
#} else {
#  query1 <- query1[query1$Path.Number == uniq_path_frame[w,1] & query1$Frame.Number == uniq_path_frame[w,2],]
#}

bbox_image <- apply(matrix(query1[1,16:23], ncol = 2, byrow = T)[,2:1], 2, function(x) as.numeric(x))
pol <- as(extent(bbox_sen1[c(1,3,2,4)]), "SpatialPolygons")
crs(pol) <- "+proj=longlat +datum=WGS84"
pol$ID <- 1

WKT_pol <- paste(c("POLYGON((", (paste(apply(pol@polygons[[1]]@Polygons[[1]]@coords, 1, function(x) paste0(x[1]," ", x[2])), collapse = ",")),"))"), collapse = "")

writeOGR(pol, "polygon.shp", "polygon.shp", "ESRI Shapefile")
plot(SpatialPoints(coords = bbox_image), pch = 19, cex = 3, col = rainbow(4))
text(bbox_image[,1], bbox_image[,2], 1:4)
points(SpatialPoints(expand.grid(bbox_sen1[c(T,F)], bbox_sen1[c(F,T)])), pch = 19, col = "red")

(list_sen1 <- query1$URL)

for (i in 1:length(list_sen1)) {
  download.file(paste0("http://",earthdatalogin,":",earthdatapw,"@", substr(list_sen1[i], 9, nchar(list_sen1[i]))), paste0(query1$Granule.Name[i],".zip"))
}
#

##### MODIS snow cover
modis_request <- paste0("https://",earthdatalogin,":",earthdatapw,"@n5eil02u.ecs.nsidc.org/egi/request?short_name=MOD10A2&version=6&time=",start_day,"T00:00:01,",end_day,"T23:59:59&bounding_box=",paste(bbox_sen1, collapse = ","),"&agent=NO&page_size=100")
dir.create("modis")
download.file(url = modis_request, destfile = 'modis/modis.zip', mode='wb')
unzip(zipfile = "modis/modis.zip", exdir = "modis")

modis_lf <- list.files(path = "modis", pattern = "*.hdf$", recursive = T, full.names = T)

modis_value <- cbind(rep(NA, length(modis_lf)), rep(NA, length(modis_lf)))
for(i in 2:length(modis_lf)) {
  # Extraction of metadata via `GDALinfo`  
  # i <- 2
  filename <- modis_lf[i]
  invisible(gdalinfo <- GDALinfo(filename, returnScaleOffset = FALSE))
  metadata <- attr(gdalinfo, "subdsmdata")
  
  # gdalinfo[grep(pattern = "DATE", x = gdalinfo)]
  
  sds <- metadata[grep("Snow", metadata)[1]]
  sds <- sapply(strsplit(sds, "="), "[[", 2)
  
  # Raster import via `readGDAL`   
  invisible(sds.rg <- readGDAL(sds, as.is = T))
  modis_r <- raster(sds.rg)
  plot(modis_r)
  pol <- spTransform(pol, crs(modis_r))
  modis_value[i,1] <- cellStats(crop(modis_r, pol), mean)
  modis_value[i,2] <- cellStats(crop(modis_r, pol), sd)
  
  #plot(pol, add = T)
  #print(paste0("Snowcover is ",(cellStats(crop(modis_r, pol), mean))))
  #plot(modis_r, main = substr(modis_lf[i], 27, 33))
  #Sys.sleep(3)
}

plot(modis_value[modis_value[,2] != 0 ,1])


if(sequent) {
  for (i in 1:(nrow(query1)-1)) {
    s1files <- list.files(pattern = "*.zip")
    (s1files <- s1files[i:(i+1)])
    
    PselectedPolarisations <- "VV"
    PdemName <- "ACE30" # "Copernicus 30m Global DEM", "Copernicus 90m Global DEM", "ACE30", "GETASSE30"
    Psubswath <- "IW1"
    dir.create(Psubswath)
    PFFTSize <- 64 # 32,62,128,256
    subset <- F
    mlook <- T
    snaphu_folder <- "se1"
    initMethod <- "MCF"
    colOverlap <- 300
    rowOverlap <- 300
    nproc <- 60
    numberOfTileCols <- 15
    numberOfTileRows <- 15
    SCM <- "DEFO"
    
    system(paste0('gpt TOPSAR-Split -t "t1_split.dim" -Ssource=',s1files[1],' -PselectedPolarisations="',PselectedPolarisations,'" -Psubswath="',Psubswath,'" -PfirstBurstIndex=1'), show.output.on.console = F)
    system(paste0('gpt TOPSAR-Split -t "t2_split.dim" -Ssource=',s1files[2],' -PselectedPolarisations="',PselectedPolarisations,'" -Psubswath="',Psubswath,'" -PfirstBurstIndex=1'), show.output.on.console = F)
    
    system(paste0('gpt Apply-Orbit-File -t "t1_split_of.dim" -Ssource="t1_split.dim" '), show.output.on.console = F)
    system(paste0('gpt Apply-Orbit-File -t "t2_split_of.dim" -Ssource="t2_split.dim" '), show.output.on.console = F)
    
    # update Java!!
    system(paste0('gpt Back-Geocoding -t "t1_split_of_bgc.dim" -PdemName="',PdemName,'" -Ssource1="t1_split_of.dim" -Ssource2="t2_split_of.dim" '), show.output.on.console = T)
    
    system(paste0('gpt Enhanced-Spectral-Diversity -t "t1_split_of_bgc_esd.dim" -Ssource="t1_split_of_bgc.dim" '), show.output.on.console = F)
    
    system(paste0('gpt Interferogram -t "t1_split_of_bgc_esd_ifg.dim" -SsourceProduct="t1_split_of_bgc_esd.dim" -PdemName="',PdemName,'" -PcohWinRg=20 '), show.output.on.console = F)
    
    # system(paste0('gpt ifg.xml'), show.output.on.console = F)
    
    system(paste0('gpt TOPSAR-Deburst -t "t1_split_of_bgc_esd_ifg_deb.dim" -Ssource="t1_split_of_bgc_esd_ifg.dim" -PselectedPolarisations="',PselectedPolarisations,'" '), show.output.on.console = F)
    
    if (mlook) {
      system(paste0('gpt TopoPhaseRemoval -t "t1_split_of_bgc_esd_ifg_deb_tpr.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb.dim" -PdemName="',PdemName,'" -PoutputTopoPhaseBand="true" '), show.output.on.console = T)
      system(paste0('gpt Multilook -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -Ssource="t1_split_of_bgc_esd_ifg_deb_tpr.dim" -PnRgLooks=8 '), show.output.on.console = F)
    } else {
      system(paste0('gpt TopoPhaseRemoval -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb.dim" -PdemName="',PdemName,'" -PoutputTopoPhaseBand="true" '), show.output.on.console = F)
    }
    
    system(paste0('gpt GoldsteinPhaseFiltering -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -PuseCoherenceMask="true" -PcoherenceThreshold=0.3 -PFFTSizeString="',PFFTSize,'" '), show.output.on.console = F)
    
    system(paste0('gpt SnaphuExport -PtargetFolder="',snaphu_folder,'" -Ssource="t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" -PinitMethod="',initMethod,'" -PcolOverlap=',colOverlap,' -ProwOverlap=',rowOverlap,' -PnumberOfProcessors=',nproc,' -PnumberOfTileCols=',numberOfTileCols,' -PnumberOfTileRows=',numberOfTileRows,' -PstatCostMode="',SCM,'" -PtileCostThreshold=500 '), show.output.on.console = F)
    
    (magic_str <- readLines(list.files(snaphu_folder, pattern = "snaphu.conf", recursive = T, full.names = T))[7])
    wd1 <- getwd()
    setwd(file.path(pf_folder,snaphu_folder,"/t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp/"))
    system(paste0("C:/Users/Professional/.snap/auxdata/snaphu-v2.0.4_win64/bin/snaphu.exe",substr(magic_str, start = 15, stop = nchar(magic_str))))
    setwd(wd1)
    
    unwp <- strsplit(intersect(list.files(pattern = "UnwPha",recursive = T), list.files(pattern = "*.img",recursive = T)), split = "/")[[1]]
    unwp <- unwp[length(unwp)]
    system(paste(c('gpt SnaphuImport -t "si.dim" "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" "',pf_folder,snaphu_folder,'/t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp/',unwp), collapse = ""), show.output.on.console = F)
    
    system(paste0('gpt PhaseToDisplacement -t "disp.dim" -Ssource="si.dim" '), show.output.on.console = F)
    
    if (subset) {
      system(paste0('gpt Subset -t "disp.dim" -SsourceProduct="disp.dim" -PgeoRegion="',WKT_pol,'"'), show.output.on.console = T)
    }
    
    system(paste0('gpt Terrain-Correction -t "disp_tc',i,'.dim" -PdemName="',PdemName,'" -Ssource="disp.dim"  '), show.output.on.console = F)
    system(paste0('gpt Terrain-Correction -t "disp_tc.dim" -PdemName="',PdemName,'" -Ssource="disp.dim"  '))
    
    
    
    system(paste0('gpt write.xml'))
    
    library(raster)
    disp_tc <- raster("disp_tc.tif")
    # ifg <- raster("ifg.tif", band = 3)
    disp_tc[disp_tc > -0.01 & disp_tc < 0.01] <- NA
    writeRaster(disp_tc, paste0(Psubswath,"/disp_tc_last",i,".tif"), "GTiff", overwrite = T)
    unlink(c("disp_tc.tif", "disp_tc.dim", "disp.dim", "si.dim", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim", "t1_split_of_bgc_esd_ifg_deb_tpr.dim", "t1_split_of_bgc_esd_ifg_deb.dim", "t1_split_of_bgc_esd_ifg.dim", "t1_split_of_bgc_esd.dim", "t1_split_of_bgc.dim", "t2_split_of.dim", "t1_split_of.dim", "t1_split.dim", "t2_split.dim"))
    unlink(c("disp_tc.data", "disp.data", "si.data", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.data", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.data", "t1_split_of_bgc_esd_ifg_deb_tpr.data", "t1_split_of_bgc_esd_ifg_deb.data", "t1_split_of_bgc_esd_ifg.data", "t1_split_of_bgc_esd.data", "t1_split_of_bgc.data", "t2_split_of.data", "t1_split_of.data", "t1_split.data", "t2_split.data", snaphu_folder), recursive = T)
    plot(disp_tc)
  }
}

if(!sequent) {
  #### ALL COMBINATIONS OF IMAGES
  all_comb <- expand.grid(1:nrow(query1), 1:nrow(query1))
  all_comb <- all_comb[apply(all_comb, 1, function(x) x[1] != x[2] & x[1] < x[2]),]
  all_comb <- all_comb[order(all_comb[,1]),]
  
  for (i in 2:nrow(all_comb)) {
    s1files <- list.files(pattern = "*.zip")
    (s1files <- s1files[c(all_comb[i,1], all_comb[i,2])])
    
    PselectedPolarisations <- "VV"
    PdemName <- "ACE30" # "Copernicus 30m Global DEM", "Copernicus 90m Global DEM", "ACE30", "GETASSE30"
    Psubswath <- "IW2"
    dir.create(Psubswath)
    PFFTSize <- 64 # 32,62,128,256
    subset <- F
    mlook <- T
    snaphu_folder <- "se1"
    initMethod <- "MCF"
    colOverlap <- 300
    rowOverlap <- 300
    nproc <- 60
    numberOfTileCols <- 15
    numberOfTileRows <- 15
    SCM <- "DEFO"
    
    system(paste0('gpt TOPSAR-Split -t "t1_split.dim" -Ssource=',s1files[1],' -PselectedPolarisations="',PselectedPolarisations,'" -Psubswath="',Psubswath,'" -PfirstBurstIndex=1'), show.output.on.console = F)
    system(paste0('gpt TOPSAR-Split -t "t2_split.dim" -Ssource=',s1files[2],' -PselectedPolarisations="',PselectedPolarisations,'" -Psubswath="',Psubswath,'" -PfirstBurstIndex=1'), show.output.on.console = F)
    
    system(paste0('gpt Apply-Orbit-File -t "t1_split_of.dim" -Ssource="t1_split.dim" '), show.output.on.console = F)
    system(paste0('gpt Apply-Orbit-File -t "t2_split_of.dim" -Ssource="t2_split.dim" '), show.output.on.console = F)
    
    # update Java!!
    system(paste0('gpt Back-Geocoding -t "t1_split_of_bgc.dim" -PdemName="',PdemName,'" -Ssource1="t1_split_of.dim" -Ssource2="t2_split_of.dim" '), show.output.on.console = T)
    
    system(paste0('gpt Enhanced-Spectral-Diversity -t "t1_split_of_bgc_esd.dim" -Ssource="t1_split_of_bgc.dim" '), show.output.on.console = F)
    
    system(paste0('gpt Interferogram -t "t1_split_of_bgc_esd_ifg.dim" -SsourceProduct="t1_split_of_bgc_esd.dim" -PdemName="',PdemName,'" -PcohWinRg=20 '), show.output.on.console = F)
    
    # system(paste0('gpt ifg.xml'), show.output.on.console = F)
    
    system(paste0('gpt TOPSAR-Deburst -t "t1_split_of_bgc_esd_ifg_deb.dim" -Ssource="t1_split_of_bgc_esd_ifg.dim" -PselectedPolarisations="',PselectedPolarisations,'" '), show.output.on.console = F)
    
    if (mlook) {
      system(paste0('gpt TopoPhaseRemoval -t "t1_split_of_bgc_esd_ifg_deb_tpr.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb.dim" -PdemName="',PdemName,'" -PoutputTopoPhaseBand="true" '), show.output.on.console = T)
      system(paste0('gpt Multilook -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -Ssource="t1_split_of_bgc_esd_ifg_deb_tpr.dim" -PnRgLooks=8 '), show.output.on.console = F)
    } else {
      system(paste0('gpt TopoPhaseRemoval -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb.dim" -PdemName="',PdemName,'" -PoutputTopoPhaseBand="true" '), show.output.on.console = F)
    }
    
    system(paste0('gpt GoldsteinPhaseFiltering -t "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" -SsourceProduct="t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim" -PuseCoherenceMask="true" -PcoherenceThreshold=0.3 -PFFTSizeString="',PFFTSize,'" '), show.output.on.console = F)
    
    system(paste0('gpt SnaphuExport -PtargetFolder="',snaphu_folder,'" -Ssource="t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" -PinitMethod="',initMethod,'" -PcolOverlap=',colOverlap,' -ProwOverlap=',rowOverlap,' -PnumberOfProcessors=',nproc,' -PnumberOfTileCols=',numberOfTileCols,' -PnumberOfTileRows=',numberOfTileRows,' -PstatCostMode="',SCM,'" -PtileCostThreshold=500 '), show.output.on.console = F)
    
    (magic_str <- readLines(list.files(snaphu_folder, pattern = "snaphu.conf", recursive = T, full.names = T))[7])
    wd1 <- getwd()
    setwd(file.path(pf_folder,snaphu_folder,"/t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp/"))
    system(paste0("C:/Users/Professional/.snap/auxdata/snaphu-v2.0.4_win64/bin/snaphu.exe",substr(magic_str, start = 15, stop = nchar(magic_str))))
    setwd(wd1)
    
    unwp <- strsplit(intersect(list.files(pattern = "UnwPha",recursive = T), list.files(pattern = "*.img",recursive = T)), split = "/")[[1]]
    unwp <- unwp[length(unwp)]
    system(paste(c('gpt SnaphuImport -t "si.dim" "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim" "',pf_folder,snaphu_folder,'/t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp/',unwp), collapse = ""), show.output.on.console = F)
    
    system(paste0('gpt PhaseToDisplacement -t "disp.dim" -Ssource="si.dim" '), show.output.on.console = F)
    
    if (subset) {
      system(paste0('gpt Subset -t "disp.dim" -SsourceProduct="disp.dim" -PgeoRegion="',WKT_pol,'"'), show.output.on.console = T)
    }
    
    system(paste0('gpt Terrain-Correction -t "disp_tc',i,'.dim" -PdemName="',PdemName,'" -Ssource="disp.dim"  '), show.output.on.console = F)
    system(paste0('gpt Terrain-Correction -t "disp_tc.dim" -PdemName="',PdemName,'" -Ssource="disp.dim"  '))
    
    
    
    system(paste0('gpt write.xml'))
    
    library(raster)
    disp_tc <- raster("disp_tc.tif")
    # ifg <- raster("ifg.tif", band = 3)
    disp_tc[disp_tc > -0.01 & disp_tc < 0.01] <- NA
    writeRaster(disp_tc, paste0(Psubswath,"/disp_tc_last",i,"_",paste(unlist(lapply(strsplit(x = s1files, split = "SDV_|T"), function(x) x[2])), collapse = "_"),".tif"), "GTiff", overwrite = T)
    unlink(c("disp_tc.tif", "disp_tc.dim", "disp.dim", "si.dim", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.dim", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.dim", "t1_split_of_bgc_esd_ifg_deb_tpr.dim", "t1_split_of_bgc_esd_ifg_deb.dim", "t1_split_of_bgc_esd_ifg.dim", "t1_split_of_bgc_esd.dim", "t1_split_of_bgc.dim", "t2_split_of.dim", "t1_split_of.dim", "t1_split.dim", "t2_split.dim"))
    unlink(c("disp_tc.data", "disp.data", "si.data", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook_fp.data", "t1_split_of_bgc_esd_ifg_deb_tpr_mlook.data", "t1_split_of_bgc_esd_ifg_deb_tpr.data", "t1_split_of_bgc_esd_ifg_deb.data", "t1_split_of_bgc_esd_ifg.data", "t1_split_of_bgc_esd.data", "t1_split_of_bgc.data", "t2_split_of.data", "t1_split_of.data", "t1_split.data", "t2_split.data", snaphu_folder), recursive = T)
    plot(disp_tc)
  }
}