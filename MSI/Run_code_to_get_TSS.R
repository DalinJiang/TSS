#******************************************************************************************************************************************************************
#
# This R code is to run the TSS estimation function to get TSS estimates using the method proposed by Jiang et al. (2023)
#
#
# GENERAL INPUT/OUTPUT OF THE TSS ESTIMSTION METHOD
#---------------------------------------------------
# (1) inputs are:
# Remote sensing reflectance (Rrs, sr-1) at MSI bands
# i.e., 443,490,560,665,705,740,783,865 nm
#
#-------------
# (2) outputs are:
# 	<1> estimated absorption coefficient (a, m-1) at reference band
#	<2> particulate backscatering coefficient (bbp, m-1) at reference band
# 	<3> the reference band (band, nm) used for TSS estimation
#	<4> estimated TSS (g/m3)
# 
#-------------
# for more informaiton, please refer to the publication:
# Jiang et al. (2023). Estimating the concentration of total suspended solids in inland and coastal waters from Sentinel-2 MSI: a semi-analytical approach. 
#						ISPRS International Journal of Photogrammetry and Remote sensing.
#
#
#
#
#
# TWO EXAMPLES OF DATA PROCESSING
#---------------------------------------------------
# (1) For in situ measured Rrs data (case 1 as below):
# 		<1> save your Rrs in a csv file, and read the csv data into a dataframe or a 2D matrix
#		<2> apply the TSS estimation function ('Estimate_TSS_Jiang_MSI') to the 2D matrix
#		<3> export the estimated result to a csv file
#
# (2) For satellite images (after atmospheric correction, case 2 as bellow):
#		<1> read the satellite data into a 3D matrix, e.g., nrow*ncol*nbands
#		<2> apply the TSS estimation function ('Estimate_TSS_Jiang_MSI') to the 3D matrix
#		<3> export the estimated result to a GeoTiff file
#
#
#
#
# Version = V1.0
# Dalin JIANG, University of Stirling, UK
# 21 August, 2023
#
#
#******************************************************************************************************************************************************************



#---------------------------Case 1: estimate TSS from in situ measured Rrs spectra-------------------
#
# change the path below to the folder where your code is
setwd('/Volumes/R/TSS_MSI/')

# source the TSS estimation function
source('TSS_function_Jiang_MSI.R')

# read the input data
in_file <- "example_insitu_Rrs_input.csv"
in_dt <- read.csv(in_file,header=TRUE)

MSI_wave <- c(443,490,560,665,705,740,783,865)
MSI_band <- paste("Rrs",MSI_wave,sep="")

# do the TSS estimation
val_Rrs <- in_dt[,MSI_band]
est_tss <- data.frame(t(apply(val_Rrs,1,Estimate_TSS_Jiang_MSI)))
final_result <- cbind(in_dt,est_tss)

# export the results to a csv file
out_file <- "example_insitu_Rrs_output.csv"
write.csv(final_result, out_file, row.names=FALSE)


rm(list=ls())
gc()

















#----------------------------- Case 2: estimate TSS from MSI image------------------------------------------
#
# Note: the example MSI image below is atmospherically corrected using C2RCC in SNAP
# The 'ncdf4','raster','abind' packages are required to open nc file
# to install the package, use 'install.packages("raster")' 
#
# output file include 4 bands: 
#		<1> band1: a
#		<2> band2: bbp
#		<3> band3: reference band
#		<4> band4: TSS
#
#

library(raster)    
library(ncdf4)
library(abind)

# change the path below to the folder where your code is
setwd('/Volumes/R/TSS_MSI/')

# source the TSS estimation function
source('TSS_function_Jiang_MSI.R')


# function to open file and do the calculation
estimate_TSS <- function(nc_img){
	nc_dt <- nc_open(nc_img)
	all_bd <- c("rrs_B1","rrs_B2","rrs_B3","rrs_B4","rrs_B5","rrs_B6","rrs_B7","rrs_B8A")
	lon <- ncvar_get(nc_dt,"lon")
	lat <- ncvar_get(nc_dt,"lat")
	img_rg <- c(min(lon),max(lon),min(lat),max(lat))
	img_crs <- "EPSG:4326"   # change this according to your data
			
	ext_one <- function(nc_data,bd_name){
		one_rrs <- ncvar_get(nc_data,bd_name)
		return(one_rrs) 
	}
	
	all_rrs <- sapply(all_bd,ext_one,nc_data=nc_dt,simplify=FALSE)  
	rrs_dt <- abind(all_rrs,along=3)  
	
	
	xx <- rrs_dt	
	tmp_result <- apply(xx, c(1,2),Estimate_TSS_Jiang_MSI)
	
	est_a <- as.matrix(tmp_result[1,,])
	est_bbp <- as.matrix(tmp_result[2,,])
	ref_band <- as.matrix(tmp_result[3,,])
	tss_jiang <- as.matrix(tmp_result[4,,])
	
	a_img <- raster(t(est_a),xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	bbp_img <- raster(t(est_bbp),xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	refband_img <- raster(t(ref_band),xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	tss_img <- raster(t(tss_jiang),xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	
	out_img <- stack(a_img,bbp_img,refband_img,tss_img)
	
	nc_close(nc_dt)
	return(out_img)
}


#-------------------------------------
input_img <- "example_MSI_image_input.nc"
tss_result <- estimate_TSS(input_img)

output_img <- "example_MSI_image_output.tif"
writeRaster(tss_result,output_img,format="GTiff")   

rm(list=ls())
gc()
