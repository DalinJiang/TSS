#******************************************************************************************************************************************************************
#
# This R code is to run the TSS estimation function to get TSS estimates using the method proposed by Jiang et al. (2021)
#
#
# GENERAL INPUT/OUTPUT OF THE TSS ESTIMSTION METHOD
#---------------------------------------------------
# (1) inputs are:
# Remote sensing reflectance (Rrs, sr-1) at MERIS or OLCI bands
# i.e., 443,490,560,620,665,754,865 nm
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
# Jiang et al.(2021), Remotely estimating total suspended solids concentration in clear to extremely turbid waters using a novel semi-analytical method. 
#					  Remote Sensing of Environment.
#
#
#
#
#
# TWO EXAMPLES OF DATA PROCESSING
#---------------------------------------------------
# (1) For in situ measured Rrs data (case 1 as below):
# 		<1> save your Rrs in a csv file, and read the csv data into a dataframe or a 2D matrix
#		<2> apply the TSS estimation function ('Estimate_TSS_Jiang_MERIS_OLCI') to the 2D matrix
#		<3> export the estimated result to a csv file
#
# (2) For satellite images (after atmospheric correction, case 2 as bellow):
#		<1> read the satellite data into a 3D matrix, e.g., nrow*ncol*nbands
#		<2> apply the TSS estimation function ('Estimate_TSS_Jiang_MERIS_OLCI') to the 3D matrix
#		<3> export the estimated result to a GeoTiff file
#
#
#
#
# Version = V1.0
# Dalin JIANG, University of Tsukuba, Japan
# 09 November, 2021
#
#
#******************************************************************************************************************************************************************



#---------------------------Case 1: estimate TSS from in situ measured Rrs spectra-------------------
#
# change the path below to the folder where your code is
setwd('/Volumes/R/TSS_OLCI/')

# source the TSS estimation function
source('TSS_function_Jiang_MERIS_OLCI.R')

# read the input data
in_file <- "example_insitu_Rrs_input.csv"
in_dt <- read.csv(in_file,header=TRUE)

OLCI_wave <- c(443,490,560,620,665,754,865)
OLCI_band <- paste("Rrs",OLCI_wave,sep="")

# do the TSS estimation
val_Rrs <- in_dt[,OLCI_band]
est_tss <- data.frame(t(apply(val_Rrs,1,Estimate_TSS_Jiang_MERIS_OLCI)))
final_result <- cbind(in_dt,est_tss)

# export the results to a csv file
out_file <- "example_insitu_Rrs_output.csv"
write.csv(final_result, out_file, row.names=FALSE)


rm(list=ls())
gc()

















#----------------------------- Case 2: estimate TSS from MEIRS/OLCI image------------------------------------------
#
# Note: the example OLCI image below is atmospherically corrected using C2RCC in SNAP
# The 'raster' package is required to open GeoTiff file
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

# change the path below to the folder where your code is
setwd('/Volumes/R/TSS_OLCI/')

# source the TSS estimation function
source('TSS_function_Jiang_MERIS_OLCI.R')

# function to open file and do the calculation
estimate_TSS <- function(input_img){
	img_dt <- stack(input_img)           # open all the bands of the image file
	img_crs <- crs(img_dt)               # get the projection of the image
	img_rg <- extent(img_dt)             # get the spatial extent of the image
	
	Rrs_dt <- img_dt[[c(4,5,7,8,9,13,15)]]              # Rrs443 to Rrs865, from C2RCC for OLCI; if for MERIS or other atmospheric correction method, change this line here
	xx <- as.array(Rrs_dt) 
	
	tmp_result <- apply(xx, c(1,2),Estimate_TSS_Jiang_MERIS_OLCI)
	
	est_a <- as.matrix(tmp_result[1,,])
	est_bbp <- as.matrix(tmp_result[2,,])
	ref_band <- as.matrix(tmp_result[3,,])
	tss_jiang <- as.matrix(tmp_result[4,,])
	
	a_img <- raster(est_a,xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	bbp_img <- raster(est_bbp,xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	refband_img <- raster(ref_band,xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	tss_img <- raster(tss_jiang,xmn=img_rg[1],xmx=img_rg[2],ymn=img_rg[3],ymx=img_rg[4],crs=img_crs) 
	
	out_img <- stack(a_img,bbp_img,refband_img,tss_img)
	
	return(out_img)
}


#-------------------------------------
input_img <- "example_OLCI_image_input.tif"
tss_result <- estimate_TSS(input_img)

output_img <- "example_OLCI_image_output.tif"
writeRaster(tss_result,output_img,format="GTiff")   

rm(list=ls())
gc()
