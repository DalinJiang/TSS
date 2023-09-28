# Readme

## Introduction

Those R code are for estimating total suspended solids (TSS) concentration from remote sensing data

There are two folders including code for:

(1) MERIS, OLCI image

(2) MSI image

in each folder, there are example input/output files, which can be used to test running the code

input: Remote sensing reflectance (sr-1), from either in situ measured spectra or atmospheric corrected satellite image

output: (i)absorption coefficient, (ii)backscattering coefficient, (iii)reference band (indicating water type), (iv)TSS
#
#

## How to use

(1) copy all the files from Github to your computer;

(2) open the R file "Run_code_to_get_TSS.R", change the working directory to the path you save the code;

(3) save your input file (either csv or tiff file) in the same folder;

(4) run the code "Run_code_to_get_TSS.R" to get TSS results.
#
#

## References

(1) Jiang, D., Matsushita, B., Pahlevan, N., Gurlin, D., Lehmann, M. K., Fichot, C. G., ... & O'Donnell, D. (2021). Remotely estimating total suspended solids concentration in clear to extremely turbid waters using a novel semi-analytical method. Remote sensing of environment, 258, 112386. DOI: https://doi.org/10.1016/j.rse.2021.112386

(2) Jiang, D., Matsushita, B., Pahlevan, N., Gurlin, D.,Fichot, C. G., ... & Spyrakos, E. (2023). Estimating the concentration of total suspended solids in inland and coastal waters from Sentinel-2 MSI: A semi-analytical approach. ISPRS Journal of Photogrammetry and Remote Sensing, 204, 362-377. DOI: https://doi.org/10.1016/j.isprsjprs.2023.09.020



![Picture 1](https://github.com/DalinJiang/TSS/assets/117453464/4550236e-f43c-4017-9272-69f9e9285bcd)


