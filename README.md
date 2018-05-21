Custom codes for identifying the association between urea cycle dysregulation and mutational bias

This repository provides the code used for our manuscript J.S. Lee, L. Adler et al, Urea cycle dysregulation in cancer results in a mutation bias associated with enhanced response to immune therapy (submitted to Cell), a collaboration between Eytan Ruppin and Ayelet Erez lab. The code was tested using R version 3.3.1 (2016-06-21), on a x86 64-pc-linux-gnu (64-bit) platform, using libraries caTools v1.17.1, ROCR v1.0-7, and data.table v1.10.4.

The codes involve two analyses by and large, immune checkpoint therapy analysis to evaluate the association of UCD and PTMB with enhanced response to checkpoint therapies (UCD.ICT.r), and peptidomics analysis to establish the immunogenecity of antigens arising from PTMB or UCD (UCD.peptidomics.r).

1. Download the code
```
get clone https://github/jooslee/UCD.git
```
2. Move to the root folder and run R
```
cd UCD
R
```
3. Install required R libraries
```
> install.packages("data.table")
> install.packages("ROCR")
> install.packages("caTools")
```
4. Launch the code of interest
```
> setwd("UCD")
> source("./R/UCD.ICT.r")
> source("./R/UCD.peptidomics.r")
