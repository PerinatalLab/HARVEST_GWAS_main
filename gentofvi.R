#!/usr/bin/Rscript

library(DatABEL)
library(GenABEL)
source("impute2databel.R")
args=commandArgs(TRUE)

cv = impute2databel(genofile=args[1], samplefile=args[2], outfile=args[3], makeprob=FALSE)
