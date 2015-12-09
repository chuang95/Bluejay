Drug response classier

R package use for classify ovarian cancer patients drug response by CA125 profile


This function classifies drug response by using CA125 history profile


Requirement

library ggplot2
library grid


Usage

Bluejay.classify(input file, Day_limit=200,CA125_limit=500,CA125_bound=35)


Arguments

inputfile       CA125 pd file
Day_limit       time period, this function will classify patient base on the CA125 history of these days + 21 (3 weeks) after surgery
CA125_limit     CA125 upper limit when plotting
CA125_bound     CA125 value lower bound, we will consider this patient temporary cured if her CA125 value lower than this number

return a drug response label: non-determined, sensitive, resistant, CA125 decrease caused by surgery


Install

Download zip file and use install.package function in R to install
R CMD INSTALL --merge-multiarch Bluejay_*.tar.gz


Example

library(ggplot2)
library(grid)
library(Bluejay)
Bluejay.classify("OCI 1021.pd")


pipeline bash script

R -q -e "library(Bluejay);library(ggplot2);library(grid);Bluejay.classify("OCI 1021.pd")"


