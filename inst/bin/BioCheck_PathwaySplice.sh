#!/bin/bash

R CMD build --resave-data PathwaySplice

R CMD check --as-cran PathwaySplice_0.99.0.tar.gz
R CMD BiocCheck PathwaySplice_0.99.0.tar.gz

#check the content of PathwaySplice
#tar -ztvf PathwaySplice_0.99.0.tar.gz

#install PathwaySplice
#R CMD INSTALL PathwaySplice_0.99.0.tar.gz
