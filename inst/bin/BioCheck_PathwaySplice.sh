#!/bin/bash

R CMD build --resave-data PathwaySplice
#tar -ztvf PathwaySplice_0.1.0.tar.gz
#R CMD check --as-cran PathwaySplice_0.1.0.tar.gz
R CMD BiocCheck PathwaySplice_0.0.1.tar.gz
#R CMD INSTALL PathwaySplice_0.1.0.tar.gz
