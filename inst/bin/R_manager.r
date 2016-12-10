#!/usr/bin/Rscript 

cat("Do you want to build this package?\n")

input<-file('stdin', 'r')
row <- readLines(input, n=1)

print(row)

if(row=="Yes") {
  
  cat("Starting build...\n")
  system("R CMD build --resave-data PathwaySplice")
  cat("Finished build ...\n")
  
  cat("Do you want to perform CRAN check for this package?\n")
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("Starting check as CRAN...\n")
    system("R CMD check --as-cran PathwaySplice_0.99.0.tar.gz")
    cat("Finished check as CRAN...\n")
  }
  
  cat("Do you want to perform BiocCheck for this package?\n")
  
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("Starting BiocCheck...\n")
    system("R CMD BiocCheck PathwaySplice_0.99.0.tar.gz")
    cat("Finished BiocCheck...\n")
  }
  
  
  cat("Do you want to install this package?\n")
  
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("Starting install...\n")
    system("R CMD INSTALL PathwaySplice_0.99.0.tar.gz -l /Users/axy148/Library/R/3.3/library")
    cat("Finished install...\n")
  }
  
}
