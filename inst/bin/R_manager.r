#!/usr/bin/Rscript 

#Usage: 
 
#Rscript  Rscript ~/PathwaySplice/inst/bin/R_manager.r

cat("Do you want to build this package?\n")

input<-file('stdin', 'r')
row <- readLines(input, n=1)

print(row)

if(row=="Yes") {
  
  cat("Starting build...\n")
  system("R CMD build --resave-data --no-build-vignettes PathwaySplice")
  cat("Finished build ...\n")
  
  cat("Do you want to perform CRAN check for this package?\n")
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("Starting check as CRAN...\n")
    system("R CMD check --no-build-vignettes --as-cran PathwaySplice_0.99.0.tar.gz")
    cat("Finished check as CRAN...\n")
  }else{
    cat("You decide Not check this package as CRAN at the moment\n")
    quit()
  }
  
  cat("Do you want to perform BiocCheck for this package?\n")
  
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("Starting BiocCheck...\n")
    system("R CMD BiocCheck PathwaySplice_0.99.0.tar.gz")
    cat("Finished BiocCheck...\n")
  }else{
    cat("You decide Not BiocCheck for this package at the moment\n")
    quit()
  }
  
  cat("Do you want to install this package?\n")
  
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    cat("please defne the library path from the following list:\n")
    print(.libPaths())
    
    input<-file('stdin', 'r')
    row <- readLines(input, n=1)
    
    if(row==1){
      R_lib=.libPaths()[1]
      cat(paste0("Starting install at ",R_lib,"...\n"))      
      cmd=paste0("R CMD INSTALL PathwaySplice_0.99.0.tar.gz -l ",R_lib) 
      system(cmd)
    }
    cat("Finished install...\n")
  }else{
    cat("You decide Not install this package at the moment\n")
    quit()
  }
  
}else{
  cat("You decide Not build this package at the moment\n")
  quit()
}