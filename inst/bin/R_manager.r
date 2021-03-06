#!/usr/bin/Rscript 

#Usage: 
 
#Rscript  Rscript ~/PathwaySplice/inst/bin/R_manager.r

cat("Do you want to build the manual of this package?\n")

input<-file('stdin', 'r')
row <- readLines(input, n=1)
print(row)

if(row=="Yes") {
  
  cat("Starting build...\n")
  pack <- "PathwaySplice"
  path <- find.package(pack)
  system(paste(shQuote(file.path(R.home("bin"), "R")),
               "CMD", "Rd2pdf", shQuote(path)))
  cat("Finished build ...\n")
}

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
    
    cat("Please specify the latest package that you build...\n")
    input<-file('stdin', 'r')
    
    input_pkg <- readLines(input, n=1)
    
    cat("Starting check as CRAN...\n")
    cmd=paste0("R CMD check --no-build-vignettes --as-cran ",input_pkg)
    system(cmd)
    cat("Finished check as CRAN...\n")
  }else{
    cat("You decide Not check this package as CRAN at the moment\n")
    quit()
  }
  
  cat("Do you want to perform BiocCheck for this package?\n")
  
  input<-file('stdin', 'r')
  row <- readLines(input, n=1)
  
  if(row=="Yes"){
    
    cat("Please specify the latest package that you build...\n")
    input<-file('stdin', 'r')
    
    input_pkg <- readLines(input, n=1)
    
    cmd=paste0("R CMD BiocCheck ",input_pkg)
    system(cmd)
    
    cat("Starting BiocCheck...\n")
    system(cmd)
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
    
    cat("Please specify the latest package that you build...\n")
    input<-file('stdin', 'r')
    
    input_pkg <- readLines(input, n=1)
  
    if(row==1){
    
      R_lib=.libPaths()[1]
      cat(paste0("Starting install at ",R_lib,"...\n"))
      cmd=paste0("R CMD INSTALL ",input_pkg," -l ",R_lib) 
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
