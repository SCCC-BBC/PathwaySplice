#' Generate simulated count based on NB distribution
#'
#' @param dataset
#' @param group
#' @param nTags
#' @param nlibs
#' @param fix.dispersion
#' @param lib.size
#' @param drop.low.lambda
#' @param drop.extreme.dispersion
#' @param add.outlier
#' @param outlierMech
#' @param pOutlier
#' @param min.factor
#' @param max.factor
#' @param pDiff
#' @param pUp
#' @param foldDiff
#' @param name
#' @param save.file
#' @param file
#' @param only.add.outlier
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
NBsim <-
  function(dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 10, pDiff=.1, pUp=.5, foldDiff=3, name = NULL, save.file = FALSE, file = NULL, only.add.outlier = FALSE, verbose=TRUE)

  {
    ## NBsim generate simulated count from the real dataset followed by the NB model ##
    require(edgeR)
    group = as.factor(group)

    sample.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it samples from the real dataset ##
      nlibs <- object$nlibs
      nTags <- object$nTags
      AveLogCPM <-object$dataset$dataset.AveLogCPM
      dispersion <- object$dataset$dataset.dispersion

      id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
      Lambda <- 2^(AveLogCPM[id_r])
      Lambda <- Lambda/sum(Lambda)
      Dispersion <- dispersion[id_r]
      id_0<- Lambda == 0
      Lambda <- Lambda[!id_0]
      Dispersion <- Dispersion[!id_0]
      Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
      object$Lambda <- Lambda
      if(!is.na(fix.dispersion))
        Dispersion <- expandAsMatrix(fix.dispersion, dim = c(nTags, nlibs))
      else Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
      object$Dispersion <- Dispersion
      object

    }
    diff.fun <- function(object)
    {

      ## it is low-level function of NBsim ##
      ## it creates diff genes according to foldDiff ##
      group <- object$group
      pDiff <- object$pDiff
      pUp <-  object$pUp
      foldDiff <- object$foldDiff
      Lambda <- object$Lambda
      nTags <- object$nTags
      g <- group == levels(group)[1]
      ind <- sample(nTags, floor(pDiff*nTags))
      if(length(ind)>0 & !foldDiff == 1 ) {
        fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
        Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
        Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir))
        #Lambda <- t(t(Lambda)/colSums(Lambda))
        object$Lambda <- Lambda
        object$indDE <- ind
        object$indNonDE <- (1:nTags)[-ind]
        object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda))
        object$mask_DEup[ind[fcDir == 1], g] <- TRUE
        object$mask_DEup[ind[fcDir == -1], !g] <- TRUE
        object$mask_DEdown[ind[fcDir == -1], g] <- TRUE
        object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE
        object$mask_nonDE[-ind,] <- TRUE
        object$mask_DE <- object$mask_DEup | object$mask_DEdown}
      if(foldDiff == 1| pDiff == 0)
        object$indDE <- NA
      object
    }
    sim.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it simulate counts using rnbinom ##
      Lambda <- object$Lambda
      Dispersion <- object$Dispersion
      nTags <- object$nTags
      nlibs <- object$nlibs
      lib.size <- object$lib.size
      counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs)
      rownames(counts) <- paste("ids", 1:nTags, sep = "")
      object$counts <- counts
      object
    }

    outlier.fun <- function(object, outlierMech, pOutlier, min.factor = 2, max.factor = 5)
    {
      ## it is low-level function of NBsim ##
      ## it makes outlier ##
      outlierMech <- match.arg(outlierMech, c("S", "M", "R"))
      dim <- dim(object$counts)
      outlier.factor <- function() runif(1, min.factor, max.factor)
      countAddOut <- object$counts
      LambdaAddOut <- object$Lambda
      DispersionAddOut <- object$Dispersion
      switch(outlierMech,
             S = {
               mask_outlier <- expandAsMatrix(FALSE, dim = dim)
               id_r <- which(runif(dim[1]) < pOutlier)
               id_c <- sample(dim[2], length(id_r), replace = TRUE)
               for(i in seq(id_r))
                 mask_outlier[id_r[i], id_c[i]] <- TRUE
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },
             R = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },

             M = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               LambdaAddOut[mask_outlier] <- sapply(LambdaAddOut[mask_outlier], function(z) z*outlier.factor())
               countAddOut[mask_outlier] <- rnbinom(sum(mask_outlier), mu = t(t(LambdaAddOut)*object$lib.size)[mask_outlier], size = 1/DispersionAddOut[mask_outlier])
             }
      )
      if(!object$foldDiff == 1 & !pDiff == 0)
      {
        indDEupOutlier <- which(apply(object$mask_DEup & mask_outlier, 1, any))
        indDEdownOutlier <- which(apply(object$mask_DEdown & mask_outlier, 1, any))
        indDEnoOutlier <- which(apply((object$mask_DE & !mask_outlier) , 1, all))
        indNonDEOutlier <- which(apply(object$mask_nonDE & mask_outlier, 1, any))
        indNonDEnoOutlier <- which(apply((object$mask_nonDE & !mask_outlier) , 1, all))
        indDEbothOutlier <- NA
        o <- indDEupOutlier %in% indDEdownOutlier
        q <-  indDEdownOutlier %in% indDEupOutlier
        if(any(o))
        {
          indDEupOutlier <- indDEupOutlier[!o]
          indDEbothOutlier <- indDEupOutlier[o]
          indDEdownOutlier <- indDEdownOutlier[!q]
        }
      }
      else
      {
        indDEupOutlier <- indDEdownOutlier <- indDEnoOutlier <- indNonDEOutlier <- indNonDEnoOutlier <- indDEbothOutlier <- NA
      }
      out <- list(countAddOut = countAddOut, outlierMech = outlierMech, pOutlier = pOutlier, mask_outlier = mask_outlier, indDEupOutlier = indDEupOutlier,
                  indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier,
                  indNonDEnoOutlier = indNonDEnoOutlier, LambdaAddOut = LambdaAddOut, DispersionAddOut = DispersionAddOut)

    }

    calProb <- function(x, l) round(1 -(1 - x)^(1/l), 2) ## calculate probability to make sure all the outlierMech produce the same amount of outliers ##


    if(verbose) message("Preparing dataset.\n")
    if(class(dataset) == "DGEList")
    {
      dat <- dataset
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.character(dataset))
    {
      load(dataset)
      dat <- get(gsub("(\\.)(\\w+)", "", basename(dataset)))
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.matrix(dataset))
    {
      if(is.null(name)) name <- deparse(substitute(dataset))
      dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))
    }
    else
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))

    if(!only.add.outlier)
    {
      if(is.null(lib.size))
        dat$lib.size <- runif(nlibs, min = 0.7*median(dat$dataset$dataset.lib.size), max = 1.3*median(dat$dataset$dataset.lib.size))

      if(is.null(nTags))
        dat$nTags <- dat$dataset$dataset.nTags
      if(verbose) message("Sampling.\n")
      dat <- sample.fun(dat)
      if(verbose) message("Calculating differential expression.\n")
      dat <- diff.fun(dat)
      if(verbose) message("Simulating data.\n")
      dat <- sim.fun(dat)
    }
    if(add.outlier){
      outlierMech <- match.arg(outlierMech,  c("S", "R", "M"), several.ok = TRUE)
      if(length(pOutlier)== 1 & length(outlierMech) > 1 & any(outlierMech == "S"))
      {
        prob <- calProb(pOutlier, length(group))
        pOutlier <- rep(pOutlier, length = length(outlierMech))
        pOutlier[!outlierMech == "S"] <- prob
      }
      else if(!length(pOutlier) == length(outlierMech))
        stop("pOutlier is not equal to outlierMech")
      if(verbose) message("Adding outliers.\n")
      dat[outlierMech] <- mapply(function(x, y) outlier.fun(dat, outlierMech = x, pOutlier = y, min.factor = min.factor, max.factor = max.factor), x = outlierMech, y = pOutlier, SIMPLIFY = FALSE)
      dat$pOutlier <- pOutlier
    }

    if(save.file)
    {

      ## save file for shiny app ##
      if(verbose) message("Saving file.\n")
      if(is.null(file))
      { g <- paste0("g", sum(levels(group)[1] == group), "v", sum(levels(group)[2] == group))
      f <- paste0("f", foldDiff)
      if(add.outlier) o <- paste0("o", sprintf( "%02d",100*pOutlier[1L]))
      else o <- paste0("o", sprintf( "%02d", 0 ))
      file <- paste0(dat$name, "/", g, f, o, ".Rdata")
      dir.create(dat$name, showWarnings = FALSE)
      }
      filenm <- eval(gsub("(\\.)(\\w+)", "", basename(file)))
      assign(filenm, dat)
      save(list = filenm, file = file)
    }
    dat
  }
