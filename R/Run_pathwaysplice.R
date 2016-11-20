#' Run_pathwaysplice
#'
#' Perform pathwaysplice in one step
#'
#' @param re.gene.based gene based results
#' @param ad bias factor to be adjusted
#' @param sub_feature feature to be checked
#' @param threshold threshold to be used for adjustment
#' @param genomeID gene to be used
#' @param geneID geneID to be used
#' @param gene_model gene model to be used
#' @param method method to be used
#'
#' @return a list that has gene set enrichment analysis results
#' @export
#'
#' @examples
#'
#' data(mds)
#' data(hg19)
#'
#' Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",
#' 0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Wallenius")
#'
#' Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",
#' 0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Sampling")
#'
#' Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",
#' 0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Hypergeometric")
#'
#'
Run_pathwaysplice <-
  function(re.gene.based,
           ad = "GL",
           sub_feature = NULL,
           threshold,
           genomeID,
           geneID,
           gene_model,
           method) {
    #Data4Goterm<-pData(re.gene.based)
    
    Data4Goterm <- re.gene.based
    
    if (is.null(sub_feature)) {
      Data4Goterm.sub_feature <- Data4Goterm
    }
    else{
      Data4Goterm.sub_feature <-
        Data4Goterm[grep(sub_feature, Data4Goterm[, 8]), ]
    }
    
    #print(dim(Data4Goterm.sub_feature))
    
    if (sub_feature == "J") {
      Data4Goterm.sub_feature.geneID.NumOfJunctions <-
        Data4Goterm.sub_feature[, c(1, 11)]
    } else{
      Data4Goterm.sub_feature.geneID.NumOfJunctions <-
        Data4Goterm.sub_feature[, c(1, 10)]
    }
    #print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))
    
    Data4Goterm.sub_feature.Sig <-
      Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[, 7] < threshold), ]
    
    #GO term analysis using GOSeq
    All.gene.id.based.on.sub_feature <-
      unique(Data4Goterm.sub_feature[, 1])
    #length(All.gene.id.based.on.sub_feature)
    All.gene.id.index <-
      rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) = All.gene.id.based.on.sub_feature
    
    All.genes.based.on.Sig.sub_feature <-
      unique(Data4Goterm.sub_feature.Sig[, 1])
    gene.DE_interest <-
      as.integer(
        which(
          All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature
        )
      )
    
    All.gene.id.index[gene.DE_interest] <- 1
    
    #print(length(All.gene.id.index))
    
    gene.with.matched.junction <-
      which(Data4Goterm.sub_feature.geneID.NumOfJunctions[, 1] %in% c(names(All.gene.id.index)))
    num.junction.4.matched.gene <-
      as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction, 2])
    
    #names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]
    
    #All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]
    
    #print(length(All.gene.id.index.2))
    
    All.gene.id.index.2 <- All.gene.id.index
    
    #print(All.gene.id.index.2)
    
    if (ad == "GL") {
      pwf.DE_interest = nullp(All.gene.id.index.2, genomeID, geneID, plot.fit = TRUE)
    }
    else
    {
      pwf.DE_interest = nullp(
        All.gene.id.index.2,
        genomeID,
        geneID,
        bias.data = num.junction.4.matched.gene,
        plot.fit = TRUE
      )
    }
    
    if (method == "Hypergeometric") {
      GO.wall.DE_interest = pathwaysplice(
        pwf.DE_interest,
        genomeID,
        geneID,
        gene.model = gene_model,
        method = "Hypergeometric",
        use_genes_without_cat = TRUE
      )
    } else if (method == "Sampling") {
      GO.wall.DE_interest = pathwaysplice(
        pwf.DE_interest,
        genomeID,
        geneID,
        gene.model = gene_model,
        method = "Sampling",
        use_genes_without_cat = TRUE
      )
    } else{
      GO.wall.DE_interest = pathwaysplice(
        pwf.DE_interest,
        genomeID,
        geneID,
        gene.model = gene_model,
        use_genes_without_cat = TRUE
      )
    }
    #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
    #enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
    re <-
      list(GO.wall.DE_interest = GO.wall.DE_interest, pwf.DE_interest)
    
    #re[[1]]<-GO.wall.DE_interest
    #re[[2]]<-pwf.DE_interest
    
    return(re)
  }


pathwaysplice = function(pwf,
                         genome,
                         id,
                         gene.model,
                         gene2cat = NULL,
                         test.cats = c("GO:CC", "GO:BP", "GO:MF"),
                         method = "Wallenius",
                         repcnt = 2000,
                         use_genes_without_cat = FALSE) {
  ################# Input pre-processing and validation ###################
  #Do some validation of input variables
  if (any(!test.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
    stop("Invalid category specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
  }
  if ((missing(genome) | missing(id))) {
    if (is.null(gene2cat)) {
      stop(
        "You must specify the genome and gene ID format when automatically fetching gene to GO category mappings."
      )
    }
    #If we're using user specified mappings, this obviously isn't a problem
    genome = 'dummy'
    id = 'dummy'
  }
  if (!any(method %in% c("Wallenius", "Sampling", "Hypergeometric"))) {
    stop(
      "Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric."
    )
  }
  if (!is.null(gene2cat) &&
      (!is.data.frame(gene2cat) & !is.list(gene2cat))) {
    stop(
      "Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again."
    )
  }
  
  #Factors are evil
  pwf = unfactor(pwf)
  gene2cat = unfactor(gene2cat)
  
  ###################### Data fetching and processing ########################
  if (is.null(gene2cat)) {
    #When we fetch the data using getgo it will be in the list format
    message("Fetching GO annotations...")
    gene2cat = getgo3(rownames(pwf), genome, id, fetch.cats = test.cats)
    names(gene2cat) = rownames(pwf)
    
    #cat("OK")
    #Do the two rebuilds to remove any nulls
    cat2gene = reversemapping(gene2cat)
    gene2cat = reversemapping(cat2gene)
    
    #print(cat2gene)
    #print(gene2cat)
    
  } else{
    #The gene2cat input accepts a number of formats, we need to check each of them in term
    message("Using manually entered categories.")
    #The options are a flat mapping (that is a data frame or matrix) or a list, where the list can be either gene->categories or category->genes
    if (class(gene2cat) != "list") {
      #it's not a list so it must be a data.frame, work out which column contains the genes
      genecol_sum = as.numeric(apply(gene2cat, 2, function(u) {
        sum(u %in% rownames(pwf))
      }))
      genecol = which(genecol_sum != 0)
      if (length(genecol) > 1) {
        genecol = genecol[order(-genecol_sum)[1]]
        warning(
          paste(
            "More than one possible gene column found in gene2cat, using the one headed",
            colnames(gene2cat)[genecol]
          )
        )
      }
      if (length(genecol) == 0) {
        genecol = 1
        warning(
          paste(
            "Gene column could not be identified in gene2cat conclusively, using the one headed",
            colnames(gene2cat)[genecol]
          )
        )
      }
      othercol = 1
      if (genecol == 1) {
        othercol = 2
      }
      #Now put it into our delicious listy format
      gene2cat = split(gene2cat[, othercol], gene2cat[, genecol])
      #Do the appropriate builds
      cat2gene = reversemapping(gene2cat)
      gene2cat = reversemapping(cat2gene)
    }
    #!!!!
    #The following conditional has been flagged as a potential issue when using certain
    #types of input where the category names are the same as gene names (which seems like
    #something you should avoid anyway...).  Leave it for now
    #!!!!
    #We're now garunteed to have a list (unless the user screwed up the input) but it could
    #be category->genes rather than the gene->categories that we want.
    if (sum(unique(unlist(gene2cat, use.names = FALSE)) %in% rownames(pwf)) >
        sum(unique(names(gene2cat)) %in% rownames(pwf))) {
      gene2cat = reversemapping(gene2cat)
    }
    #Alright, we're garunteed a list going in the direction we want now.  Throw out genes which we will not use
    gene2cat = gene2cat[names(gene2cat) %in% rownames(pwf)]
    
    #Rebuild because it's a fun thing to do
    cat2gene = reversemapping(gene2cat)
    gene2cat = reversemapping(cat2gene)
    
    ## make sure we remove duplicate entries .. e.g. see
    ## http://permalink.gmane.org/gmane.science.biology.informatics.conductor/46876
    cat2gene = lapply(cat2gene, function(x) {
      unique(x)
    })
    gene2cat = lapply(gene2cat, function(x) {
      unique(x)
    })
  }
  
  nafrac = (sum(is.na(pwf$pwf)) / nrow(pwf)) * 100
  if (nafrac > 50) {
    warning(
      paste(
        "Missing length data for ",
        round(nafrac),
        "% of genes.  Accuarcy of GO test will be reduced.",
        sep = ''
      )
    )
  }
  #Give the genes with unknown length the weight used by the median gene (not the median weighting!)
  pwf$pwf[is.na(pwf$pwf)] = pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data)) /
                                                                                               2)], pwf$bias.data)]
  
  ###################### Calculating the p-values ########################
  # Remove all the genes with unknown GOterms
  unknown_go_terms = nrow(pwf) - length(gene2cat)
  if ((!use_genes_without_cat) && unknown_go_terms > 0) {
    message(
      paste(
        "For",
        unknown_go_terms,
        "genes, we could not find any categories. These genes will be excluded."
      )
    )
    message(
      "To force their use, please run with use_genes_without_cat=TRUE (see documentation)."
    )
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf = pwf[rownames(pwf) %in% names(gene2cat), ]
  }
  #A few variables are always useful so calculate them
  cats = names(cat2gene)
  DE = rownames(pwf)[pwf$DEgenes == 1]
  num_de = length(DE)
  num_genes = nrow(pwf)
  pvals = data.frame(
    category = cats,
    over_represented_pvalue = NA,
    under_represented_pvalue = NA,
    stringsAsFactors = FALSE,
    numDEInCat = NA,
    numInCat = NA
  )
  if (method == "Sampling") {
    #We need to know the number of DE genes in each category, make this as a mask that we can use later...
    num_DE_mask = rep(0, length(cats))
    a = table(unlist(gene2cat[DE], FALSE, FALSE))
    
    num_DE_mask[match(names(a), cats)] = as.numeric(a)
    num_DE_mask = as.integer(num_DE_mask)
    #We have to ensure that genes not associated with a category are included in the simulation, to do this they need an empty entry in the gene2cat list
    gene2cat = gene2cat[rownames(pwf)]
    names(gene2cat) = rownames(pwf)
    message("Running the simulation...")
    #Now do the actual simulating
    lookup = matrix(0, nrow = repcnt, ncol = length(cats))
    for (i in 1:repcnt) {
      #A more efficient way of doing weighted random sampling without replacment than the built in function
      #The order(runif...)[1:n] bit picks n genes at random, weighting them by the PWF
      #The table(as.character(unlist(...))) bit then counts the number of times this random set occured in each category
      a = table(as.character(unlist(gene2cat[order(runif(num_genes) ^ (1 /
                                                                         pwf$pwf), decreasing = TRUE)[1:num_de]], FALSE, FALSE)))
      lookup[i, match(names(a), cats)] = a
      pp(repcnt)
    }
    message("Calculating the p-values...")
    #The only advantage of the loop is it uses less memory...
    #for(i in 1:length(cats)){
    #	pvals[i,2:3]=c((sum(lookup[,i]>=num_DE_mask[i])+1)/(repcnt+1),(sum(lookup[,i]<=num_DE_mask[i])+1)/(repcnt+1))
    #	pp(length(cats))
    #}
    pvals[, 2] = (colSums(lookup >= outer(rep(1, repcnt), num_DE_mask)) +
                    1) / (repcnt + 1)
    pvals[, 3] = (colSums(lookup <= outer(rep(1, repcnt), num_DE_mask)) +
                    1) / (repcnt + 1)
  }
  if (method == "Wallenius") {
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum = which(pwf$DEgenes == 1)
    #Turn all genes into a reference to the pwf object
    
    cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
    #This value is used in every calculation, by storing it we need only calculate it once
    alpha = sum(pwf$pwf)
    
    #Each category will have a different weighting so needs its own test
    pvals[, 2:3] = t(sapply(cat2genenum, function(u) {
      #The number of DE genes in this category
      num_de_incat = sum(degenesnum %in% u)
      
      #The total number of genes in this category
      num_incat = length(u)
      
      #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
      avg_weight = mean(pwf$pwf[u])
      weight = (avg_weight * (num_genes - num_incat)) / (alpha - num_incat *
                                                           avg_weight)
      if (num_incat == num_genes) {
        weight = 1
      } #case for the root GO terms
      #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
      
      #cat(num_de_incat,"\t",num_incat,"\t",num_genes,"\t",num_de,"\t",weight,"\n")
      
      c(
        dWNCHypergeo(
          num_de_incat,
          num_incat,
          num_genes - num_incat,
          num_de,
          weight
        )
        + pWNCHypergeo(
          num_de_incat,
          num_incat,
          num_genes - num_incat,
          num_de,
          weight,
          lower.tail = FALSE
        ),
        pWNCHypergeo(
          num_de_incat,
          num_incat,
          num_genes - num_incat,
          num_de,
          weight
        )
      )
    }))
  }
  if (method == "Hypergeometric") {
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum = which(pwf$DEgenes == 1)
    #Turn all genes into a reference to the pwf object
    cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
    #Simple hypergeometric test, one category at a time
    pvals[, 2:3] = t(sapply(cat2genenum, function(u) {
      #The number of DE genes in this category
      num_de_incat = sum(degenesnum %in% u)
      #The total number of genes in this category
      num_incat = length(u)
      #Calculate the sum of the tails of the hypergeometric distribution (the p-values)
      c(
        dhyper(num_de_incat, num_incat, num_genes - num_incat, num_de) + phyper(
          num_de_incat,
          num_incat,
          num_genes - num_incat,
          num_de,
          lower.tail = FALSE
        ),
        phyper(num_de_incat, num_incat, num_genes - num_incat, num_de)
      )
    }))
  }
  
  #Populate the count columns...
  degenesnum = which(pwf$DEgenes == 1)
  cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
  pvals[, 4:5] = t(sapply(cat2genenum, function(u) {
    c(sum(degenesnum %in% u), length(u))
  }))
  
  #Got the name of DE gene in each GO
  #print(head(cat2gene))
  
  DE_pwf = rownames(pwf[degenesnum, ])
  
  #cat2degenenum=relist(match(unlist(cat2gene),rownames(pwf[degenesnum,])),cat2gene)
  #print(cat2degenenum)
  
  pvals.6 <- sapply(cat2gene, function(u, DE_pwf) {
    #c(sum(degenesnum%in%u),length(u))
    #c(rownames(pwf)[u[-which(is.na(u))]])
    x <- u[which(u %in% DE_pwf)]
    x
  }, DE_pwf)
  
  #cat("After matching\n")
  
  #cat("DE_ensemble\n")
  
  #print(unique(as.character(unlist2(pvals.6))))
  #cat(length(unique(as.character(unlist2(pvals.6)))),"\n")
  
  pvals.6.gene.symbol <- sapply(pvals.6, function(u, gene.model) {
    y <- gene.model[which(as.character(gene.model[, 3]) %in% u), 1]
    y
  }, gene.model)
  
  #cat("DE_symbol\n")
  #print(unique(as.character(unlist2(pvals.6.gene.symbol))))
  #cat(length(unique(as.character(unlist2(pvals.6.gene.symbol)))),"\n")
  
  #Convert list to data frame
  pvals.6.df <- list_to_df(pvals.6)
  #cat("pvals_6_df","\n")
  #cat(dim(pvals.6.df)[1],"\n")
  
  pvals.6.gene.symbol.df <- list_to_df(pvals.6.gene.symbol)
  #cat("pvals_6_gene_symbol_df","\n")
  #cat(dim(pvals.6.gene.symbol.df)[1],"\n")
  
  #cat(length(unique(as.character(pvals.6.gene.symbol.df[,2]))),"\n")
  #print(unique(as.character(pvals.6.gene.symbol.df[,2])))
  
  dataset2 <- pvals.6.gene.symbol.df
  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)
             sapply(x, function(y)
               paste(unlist(y), collapse = ", ")))
  
  #print(dataset2)
  #df_args <- c(pvals.6.gene.symbol.df[,2], sep=",")
  #df_args.data<-do.call(paste, df_args)
  
  #print(class(dataset2))
  #cat(dim(dataset2),"\n")
  
  temp.gene.name = unique(apply(dataset2[, 2], 1, c))
  temp.gene.name.2 = unique(gdata::trim(unlist(strsplit(
    temp.gene.name, split = ","
  ))))
  #print(class(temp.gene.name.2))
  
  #cat(length(temp.gene.name.2),"\n")
  #print(temp.gene.name.2)
  DE_from_GO <- temp.gene.name.2
  
  colnames(pvals.6.df) = c("category", "DEgene_ID")
  colnames(pvals.6.gene.symbol.df) = c("category", "DEgene_symbol")
  
  #Finally, sort by p-value
  pvals = pvals[order(pvals$over_represented_pvalue), ]
  
  # Supplement the table with the GO term name and ontology group
  # but only if the enrichment categories are actually GO terms
  if (any(grep("^GO:", pvals$category))) {
    GOnames = select(GO.db,
                     keys = pvals$category,
                     columns = c("TERM", "ONTOLOGY"))[, 2:3]
    colnames(GOnames) <- tolower(colnames(GOnames))
    pvals = cbind(pvals, GOnames)
  }
  
  # And return
  pvals.2 <- merge(pvals, pvals.6.df, by = "category", sort = FALSE)
  
  pvals.3 <-
    merge(pvals.2,
          pvals.6.gene.symbol.df,
          by = "category",
          sort = FALSE)
  
  pvals.4 <- list(GO = pvals.3, DE_GO = DE_from_GO)
  
  return(pvals.4)
  
}

getgo3 = function(genes,
                  genome,
                  id,
                  fetch.cats = c("GO:CC", "GO:BP", "GO:MF")) {
  #Check for valid input
  if (any(!fetch.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
    stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
  }
  #Convert from genome ID to org.__.__.db format
  orgstring = as.character(.ORG_PACKAGES[match(gsub("[0-9]+", '', genome), names(.ORG_PACKAGES))])
  #Multimatch or no match
  if (length(orgstring) != 1) {
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  #Load the library
  library(paste(orgstring, "db", sep = '.'), character.only = TRUE)
  #What is the default ID that the organism package uses?
  coreid = strsplit(orgstring, "\\.")[[1]][3]
  
  #Now we need to convert it into the naming convention used by the organism packages
  userid = as.character(.ID_MAP[match(id, names(.ID_MAP))])
  #Multimatch or no match
  if (is.na(userid) | (length(userid) != 1)) {
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  #The (now loaded) organism package contains a mapping between the internal ID and whatever
  #the default is (usually eg), the rest of this function is about changing that mapping to
  #point from categories to the ID specified
  #Fetch the mapping in its current format
  #Because GO is a directed graph, we need to get not just the genes associated with each ID,
  #but also those associated with its children.  GO2ALLEGS does this.
  core2cat = NULL
  if (length(grep("^GO", fetch.cats)) != 0) {
    #Get the name of the function which maps gene ids to go terms
    #usually this will be "GO2ALLEG"
    gomapFunction = .ORG_GOMAP_FUNCTION[orgstring]
    if (is.na(gomapFunction))
      gomapFunction = .ORG_GOMAP_FUNCTION["default"]
    x = toTable(get(paste(orgstring, gomapFunction, sep = '')))
    #Keep only those ones that we specified and keep only the names
    #		core2cat=x[x$Ontology%in%gsub("^GO:",'',fetch.cats),1:2]
    x[!x$Ontology %in% gsub("^GO:", '', fetch.cats), 2] <- "Other"
    core2cat = x[, 1:2]
    colnames(core2cat) = c("gene_id", "category")
  }
  if (length(grep("^KEGG", fetch.cats)) != 0) {
    x = toTable(get(paste(orgstring, "PATH", sep = '')))
    #Either add it to existing table or create a new one
    colnames(x) = c("gene_id", "category")
    if (!is.null(core2cat)) {
      core2cat = rbind(core2cat, x)
    } else{
      core2cat = x
    }
  }
  
  #Now we MAY have to convert the "gene_id" column to the format we are using
  if (coreid != userid) {
    #The mapping between user id and core id, don't use the <USER_ID>2<CORE_ID> object as the naming is not always consistent
    user2core = toTable(get(paste(orgstring, userid, sep = '')))
    #Throw away any user ID that doesn't appear in core2cat
    user2core = user2core[user2core[, 1] %in% core2cat[, 1], ]
    #Make a list version of core2cat, we'll need it
    list_core2cat = split(core2cat[, 2], core2cat[, 1])
    #Now we need to replicate the core IDs that need replicating
    list_core2cat = list_core2cat[match(user2core[, 1], names(list_core2cat))]
    #Now we can replace the labels on this list with the user ones from user2core,
    #but there will be duplicates, so we have to unlist, label, then relist
    user2cat = split(unlist(list_core2cat, FALSE, FALSE),
                     rep(user2core[, 2], sapply(list_core2cat, length)))
    #Now we only want each category listed once for each entry...
    user2cat = sapply(user2cat, unique)
    ###In case you don't believe that this works as it should, here is the slow as all hell way for comparison...
    ###Make first list
    ##list_user2core=split(user2core[,1],user2core[,2])
    ###Make the second
    ##list_core2cat=split(core2cat[,2],core2cat[,1])
    ###Go through each entry in first list and expand using second...
    ##user2cat=sapply(list_user2core,function(u){unique(unlist(list_core2cat[u],FALSE,FALSE))})
    
  } else{
    #We don't need to convert anything (WOO!), so just make it into a list
    user2cat = split(core2cat[, 2], core2cat[, 1])
    user2cat = sapply(user2cat, unique)
  }
  #remove any empty strings
  user2cat = lapply(user2cat, function(x) {
    if (length(x) > 1)
      x = x[x != "Other"]
    x
  })
  
  ## we don't like case sensitivity
  names(user2cat) <- toupper(names(user2cat))
  
  #Now look them up
  return(user2cat[toupper(genes)])
}

#Description: Prints progress through a loop
#copy from Matthew Young's goseq

pp = function(total, count, i = i) {
  if (missing(count)) {
    count = evalq(i, envir = parent.frame())
  }
  if (missing(total)) {
    total = evalq(stop, envir = parent.frame())
  }
  cat(round(100 * (count / total)), "%   \r")
}
