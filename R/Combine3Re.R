#' Combine3Re
#'
#' This function combines feature based differential usage results, the derived gene based results from feature based results and
#' the resutls from rMAT to identify the common results between 3 methods
#'
#' @param Re.PJ: Output from JunctionSeq
#' @param re.PJ.gene.based: Gene based results from JunctionSeq
#' @param re.rMAT: Results from rMAT
#' @param outputfile_DGE_2_FC_P_geneWise
#' @param gene.model
#' @param Output_GO_BP_feature.gene.rMAT.2.FC
#' @param Re.PJ.selected.feature.FC.p
#' @param Output_GO_BP_gene
#' @param Output_GO_BP_feature
#' @param Output_GO_BP_feature_gene
#' @param Output_GO_BP_feature_gene_redefined
#'
#' @return
#' @export
#'
#' @examples
#' Re.combine.3.methods<-Combine3Re(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model,output.dir.name.new)
#'
#' Re.combine.3.methods.rMAT.SE<-Combine3Re(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model,output.dir.name.rMAT_SE_based)
#'
#'
Combine3Re <- function(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model, Output_files_dir) {

  Re.PJ.selected.feature.2.FC.p<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,"SE",1,0.05,Output_files_dir)

  DE_type<-c("Feature","GeneWise","rMAT","FeatureGeneWise","FeaturerMAT","GeneWiserMAT","FeatureGeneWiseRMAT")

  Re.combine.1<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"Feature",gene.model,Output_files_dir)
  Re.combine.2<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"GeneWise",gene.model,Output_files_dir)
  Re.combine.3<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"rMAT",gene.model,Output_files_dir)
  Re.combine.4<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeatureGeneWise",gene.model,Output_files_dir)
  Re.combine.5<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeaturerMAT",gene.model,Output_files_dir)
  Re.combine.6<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"GeneWiserMAT",gene.model,Output_files_dir)
  Re.combine.7<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeatureGeneWiseRMAT",gene.model,Output_files_dir)

  #Re.combine<-list(rMAT=Re.combine.3)

   Re.combine<-list(Feature=Re.combine.1,
                    GeneWise=Re.combine.2,
                    rMAT=Re.combine.3,
                    FeatureGeneWise=Re.combine.4,
                    FeaturerMAT=Re.combine.5,
                    GeneWiserMAT=Re.combine.6,
                    FeatureGeneWiseRMAT=Re.combine.7)

  return(Re.combine)

}
