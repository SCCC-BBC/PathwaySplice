#' Title
#'
#' @param jscs.object
#' @param gene.name
#' @param outfile_prefix
#'
#' @return
#' @export
#'
#' @examples
#' OutputFigure(jscs.object=Re.PJ,gene.name="ENSMUSG00000000028+ENSMUSG00000005262",
#' outfile_prefix="/media/H_driver/PJ/plot_ENSMUSG00000000028_ENSMUSG00000005262/")
#'
OutputFigure <- function(jscs.object,gene.name,outfile_prefix) {
  buildAllPlots(jscs=Re.PJ,outfile.prefix=outfile_prefix,gene.list=gene.name,
                use.plotting.device="png",plot.gene.level.expression=TRUE);
}
