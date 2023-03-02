# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get genes in pathways
#'
#' @description Get a list of all genes found in a list of pathways (using HumanMine via InterMineR package).
#' This function take a vector of Reactome pathway identifiers and calls HumanMine's template query 'PathwayGenes' which returns all genes which are found in each of the inputted pathways.
#' 
#' @param pathways Character vector of Reactome pathways identifiers to be searched
#'
#' @import InterMineR
#' @importFrom glue glue
#' @importFrom logger log_fatal log_warn
#' @include global_functions.R
#' @return A data.frame of genes symbol and gene HumanMine identifiers with pathway identifier which the genes belong to. 
#' @export
getGenes <- function(pathways) {
  checkClass(parameter = "pathways", value = pathways, classType = "character", functionName = "getGenes")
  pathToGeneResult <- runHumanMineQuery(templateName = "PathwayGenes",
                                        constraintPath = c("Pathway.identifier", "Pathway.genes.organism.shortName"),
                                        constraintOperators = c("=", "="),
                                        constraintValues = c(list(as.character(pathways)), "H. sapiens"),
                                        constraintm.index = c(1, 2),
                                        select = c("Pathway.genes.primaryIdentifier",
                                                   "Pathway.genes.symbol",
                                                   "Pathway.identifier")
  )
  if(is.null(pathToGeneResult)) {
    e <- "None of the inputted pathways were found as Reactome pathways in HumanMine. Please double check the identifiers"
    logger::log_fatal(e)
    stop(e)
  } else {
    colnames(pathToGeneResult) <- c("Gene Identifier", "Gene Symbol", "Pathway Identifier")
    pathsNotFound <- pathways[!pathways %in% pathToGeneResult$`Pathway Identifier`]
    if(length(pathsNotFound) > 0) {
      e <- glue::glue()
      logger::log_warn("{pathsNotFound} input pathway(s) not mapped to any Human Reactome Pathways. Please double check pathway identifier.")
    }
  }
  return(pathToGeneResult)
}


