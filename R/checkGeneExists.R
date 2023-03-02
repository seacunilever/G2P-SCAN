# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Check Gene Symbols Exist
#'
#' @description Search HumanMine to determine if gene symbols are found or if synonyms exist, 
#'
#' @details
#' This function calls HumanMines Gene_description template to search for a list of input genes. The template will search using a lookup so if a gene symbol is a synonym of another gene a result will still be return.
#' This function will identify if the input genes are found as the official gene symbol, a synonym of another gene (if so return the synonym) or if the gene is not found in form of a data.frame output.
#' @param genes character vector of gene symbols to be searched
#' @importFrom logger log_info
#' @importFrom glue glue
#' @include global_functions.R
#' @return data.frame of search gene, column Found (TRUE, FALSE, Synonym) and column Synonym giving the synonym gene symbol if available.
#' @export
checkGenesExists <- function(genes) {
  checkClass(parameter = "genes", value = genes, classType = "character", functionName = "checkGenesExists")
      gene_search <- runHumanMineQuery(templateName = "Gene_description",
                    constraintPath = "Gene",
                    constraintOperators = "LOOKUP",
                    constraintValues = list(genes),
                    constraintm.index = 2,
                    select = c("Gene.primaryIdentifier",
                               "Gene.symbol",
                               "Gene.synonyms.value")
                      )
    
    gene_reduced <- unique(gene_search[gene_search$Gene.synonyms.value %in% genes, ])
    genes_synonyms_only <- gene_reduced[gene_reduced$Gene.symbol!= gene_reduced$Gene.synonyms.value, ]
    found_genes <-  gene_reduced[gene_reduced$Gene.symbol== gene_reduced$Gene.synonyms.value, ]
    not_found_genes <- genes[!genes %in% gene_reduced$Gene.synonyms.value]
    
    if(is.null(gene_search)) {
      e <- glue::glue("All input genes - {paste0(genes, collapse = ', ' )} - are NOT found gene symbols.")
      logger::log_fatal(e)
      stop(e)
    } else {
    if(nrow(genes_synonyms_only) > 0) {
    logger::log_info("{genes_synonyms_only$Gene.synonyms.value} is not a found gene as it's a SYNONYM of {genes_synonyms_only$Gene.symbol} and this gene should be searched instead.")
    } 
    if(nrow(found_genes) > 0) {
      logger::log_info("{found_genes$Gene.symbol} is a FOUND gene symbol.")
    } 
    if(length(not_found_genes) > 0) {
      logger::log_info("{not_found_genes} is NOT a found gene symbol.")
    }
    }
    
    found_dt <- data.frame("Gene" = genes, "Found" = FALSE, "Synonym" = NA)
    for(gene in genes) {
      if(gene %in% genes_synonyms_only$Gene.synonyms.value) {
        found_dt[found_dt$Gene == gene, "Found"] <- "Synonym"
        found_dt[found_dt$Gene == gene, "Synonym"] <- genes_synonyms_only[genes_synonyms_only$Gene.synonyms.value == gene, "Gene.symbol"]
      } else if(gene %in% found_genes$Gene.symbol) {
        found_dt[found_dt$Gene == gene, "Found"] <- TRUE
      }
    }
    return(found_dt)
}

