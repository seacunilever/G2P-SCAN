# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Check data is expected orthologue data
#'
#' @description This function checks data is a data.frame of following columns \itemize{
#' \item {Gene.primaryIdentifier}
#' \item {Gene.symbol}
#' \item {Gene.homologues.homologue.primaryIdentifier}
#' \item {Gene.homologues.homologue.symbol}
#' \item {Gene.homologues.homologue.organism.shortName}
#' \item {Gene.homologues.type}
#' }
#' @param orthologueData Data to be checked for class and columns.
#' @include global_functions.R
checkOrthologueDataFrame <- function(orthologueData, functionName) {
  checkClass(parameter = "orthologueData", value = orthologueData, classType = "data.frame", functionName = functionName)
  needed_cols <- c("Gene.primaryIdentifier",
                   "Gene.symbol",
                   "Gene.homologues.homologue.primaryIdentifier",
                   "Gene.homologues.homologue.symbol",
                   "Gene.homologues.homologue.organism.shortName",
                   "Gene.homologues.type")
  checkColumns(dt = orthologueData, columns = needed_cols, dtParameter = orthologueData, functionName = functionName)
}

#' Filter orthologues for species of interest
#'
#' @description Reads in data.frame of orthologues and filters the "Gene.homologues.homologue.organism.shortName" column on species input.
#' @param species List of species to select to run analysis on and filter for, out of: \itemize{
#' \item {R. norvegicus}
#' \item {M. musculus}
#' \item {D. rerio}
#' \item {C. elegans}
#' \item {D. melanogaster}
#' \item {S. cerevisiae}
#' } 
#' When NULL no filter is applied and all species are returned
#' @param dataToFilter Orthologue data.frame to filter. Columns should include all of the following \itemize{
#' \item {Gene.primaryIdentifier}
#' \item {Gene.symbol}
#' \item {Gene.homologues.homologue.primaryIdentifier}
#' \item {Gene.homologues.homologue.symbol}
#' \item {Gene.homologues.homologue.organism.shortName}
#' \item {Gene.homologues.type}
#' }
#' 
#' @importFrom logger log_fatal
#' @include global_functions.R
#' 
#' @return Filtered orthologue data.frame
filterOrthologueSpecies <- function(dataToFilter, species = NULL) {
  checkOrthologueDataFrame(dataToFilter, functionName = "filterOrthologueSpecies")
  species <- checkSpecies(species = species, functionName = "filterOrthologueSpecies")
  data_filtered <- dataToFilter[dataToFilter$"Gene.homologues.homologue.organism.shortName" %in% species, ]
  if(nrow(data_filtered) == 0) {
    logger::log_fatal("No orthologue results returned for the species input.")
    stop("No orthologue results returned.")
  }
  return(data_filtered)
}

#' Filter orthologues for least divergent orthologues only
#'
#' @description Reads in data.frame of orthologues and filters the "Gene.homologues.type" column to only include "least divergent orthologues".
#' @inheritParams filterOrthologueSpecies 
#' @importFrom logger log_fatal
#' 
#' @return Filtered orthologue data.frame
filterOrthologueLDO <- function(dataToFilter, species = NULL) {
  checkOrthologueDataFrame(dataToFilter)
  
  data_filtered <- dataToFilter[dataToFilter$Gene.homologues.type == "least diverged orthologue", ]
  if(nrow(data_filtered) == 0) {
    logger::log_fatal("No orthologue results returned for the species with orthologue filter.")
    stop("No orthologue results returned.")
  }
  return(data_filtered)
}

#' Get all orthologues
#'
#' @description This function get all orthologues for a list of genes. Species to get orthologues for can be chosen from \itemize{
#' \item{R. norvegicus}
#' \item{M. musculus}
#' \item{D. rerio}
#' \item{C. elegans}
#' \item{D. melanogaster}
#' \item{S. cerevisiae}}
#' Orthologues are pulled from Panther via HumanMine using HumanMine query template "Gene_Orth". Orthologues which are included can be chosen from the following: \itemize{
#' \item {All Orthologues - where all orthologues are selected}
#' \item {Least Divergent Orthologues (LDO) - where only the most nearly 'equivalent' gene in another organism is selected (in case of gene duplications for instance)}
#' }
#' 
#' @param genes List of Human gene symbol of which orthologues should be searched
#' @inheritParams filterOrthologueSpecies
#' @param orthologueFilter "All" (include all orthologues) or "LDO" (filter to least divergent orthologues)
#'
#' @importFrom logger log_fatal
#' @include global_functions.R
#' @return Data.frame of Human genes and species orthologues. Columns :\itemize{
#' \item {Gene.primaryIdentifier}
#' \item {Gene.symbol}
#' \item {Gene.homologues.homologue.primaryIdentifier}
#' \item {Gene.homologues.homologue.symbol}
#' \item {Gene.homologues.homologue.organism.shortName}
#' \item {Gene.homologues.type}
#' }
#' @export
getOrthologueGenes <- function(genes, species = NULL, orthologueFilter = "LDO") {
  checkClass(parameter = "genes", value = genes, classType = "character", functionName = "getOrthologueGenes")
  checkClass(parameter = "orthologueFilter", value = orthologueFilter, classType = "character", functionName = "getOrthologueGenes")
  if(!any(tolower(orthologueFilter) %in% c("ldo", "all"))) {
    logger::log_fatal("OrthologueFilter passed to getOrthologueGenes is neither 'LDO' or 'All'. Value passed: {orthologueFilter}.")
    stop("OrthologueFilter passed to getOrthologueGenes is neither 'LDO' or 'All'")
  }
  
  genes <- unique(genes)
  
  gene_to_ortho_result <- unique( #unique as each gene is run for upper and lower case and likely to bring back same result
    runHumanMineQuery(templateName = "Gene_Orth",
                      constraintPath = "Gene.symbol",
                      constraintOperators = "=",
                      constraintValues = list(genes),
                      constraintm.index = 1,
                      select = c("Gene.primaryIdentifier",
                                 "Gene.symbol",
                                 "Gene.homologues.homologue.primaryIdentifier",
                                 "Gene.homologues.homologue.symbol",
                                 "Gene.homologues.homologue.organism.shortName",
                                 "Gene.homologues.type")
    )
  )
  if(!is.null(gene_to_ortho_result)){
    # Filter results for only orthologues of selected organisms
    gene_to_ortho_species_filtered <- filterOrthologueSpecies(dataToFilter = gene_to_ortho_result, species = species)
    # Filter for orthologues is needed.
    if(tolower(orthologueFilter) == "ldo") {
      gene_to_ortho_species_ortho_filtered <- filterOrthologueLDO(dataToFilter = gene_to_ortho_species_filtered)
    } else {
      gene_to_ortho_species_ortho_filtered <- gene_to_ortho_species_filtered
    }
    # Set 'Gene symbol not found' when no orthologue is found for a gene-species.
    gene_to_ortho_species_ortho_filtered[(is.na(gene_to_ortho_species_ortho_filtered[, "Gene.homologues.homologue.symbol"]) |
                                            trimws(gene_to_ortho_species_ortho_filtered[, "Gene.homologues.homologue.symbol"]) == ""),
                                         "Gene.homologues.homologue.symbol"] <-"Gene symbol not found"
    return(gene_to_ortho_species_ortho_filtered)
  } else {
    logger::log_fatal("No orthologue results returned for gene input.")
    stop("No orthologue results returned.")
  }
}


#' Format orthologues for one species and one human gene into single string
#'
#' @description Create a single orthologue string for use in output representation using the following orthologue information: HumanMine primaryidentifier, gene symbol or both.
#' This function takes the output of the function 'getOrthologues' and select the singular species orthologue for the given human genes and returns a formatted string with select information type.
#' @param orthologueData The data.frame from the output of function getOrthologueGenes
#' @param humanGene Human gene to get orthologue of
#' @param species Species to the the orthologue gene for
#' @param orthologueOutput Whether orthologue matrix should use the gene HumanMine "primaryIdentifier" or "geneSymbol" or "both"(default) (not case sensitive). When "both" an orthologue will be outputed as 'primaryIdentifier(geneSymbol)'. Orthologues regularly do not have gene symbols associated with them. 
#' 
#' @importFrom logger log_warn log_fatal
#' @importFrom R.utils getOption
#' @include global_functions.R
#' @return Formated orthologue for single human gene and species
getFormatedOrthologueString <- function(orthologueData, humanGene, species, orthologueOutput) {
  checkClass(parameter = "humanGene", value = humanGene, classType = "character", functionName = "getFormatedOrthologueString")
  checkOrthologueDataFrame(orthologueData, functionName = "getFormatedOrthologueString")
  checkSpecies(species = species, functionName = "getFormatedOrthologueString")
  
  checkClass(parameter = "orthologueOutput", value = orthologueOutput, classType = "character", functionName = "getFormatOrthologue")
  if(!any(tolower(orthologueOutput) %in% c("primaryidentifier", "genesymbol", "both"))) {
    logger::log_fatal("orthologueOutput passed to getFormatOrthologue is neither 'primaryIdentifier', 'geneSymbol' or 'both'. Value passed: {orthologueOutput}.")
    stop("orthologueOutput passed to getFormatOrthologue is neither 'primaryIdentifier', 'geneSymbol' or 'both'")
  }
  
  gene_species_orthologues_DT <- orthologueData[orthologueData$Gene.symbol == humanGene & orthologueData$Gene.homologues.homologue.organism.shortName == species, ]
  if(nrow(gene_species_orthologues_DT) > 0) {
    if(tolower(orthologueOutput) == "both") {
      output <- paste0(gene_species_orthologues_DT$Gene.homologues.homologue.primaryIdentifier, "(", gene_species_orthologues_DT$Gene.homologues.homologue.symbol, ")")
    }
    else if(tolower(orthologueOutput) == "primaryidentifier") {
      output <- gene_species_orthologues_DT$Gene.homologues.homologue.primaryIdentifier
    }
    else {
      output <- gene_species_orthologues_DT$Gene.homologues.homologue.symbol
    }
  } else {
    output <- NA
  }
  output_string <- paste0(output, collapse = ", ")
  return(output_string)
}

#' Get matrix of orthologues
#'
#' @description Create data.frame matrix of human genes as rows and columns as species where orthologues are listed as values.
#' @param orthologueData Data from output of function getOrthologueGenes
#' @param allHumanGenes Character vector of all human genes symbols to be included in matrix, even if there is no orthologue result in orthologueData. If NULL only genes found in orthologueData will be used.
#' @inheritParams getFormatedOrthologueString
#' @inheritParams getOrthologueGenes
#' @importFrom logger log_warn log_fatal
#' @importFrom R.utils getOption
#' @include global_variables.R
#' @return data.frame matrix of human genes as rows and columns as species where orthologues are listed as values.
#' @export
getOrthologueMatrix <- function(orthologueData, species, allHumanGenes = NULL, orthologueOutput = "both") {
  checkOrthologueDataFrame(orthologueData, functionName = "getOrthologueMatrix")
  orthologue_data_human_genes <- unique(orthologueData$Gene.symbol)
  
  if(!is.null(allHumanGenes)) {
    checkClass(parameter = "allHumanGenes", value = allHumanGenes, classType = "character", functionName = "getOrthologueMatrix")
  } else {
    allHumanGenes <- orthologue_data_human_genes
  }
  
  orthologue_not_all_human <- orthologue_data_human_genes[!orthologue_data_human_genes %in% allHumanGenes] #Human genes found in orthologue data but not in allHumanGenes vectors, suggesting allHumanGenes was not the same input for getOrthologueGenes function to produce orthologueData.
  if(length(orthologue_not_all_human) > 0) {
    logger::log_warn("Some human gene symbols found in orthologueData are not found in allHumanGenes list in function getOrthologueMatrix. This suggests that allHumanGenes is not the same input used to create orthologueData via function getOrthologueGenes(). Not found genes: {orthologue_not_all_human}.
                     Matrix output of getOrthologueMatrix will include all genes found in orthologueData and allHumanGenes.")
  }
  
  genes_no_orthologues <- unique(allHumanGenes[!allHumanGenes %in% orthologue_data_human_genes]) # genes not found in orthologueData and so have no orthologues.
  all_genes <- unique(c(orthologue_data_human_genes, allHumanGenes))
  
  # get species found in orthologueData. Get full list of species possible to maintain order of species.
  species_ordered <-  checkSpecies(species = species, functionName = "getOrthologueMatrix", includeHuman = F)
  
  # create empty matrix
  orthologue_matrix <- data.frame(matrix(nrow = length(all_genes), ncol = length(species_ordered)))
  rownames(orthologue_matrix) <- sort(all_genes)
  colnames(orthologue_matrix) <- species_ordered
  
  for(sp in species_ordered) {
    orthologue_matrix[all_genes, sp] <- sapply(FUN = getFormatedOrthologueString, X= all_genes, species = sp, orthologueData = orthologueData, orthologueOutput = orthologueOutput)
  }
  orthologue_matrix[orthologue_matrix == "NA"] <- ""   
  return(orthologue_matrix)
}


