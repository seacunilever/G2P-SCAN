# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Count unique terms in matrix columns
#'
#' @description Count the number of unique values per column (species) in a matrix for a given list of Human genes. 
#' @details  Data to be counted will be filtered to the selected list human genes by filtering on row names. Unique values in each columns are then counted for the remaining rows and
#' the count values are returned as a list in the order of the columns and do not include counts for NA and blank values.
#' 
#' This function has been built so that the data to be counted is a matrix of orthologues or proteins across species (as columns) with the human gene mapped to the orthologues or proteins as row names.
#' It has been designed so that the genes to filter to data prior counting are a list of genes within a pathway.
#' This is so that counts for each species' genes/orthologues are returned for a singular pathway.
#' @param genesToCount A vector of human genes which should be used to filter matrix prior count
#' @param dataToCount Matrix of human gene rows and with columns indicating a new entity to counts eg, species. Data can include genes or proteins.
#' @return list of counts
#' @include global_functions.R
getDataMatrixCount <- function(genesToCount, dataToCount) {
  checkClass(parameter = "dataToCount", value = dataToCount, classType = "data.frame", functionName = "getDataMatrixCount")
  checkClass(parameter = "genesToCount", value = genesToCount, classType = "character", functionName = "getDataMatrixCount")
  path_data <- dataToCount[genesToCount, ]
  counts <- lapply(X = path_data, FUN = function(data){
    if(!all(is.na(data))) {
      separated <- strsplit(data, split = ", ")
      unique_values <- unique(unlist(separated))
      unique_values <- unique_values[unique_values != "NA" & !is.na(unique_values)]
      count <- length(unique_values)
    } else {
      count <- 0
    }
    count
  })
  return(counts)
}

#' Get counts of unique Human genes which have values in each of the species columns and identify if input genes are mapped.
#'
#' @description Count the number of Human Genes (within the given list of Human genes) which have values associated in each column (species) in a matrix.
#' @details  Data to be counted will be filtered to the selected list human genes by filtering on row names. The number of human genes with values in each column are counted for the remaining rows and
#' the count values are returned as a list in the order of the columns. The input genes will be counted per column and the found input genes and count will be returned per column.
#' 
#' This function has been built so that the data to be counted is a matrix of orthologues or proteins across species (as columns) with the human gene mapped to the orthologues or proteins as row names.
#' It has been designed so that the genes to filter to data prior counting are a list of genes within a pathway.
#' This is so that counts for each species' mapped Human genes are returned for a singular pathway.
#' @param inputGenes The human input genes used to conduct the analysis
#' @inheritParams getDataMatrixCount
#' @include global_functions.R
#' @return list of lists of mapped counts, the input genes found and count of input found. Each list is return in the same order as the matrix columns.
getDataMatrixUniqueHumanCount <- function(inputGenes, genesToCount, dataToCount) {
  checkClass(parameter = "dataToCount", value = dataToCount, classType = "data.frame", functionName = "getDataMatrixUniqueHumanCount")
  checkClass(parameter = "genesToCount", value = genesToCount, classType = "character", functionName = "getDataMatrixUniqueHumanCount")
  checkClass(parameter = "inputGenes", value = inputGenes, classType = "character", functionName = "getDataMatrixUniqueHumanCount")
  path_data <- dataToCount[genesToCount, ]
  mapped_list <- lapply(X = 1:ncol(path_data), FUN = function(col, data){
    data_col <-data.frame(orthos = data[, col], Human_genes = rownames(data))
    unique_human <- data_col[data_col$orthos != "NA" & !is.na(data_col$orthos) ,]
    mapped <- unique_human$Human_genes
    mapped
  }, data = path_data)
  mapped_count <- lapply(X = mapped_list, FUN = length)
  input_found <- lapply(X = mapped_list, inputGenes = inputGenes, FUN = function(x, inputGenes) {inputGenes[inputGenes %in% x]})
  input_found_count <- lapply(X = input_found, FUN = length)
  
  names(mapped_count) <- colnames(path_data)
  names(input_found) <- colnames(path_data)
  names(input_found_count) <- colnames(path_data)
  
  return(list(mapped_count, input_found, input_found_count))
  
}


#' Get counts of mapped and unmapped  families for a set of human genes.
#'
#' @description Count the number of unique terms (families) per column (species) and count the number of proteins which do not have terms associated with them
#' @details  Filter a data matrix of Human gene rows and data types (species) as column for a selection of human genes (pathway genes) and count number of unique ids in each of the columns. Values in columns are separated by commas before counting.
#' Mapped counts are the number of unique id found per column (species) for the genesToCounts. Unmapped counts are the counts of unique proteins without a family ID assigned.
#' Counts are returned as a list and do not include counts for NA and blank values.
#'
#' This function has been built so that the data to be counted is a matrix of families across species (as columns) with the human gene mapped to the protein of which the family is derived from are as row names.
#' It has been designed so that the genes to filter to data prior counting are a list of genes within a pathway.
#' This is so that counts for each species' families are returned for a singular pathway. 
#' @inheritParams getDataMatrixCount
#' @inheritParams getOrthologueGenes
#' @param familyDT The dataframe output of getFamilies function. Columns required: 'filtered_family', 'protein.accession' and 'Gene.homologues.homologue.organism.shortName'
#' @param proteinDT The dataframe output of the getProteins function. Columns required: 'Gene.symbol', 'protein.accession' and 'Gene.homologues.homologue.organism.shortName'
#' @include global_functions.R
#' @return list of mapped and unmapped counts
getFamilyCount <- function(genesToCount, familyDT, proteinDT, species) {
  checkClass(parameter = "proteinDT", value = proteinDT, classType = "data.frame", functionName = "getFamilyCount")
  checkClass(parameter = "genesToCount", value = genesToCount, classType = "character", functionName = "getFamilyCount")
  checkClass(parameter = "familyDT", value = familyDT, classType = "data.frame", functionName = "getFamilyCount")
  checkColumns(dt = familyDT, columns = c("filtered_family", "protein.accession", "Gene.homologues.homologue.organism.shortName") , dtParameter = "familyDT", functionName = "getFamilyCount")
  checkColumns(dt = proteinDT, columns = c("Gene.symbol", "protein.accession", "Gene.homologues.homologue.organism.shortName"), dtParameter = "proteinDT", functionName = "getFamilyCount")
  path_protein_data <- proteinDT[proteinDT$Gene.symbol %in% genesToCount, ]
  path_family <- merge(x = familyDT, path_protein_data, all.x = F, all.y = T, by = c("protein.accession", "Gene.homologues.homologue.organism.shortName"))
  species_vector <- c("Human", checkSpecies(species = species, functionName = "getFamilyCount", includeHuman = F))
  #species_vector <- unique(proteinDT$Gene.homologues.homologue.organism.shortName)
  counts <- lapply(X = species_vector, FUN = function(sp, path_data){
    data_sp <- path_data[path_data$Gene.homologues.homologue.organism.shortName == sp, ]
    separated <- strsplit(data_sp$filtered_family, split = ", ")
    unique_values <- unique(unlist(separated))
    unique_values <- unique_values[!is.na(unique_values) & unique_values != ""]
    count <- length(unique_values)
    count
  }, path_data = path_family)
  unmapped <- lapply(X = species_vector, FUN = function(sp, path_data){
    data_sp <- path_data[path_data$Gene.homologues.homologue.organism.shortName == sp, ]
    na_values <- unique(data_sp[(is.na(data_sp$filtered_family) | data_sp$filtered_family == "") & !is.na(data_sp$protein.accession), "protein.accession"])
    count <- length(na_values)
    count
  }, path_data = path_family)
  names(counts) <- species_vector
  names(unmapped) <- species_vector
  return(list(counts, unmapped))
}

#' Get counts for given pathway
#'
#' @description For one given pathway, the following counts are calculated:
#' \itemize{
#' \item{Pathway Reactome ID - The pathway identifier for the data in the dataframe}
#' \item{Pathway Reactome Name - The pathway name for the data in the dataframe}
#' \item{Input Found - The input genes which have been found in that Human pathway or have a orthologue found for in other species}
#' \item{Input Found Count - The number of input genes which have been found in that Human pathway or have a orthologue found for in other species}
#' \item{Mapped Human Genes - The number of Human genes which are represented in the orthologues found for each of the species}
#' \item{Total Gene Count - The total number of unique gene/ orthologues found for the given pathway and species}
#' \item{Coverage Percentage - The percentage of input genes found in the pathway vs the number of genes in the pathway. Value only given for Human data}
#' \item{Protein Count - The number of unique proteins found in each species for the given pathway}
#' \item{Family Count - The number of unique families found in each species for the given pathway}
#' \item{Proteins Unmapped to families Count - The number of unique proteins without a family ID assigned in each species for the given pathway}
#' \item{Entity Count - The number of reactome entities associated with the given pathway for each species}
#' \item{Reaction Count - The number of reactome reactions associated with the given pathway for each species}}
#' 
#' @param pathway The Reactome pathway ID of which counts should be calculated
#' @param pathwayDT The output of the getPathways function.
#' @param inputGenes The human input genes used to conduct the analysis
#' @param genesInPath The output of the getGenes.
#' @param orthologueMatrix The output of getOrthologueMatrix. Recommended to use matrix of primary ids.
#' @param proteinDT The output of getProteins
#' @param proteinMatrix The output of getProteinMatrix
#' @param familyDT The output of getFamilies
#' @param entitiesReactions The output of getEntitiesReactionsCounts
#' @inheritParams getOrthologueGenes
#' @include global_functions.R
#' @return a dataframe of count types as rows, the species as columns and count as the data. 
getCountsForPathway <- function(pathway, species, pathwayDT, inputGenes, genesInPath, orthologueMatrix, proteinDT, proteinMatrix, familyDT, entitiesReactions) {
  checkClass(parameter = "pathway", value = pathway, classType = "character", functionName = "getCountsForPathway")
  checkClass(parameter = "inputGenes", value = inputGenes, classType = "character", functionName = "getCountsForPathway")
  checkClass(parameter = "genesInPath", value = genesInPath, classType = "data.frame", functionName = "getCountsForPathway")
  checkClass(parameter = "orthologueMatrix", value = orthologueMatrix, classType = "data.frame", functionName = "getCountsForPathway")
  checkClass(parameter = "proteinDT", value = proteinDT, classType = "data.frame", functionName = "getCountsForPathway")
  checkClass(parameter = "proteinMatrix", value = proteinMatrix, classType = "data.frame", functionName = "getCountsForPathway")
  checkClass(parameter = "familyDT", value = familyDT, classType = "data.frame", functionName = "getCountsForPathway")
  checkClass(parameter = "entitiesReactions", value = entitiesReactions, classType = "data.frame", functionName = "getCountsForPathway")
  checkColumns(dt = pathwayDT, columns = c("Pathway Identifier", "Pathways"), dtParameter = "pathwayDT", functionName = "getCountsForPathway")
  checkColumns(dt = genesInPath, columns = c("Pathway Identifier", "Gene Symbol"), dtParameter = "genesInPath", functionName = "getCountsForPathway")
  
  ordered_species <-  c("Human", checkSpecies(species = species, functionName = "getCountsForPathway", includeHuman = F))
  pathway_name <- pathwayDT[pathwayDT$`Pathway Identifier` == pathway, "Pathways"]
  pathway_genes <- unique(genesInPath[genesInPath$`Pathway Identifier` == pathway, "Gene Symbol"])
  input_mapped <- inputGenes[tolower(inputGenes) %in% tolower(pathway_genes)]
  n_input_mapped <- length(input_mapped)
  percent_input_mapped <- round(n_input_mapped/length(pathway_genes) * 100, digits = 2)
  orthologueMatrix$Human <- rownames(orthologueMatrix)
  mapped_human_genes <- getDataMatrixUniqueHumanCount(genesToCount = pathway_genes, inputGenes = inputGenes, dataToCount = orthologueMatrix)
  mapped_human_genes_count <- mapped_human_genes[[1]]
  input_found <- mapped_human_genes[[2]]
  input_found_count <- mapped_human_genes[[3]]
  gene_count <- getDataMatrixCount(genesToCount = pathway_genes, dataToCount = orthologueMatrix)
  protein_count <- getDataMatrixCount(genesToCount = pathway_genes, dataToCount = proteinMatrix)
  family_count_fun <- getFamilyCount(genesToCount = pathway_genes, familyDT = familyDT, proteinDT = proteinDT, species = species)
  family_count <- family_count_fun[[1]]
  family_proteins_unmapped <- family_count_fun[[2]]
  entities_counts <- as.list(entitiesReactions[pathway, grepl("_entities", colnames(entitiesReactions))])
  names(entities_counts) <- gsub("_entities", "", names(entities_counts))
  reactions_counts <- as.list(entitiesReactions[pathway, grepl("_reactions", colnames(entitiesReactions))])
  names(reactions_counts) <- gsub("_reactions", "", names(reactions_counts))
  counts_DT <- (do.call(rbind, list("Pathway Reactome ID"= pathway,
                                    "Pathway Name" = pathway_name,
                                    "Input Found" = input_found[ordered_species],
                                    "Input Found Count" = input_found_count[ordered_species],
                                    "Mapped Human Genes" = mapped_human_genes_count[ordered_species],
                                    "Total Gene Count" = gene_count[ordered_species],
                                    "Coverage %" = NA,
                                    "Protein Count" = protein_count[ordered_species],
                                    "Family Count" = family_count[ordered_species],
                                    "Proteins Unmapped to families Count" = family_proteins_unmapped[ordered_species],
                                    "Entity Count" = entities_counts[ordered_species],
                                    "Reaction Count" = reactions_counts[ordered_species])))

  colnames(counts_DT) <- ordered_species
  counts_DT["Coverage %", "Human"] <- percent_input_mapped
  input_null <- sapply(counts_DT["Input Found",], is.null)
  counts_DT["Input Found", input_null] <- NA
  i1 <- sapply(counts_DT, is.null) 
  counts_DT[i1] <- 0
  return(counts_DT)
}

#' Get counts for all pathways
#'
#' @description Loop through all pathways to calculate the following counts:
#' \itemize{
#' \item{Pathway Reactome ID - The pathway identifier for the data in the dataframe}
#' \item{Pathway Reactome Name - The pathway name for the data in the dataframe}
#' \item{Input Found - The input genes which have been found in that Human pathway or have a orthologue found for in other species}
#' \item{Input Found Count - The number of input genes which have been found in that Human pathway or have a orthologue found for in other species}
#' \item{Mapped Human Genes - The number of Human genes which are represented in the orthologues found for each of the species}
#' \item{Total Gene Count - The total number of unique gene/ orthologues found for the given pathway and species}
#' \item{Coverage Percentage - The percentage of input genes found in the pathway vs the number of genes in the pathway. Value only given for Human data}
#' \item{Protein Count - The number of unique proteins found in each species for the given pathway}
#' \item{Family Count - The number of unique families found in each species for the given pathway}
#' \item{Proteins Unmapped to families Count - The number of unique proteins without a family ID assigned in each species for the given pathway}
#' \item{Entity Count - The number of reactome entities associated with the given pathway for each species}
#' \item{Reaction Count - The number of reactome reactions associated with the given pathway for each species}}
#' 
#' @inheritParams getCountsForPathway
#' @inheritParams getOrthologueGenes
#' @include global_functions.R
#' @return list of dataframes of count types as rows, the species as columns and count as the data. Dataframes in the list are named with the pathway ID.
#' @export
getAllCounts <- function(inputGenes, species, pathwayDT, genesInPath, orthologueMatrix, proteinDT, proteinMatrix, familyDT, entitiesReactions) {
  species_list <-  c("Human", checkSpecies(species = species, functionName = "getAllCounts", includeHuman = F))
  
  #species_list <- colnames(proteinMatrix)
  pathway_ids <- unique(pathwayDT$`Pathway Identifier`)
  
  allPathwaysCounts <- lapply(X = pathway_ids, FUN = getCountsForPathway,
                              species = species,
                              inputGenes = inputGenes,
                              pathwayDT = pathwayDT,
                              genesInPath = genesInPath,
                              orthologueMatrix = orthologueMatrix,
                              proteinDT = proteinDT,
                              proteinMatrix = proteinMatrix,
                              familyDT = familyDT,
                              entitiesReactions = entitiesReactions)
  names(allPathwaysCounts) <- pathway_ids
  return(allPathwaysCounts)
}

#' Get counts for all pathways
#'
#' @description Add worksheets to workbook object where the pathway count matrix is given on a separate tab.
#' @param wb workbook to add worksheet to
#' @param pathwaysCountsList List of pathway count data.frames
#' @importFrom openxlsx writeDataTable addWorksheet
#' @return wb object with a sheet for each data.frame in pathwayCountsList
#' @export
getAllPathwayCountsWb <- function(wb, pathwaysCountsList) {
  checkClass(parameter = "pathwaysCountsList", value = pathwaysCountsList, class = "list", functionName = "getAllPathwayCountsWb")
  pathways <- names(pathwaysCountsList)
  for(path in pathways) {
    path_sheet <- openxlsx::addWorksheet(wb = wb, sheetName = path)
    counts <- as.data.frame(pathwaysCountsList[[path]])
    openxlsx::writeDataTable(wb = wb, path_sheet, x = counts, rowNames = T, na.string = "" )
  }
  return(wb)
}
