# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get InterPro families for a single protein accession
#'
#' @description Calls the IntePro API with the following API request 'https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/accession' for single searchID.
#' Results from the API call are filtered to only include those with an InterPro type of 'Family', excluding super families and domains. A final comma separated list of family IDs is added to a single row dataframe with the protein accession search term and species.
#' 
#' @param accession The protein accession to search
#' @param species The species to search
#' @include global_functions.R
#' @return singe row data frame of protein.accession, species and comma delimited list of InterPro ids
runInterproSearch <- function(accession, species) {
  checkClass(parameter = "accession", value = accession, classType = "character", functionName = "runInterproSearch")
  checkSpecies(species = species, functionName = "runInterproSearch", includeHuman = TRUE)
  requestQuery <- paste0("https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/", accession)
  api_result <- callApi(query = requestQuery, failOn404 = TRUE, failonErrors = FALSE, simplifyVector = T)
  api_resultDT <- api_result$result
  family_only <- api_resultDT[api_resultDT$metadata$type == "family", "metadata"]
  protein_interpro <- data.frame("protein.accession" = accession,
                                 "Gene.homologues.homologue.organism.shortName" = species,
                                 "family" = paste0(family_only$accession, collapse = ", "))
  return(protein_interpro)
}

#' Filter a single list of InterPro family IDs based on the InterPro family hierarchy
#'
#' @description Filters a list of InterPro IDs based on their hierarchy position based on the interpro_levels data.frame.
#' A list of InterPro IDs are passed to the function as a single comma separated string and only IDs which are a level 2 (first child) for a branch in the hierarchy is kept. Each branch is filtered separately.
#' If a branch has only a level 1 ID (parent) and no child in the list of InterPro IDs, the level 1 id is kept.
#' If an InterPro ID is not found in the hierarchy (usually down to different API and hierarchy versions used) then the InterPro ID is kept also.
#' 
#' @param interpro_list A single comma separated string of InterPro IDs to filter. This should a list of family IDs belonging to a single protein.
#' @param interpro_levels The InterPro hierarchy data.frame including hierarchy level column and hierarchy group (branch group). 
#' @importFrom logger log_warn
#' @include global_functions.R
#' @return singe row data frame of a comma delimited list of unfiltered InterPro IDs, and comma delimited list of filtered InterPro ids.
filterInterProPathwayResults <- function(interpro_list, interpro_levels) {
  checkClass(parameter = "interpro_list", value = interpro_list, classType = "character", functionName = "filterInterProPathwayResults")
  checkClass(parameter = "interpro_levels", value = interpro_levels, classType = "data.frame", functionName = "filterInterProPathwayResults")
  columns_needed <-  c("interpro_id", "hierarchy_level", "family_branch")
  checkColumns(dt = interpro_levels, columns = columns_needed, dtParameter = "interpro_levels", functionName = "filterInterProPathwayResults")
  interpro_list_split <- unique(unlist(strsplit(interpro_list, split = ", ")))
  interpro_levels_filtered <- interpro_levels[interpro_levels$interpro_id %in% interpro_list_split, ]
  keep_list <- c()
  not_found <- c()
  missing_id <- interpro_list_split[!interpro_list_split %in% interpro_levels_filtered$interpro_id]
  not_found <- append(not_found, missing_id)
  if(nrow(interpro_levels_filtered) >0) {
    for(branch in unique(interpro_levels_filtered$`family_branch`)){
      interpro_levels_filtered_branch <- interpro_levels_filtered[interpro_levels_filtered$`family_branch` == branch, ]
      if(2 %in% interpro_levels_filtered_branch$`hierarchy_level`) {
        keep_list <- append(keep_list,
                            interpro_levels_filtered_branch[interpro_levels_filtered_branch$`hierarchy_level` == 2, "interpro_id"])
      } else {
        keep_list  <- append(keep_list,
                             interpro_levels_filtered_branch[interpro_levels_filtered_branch$`hierarchy_level` == 1, "interpro_id"])
      }
    }
  } else {
    not_found <- append(not_found, interpro_list_split)
    
  }
  combine_list <- unique(c(keep_list, not_found))
  combine_text <- paste0(combine_list, collapse = ", ")
  return_dt <- data.frame(family = interpro_list, filtered_family = combine_text)
  return(list("return_dt" = return_dt, "not_found" = not_found))
}

#' Filter InterPro ID results (all proteins) based on hierarchy position for a dataframe of protein families.
#'
#' @description Each protein family results are filtered to only IDs which are a level 2 (first child) for a branch in the hierarchy is kept. Each branch is filtered separately.
#' If a branch has only a level 1 ID (parent) and no child in the list of InterPro IDs, the level 1 id is kept. If an InterPro ID is not found in the hierarchy then the InterPro ID is kept also.
#' 
#' Each protein result is filtered separately. 
#' 
#' @details Filtering is done using the InterPro Hierarchy file found in the interpro_hierarchy_with_levels directory of the dataLocation directory.
#' The file which matches the InterPro API version shall be used or if not available the file will be updated using this version - see ?updateInterproHierarchy.
#' @include global_functions.R
#' @param dataLocation The directory path of where the 'g2pScanData' folder should be found
#' @param familyResults A data.frame of 'protein.accession', 'Gene.homologues.homologue.organism.shortName' and 'family'. 
#' @return data.frame of 'protein.accession', 'Gene.homologues.homologue.organism.shortName' and 'filtered_family'. 
filterInterProResults <- function(familyResults, dataLocation = ".") {
  checkClass(parameter = "familyResults", value = familyResults, classType = "data.frame", functionName = "filterInterProResults")
  columns_needed <-  c("protein.accession", "Gene.homologues.homologue.organism.shortName", "family")
  checkColumns(dt = familyResults, columns = columns_needed, dtParameter = "familyResults", functionName = "filterInterProResults")
  
  update_path <- updateInterproHierarchy(dataLocation = dataLocation) # set localVars interpro_levels to get path to right hierarchy csv - update if needed.
  interpro_levels <- read.csv(update_path, row.names = NULL)
  familys_filitered_output <- lapply(unique(familyResults$family), FUN = filterInterProPathwayResults, interpro_levels = interpro_levels)
  familys_filitered_dt <- do.call(rbind, lapply(familys_filitered_output, `[[`, 1))
  family_filtered_merge <- merge(familyResults, familys_filitered_dt, by = "family", all.x = T)
  family_filtered_dt <- family_filtered_merge[, c( "protein.accession","Gene.homologues.homologue.organism.shortName", "filtered_family")]
  not_found <- unique(unlist(lapply(familys_filitered_output, `[[`, 2)))
  return(list("family_filtered_dt" = family_filtered_dt, "not_found" = not_found))
}


#' Get families for each protein accession retrieved from getProteins.
#'
#' @description For every protein row passed to this function in the proteinResult, InterPro's API will be called to pull InterPro family IDs for protein accession. 
#' Steps:
#' \itemize{
#' \item{1. Loop through each non-NA row of the input proteinResult and take protein.accession as the searchID and "Gene.homologues.homologue.organism.shortName" as species of the input to the runInterproSearch function}
#' \item{2. The runInterproSearch function to run the following API request 'https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/accession' for single searchID.}
#' \item{3. Results from the API call are filtered to only include those with an InterPro type of 'Family', excluding superfamilies and domains.}
#' \item{4. Results from each API call are combine into one large dataframe}
#' \item{5. Each protein family results are filtered to only IDs which are a level 2 (first child) for a branch in the hierarchy is kept. Each branch is filtered separately.
#'  If a branch has only a level 1 ID (parent) and no child in the list of InterPro IDs, the level 1 id is kept. If an InterPro ID is not found in the hierarchy then the InterPro id is kept also.
#'  Each protein result is filtered separately.}}
#' 
#' If cluster is passed to the function the API calls are run in parallel
#' 
#' @inheritParams getProteinMatrix
#' @inheritParams runAllProteinSearch
#' @inheritParams filterInterProResults
#' @importFrom foreach foreach %dopar%
#' @importFrom logger log_info
#' @importFrom parallel stopCluster
#' @include global_functions.R
#' @return data.frame of 'protein.accession', 'Gene.homologues.homologue.organism.shortName' and final filtered 'family'. 
#' @export
getFamilies <- function(proteinResult, cluster = NULL, stopCluster = TRUE,  dataLocation = ".") { 
  checkClass(parameter = "proteinResult", value = proteinResult, classType = "data.frame", functionName = "getFamilies")
  columns_needed <- c("protein.accession",
                      "Gene.symbol",
                      "Gene.homologues.homologue.symbol",
                      "Gene.homologues.homologue.organism.shortName"
  )
  checkColumns(dt = proteinResult, columns = columns_needed, dtParameter = "proteinResult", functionName = "getFamilies")
  
  proteinResult_noNA <- proteinResult[!is.na(proteinResult$protein.accession), ] 
  search_terms <- proteinResult_noNA[,c("protein.accession", "Gene.homologues.homologue.organism.shortName")]
  logger::log_info("Running InterPro protein API for {nrow(search_terms)} to get InterPro families.")
  if(!is.null(cluster)){
    protein_families <- foreach::foreach(X = 1:nrow(search_terms), .combine = rbind) %dopar% {
      accession <- search_terms[X, "protein.accession"]
      sp <- search_terms[X, "Gene.homologues.homologue.organism.shortName"]
      families <- runInterproSearch(accession = accession, species = sp)
      return(families)
    }
    if(stopCluster) {
      endCluster()
    }
  } else {
    protein_families_list <- lapply(X = 1:nrow(search_terms), FUN = function(X) {
      accession <- search_terms[X, "protein.accession"]
      sp <- search_terms[X, "Gene.homologues.homologue.organism.shortName"]
      families <- runInterproSearch(accession = accession, species = sp)
      return(families)
    })
    protein_families <- do.call("rbind", protein_families_list)
  }
  filtered_families_output <- filterInterProResults(protein_families, dataLocation = dataLocation)
  filtered_families <- filtered_families_output$family_filtered_dt
  not_found <-  filtered_families_output$not_found
  if(length(not_found) > 0) {
    logger::log_warn("{paste0(not_found, collapse = ', ')} family IDs could not be found in the provided family hierachy. So each of these IDs are included in the family counts as hierarchy level is not known. This can because the family is not apart of a hierarchy.")
  }
  return(filtered_families)
}

#' Get matrix of proteins
#'
#' @description Create data.frame matrix of human genes as rows and columns as species, including human, where families are listed as values.
#' @param familyResult Data from output of function getFamilies
#' @inheritParams getProteinMatrix
#' @importFrom logger log_fatal
#' @importFrom R.utils getOption
#' @include global_functions.R
#' @return data.frame matrix of human genes as rows and columns as species,including Human, where families are listed as values.
#' @export
getFamilyMatrix <- function(familyResult, proteinResult, species) {
  checkClass(parameter = "familyResult", value = familyResult, classType = "data.frame", functionName = "getFamilyMatrix")
  checkClass(parameter = "proteinResult", value = proteinResult, classType = "data.frame", functionName = "getFamilyMatrix")
 
  
  if(!"Human" %in% proteinResult$Gene.homologues.homologue.organism.shortName) {
    logger::log_fatal("No Human data found in proteinResult dataframe passed to getFamilyMatrix")
    stop("No Human data found in proteinResult dataframe passed to getFamilyMatrix")
  }
  
  columns_needed <- c("filtered_family",
                      "protein.accession")
  missing <- columns_needed[!columns_needed %in% colnames(familyResult)]
  if(length(missing) > 0){
    logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter proteinResult passed to getFamilyMatrix")
    stop("Missing columns in dataframe passed to getFamilyMatrix")
  }
  columns_needed <- c("Gene.symbol",
                      "protein.accession")
  missing <- columns_needed[!columns_needed %in% colnames(proteinResult)]
  if(length(missing) > 0){
    logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter proteinResult passed to getFamilyMatrix")
    stop("Missing columns in dataframe passed to getFamilyMatrix")
  }
  
  protien_family_merged <- merge(x = proteinResult, y = familyResult[, c("protein.accession",
                                                                          "filtered_family")],
                                by = c("protein.accession"),
                                all.x = TRUE)
  
  human_genes <- unique(protien_family_merged$Gene.symbol)
  
  proteins_sorted <- protien_family_merged[order(proteinResult$Gene.symbol), ] 
  # get species found in human_genes Get full list of species possible to maintain order of species.
  species_ordered <- c("Human", c("Human", checkSpecies(species = species, functionName = "getFamilyMatrix", includeHuman = F)))
  
  # create empty matrix
  family_matrix <- data.frame(matrix(nrow = length(human_genes), ncol = length(species_ordered)))
  rownames(family_matrix) <- sort(human_genes)
  colnames(family_matrix) <- species_ordered
  # fill matrix
  for(sp in species_ordered) {
    species_results <- protien_family_merged[protien_family_merged$Gene.homologues.homologue.organism.shortName == sp, ]
    # concat when there are more than one protein for 1 human gene
    species_results[,"families.combined.string"] <- sapply(X= species_results$Gene.symbol,
                                                           FUN = function(p, species_results){
                                                             families <- species_results[species_results$Gene.symbol == p, "filtered_family"]
                                                             families_vector <- unlist(strsplit(families, split = ", "))
                                                             if(length(families_vector)>0) {
                                                               output <- paste0(unique(families_vector), collapse = ", ")
                                                             } else {
                                                               output <- NA
                                                             }
                                                             output
                                                           }, species_results)
    species_human_genes <- species_results$Gene.symbol
    family_matrix[species_human_genes, sp] <- species_results$families.combined.string  }
  return(family_matrix)
}
