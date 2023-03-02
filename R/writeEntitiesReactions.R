# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get Human Entities Counts
#'
#' @description Run Reactome content service APIs to extract Human pathway entity counts for a single given pathway identifier.
#' The API will pulls the data of the live version of Reactome and so this function only needs to be run when the HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
#' The function calls "https://reactome.org/ContentService/data/participants/path_id/referenceEntities" and filters participants returned to only those from the following databases:
#' \itemize{
#' \item{UniProt}
#' \item{ChEBI}
#' \item{EMBL}
#' \item{ENSEMBL}
#' \item{miRBase}}
#'  Finally the number of participants returned after filtering is set as the Human entity count for the given pathway.
#'  This filtering occurs to match those counts produced in the Reactome species comparison tool.
#'  This function is intended to be run for pathways which no count data is found across other species when updating the entities and reaction file used in the tool.
#'  @param pathway Reactome pathway identifier of interested
#'  @importFrom logger log_warn
#'  @include global_functions.R
#'  @details 
#'  The function will not stop when an incorrect pathway identifier is used but instead will return NULL. If there is an error other than 404 error in the API then the function will fail and stop.
#'  
#'  @return Human entity count for given pathway
getHumanEntitiesCounts <- function(pathway) {
  # Error Handling
  checkClass(parameter = "pathway", value = pathway, classType = "character", functionName = "getHumanEntitiesCounts")
  
  human_ent_query <- paste0("https://reactome.org/ContentService/data/participants/", pathway, "/referenceEntities")
  human_ent_result <- callApi(query = human_ent_query, failOn404 = FALSE)
  ent_total <- 0
  # filter count to only include participants from the databases UniProt, ChEBI, EMBL, ENSEMBL and miRBase.
  for(entry in human_ent_result) {
    if(entry$databaseName %in% c("UniProt", "ChEBI", "EMBL", "ENSEMBL", "miRBase")) {
      ent_total = ent_total + 1
    }
  }
  return(ent_total)
}

#' Get Human Reactions Counts
#'
#' @description Run Reactome content service APIs to extract Human pathway reaction counts for a single given pathway identifier.
#' The API will pull the data of the live version of Reactome and so this function only needs to be run when the HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
#' The function calls "https://reactome.org/ContentService/data/pathway/path_id/containedEvents" and filters events returned to only those which have 'className' = Reaction.
#' This filtering occurs to match those counts produced in the Reactome species comparison tool. Finally the number of events returned after filtering is set as the Human reaction count for the given pathway.
#' This function is intended to be run for pathways which no count data is found across other species when updating the entities and reaction file used in the tool.
#' @param pathway Reactome pathway identifier of interested
#' @importFrom logger log_warn 
#' @include global_functions.R
#' @details The function will not stop when an incorrect pathway identifier is used but instead will return NULL. If there is an error other than 404 error in the API then the function will fail and stop.
#' 
#' @return Human reaction count for given pathway
getHumanReactionsCounts <- function(pathway) {
  # Error Handling
  checkClass(parameter = "pathway", value = pathway, classType = "character", functionName = "getHumanReactionsCounts")
  
  human_react_query <- paste0("https://reactome.org/ContentService/data/pathway/", pathway, "/containedEvents")
  human_react_result <- callApi(query = human_react_query, failOn404 = FALSE)
  
  react_total <- 0
  # filter count to only include events with class 'Reaction'.
  for(entry in human_react_result) {
    if(entry$className == "Reaction") {
      react_total = react_total + 1
    }
  }
  return(react_total)
  
}

#' Get A Single Species vs Human Entity and Reaction Counts
#'
#' @description Run Reactome analysis service API to extract Human and a given species' pathway entity and reaction counts.
#' The API will pull the data of the live version of Reactome and so this function only needs to be run when HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
#' The function calls "https://reactome.org/AnalysisService/species/homoSapiens/speciesCode?sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1" and takes counts of Human and species entities and reactions from the results.
#' 
#' @param speciesName Species name to be labeled on data.frame results column
#' @param speciesCode Species Reactome identifier code to be used in the API call
#' @include global_functions.R
#' 
#' @return Data.frame of pathway identifier, Human and species entities and reaction counts.
getHumanSpeciesEntitiesReactionsCounts <- function(speciesName, speciesCode) {
  
  # Error Handling
  checkClass(parameter = "speciesName", value = speciesName, classType = "character", functionName = "getHumanSpeciesEntitiesReactionsCounts")
  checkClass(parameter = "speciesCode", value = speciesCode, classType = "numeric", functionName = "getHumanSpeciesEntitiesReactionsCounts")
  
  columns <- c("pathway", "Human_entities" , paste0(speciesName, "_entities"), "Human_reactions", paste0(speciesName, "_reactions"))
  species_ent_react_data <- data.frame(matrix(data = NA, ncol = 5, nrow = 0))
  colnames(species_ent_react_data) <- columns
  
  species_query <- paste0("https://reactome.org/AnalysisService/species/homoSapiens/", speciesCode, "?sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1")
  species_result <- callApi(query = species_query)
  result_path <- species_result$pathways
  
  for(entry in result_path){
    entry_row <- data.frame(
      entry$stId,
      entry$entities$total,
      entry$entities$found,
      entry$reactions$total,
      entry$reactions$found)
    colnames(entry_row) <- columns
    species_ent_react_data <- rbind(species_ent_react_data, entry_row)
  }
  return(species_ent_react_data)
}


#' Get Human only Entities and Reactions counts and produce data.frame of results
#'
#' @description Run Reactome analysis service API to extract entity and reaction counts for given Human pathways.
#' The API will pull the data of the live version of Reactome and so this function only needs to be run when HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
#' This function is intended to be run to retrieve Human pathway entities and reaction counts for pathways which do not have data in any for the species comparison run previously. Pathways are passed to the function and for each pathway functions getHumanEntitiesCounts and getHumanReactionsCounts are run and results collated. 
#' @param pathways List of Reactome pathway identifier of interested
#' @include global_functions.R
#' @return Human Entity Reaction data.frame matrix for given pathways only
getHumanEntitiesReactionsCounts <- function(pathways) {
  checkClass(parameter = "pathways", value = pathways, classType = "character", functionName = "getHumanEntitiesReactionsCounts")
  ent_react_human <- data.frame(matrix(nrow = length(pathways), ncol = 3))
  colnames(ent_react_human) <- c("pathway", "Human_entities", "Human_reactions")
  ent_react_human$pathway <- pathways
  for(pathway in pathways) {
    ent_react_human[ent_react_human$pathway == pathway, "Human_entities"] <- getHumanEntitiesCounts(pathway)
    ent_react_human[ent_react_human$pathway == pathway, "Human_reactions"] <- getHumanReactionsCounts(pathway)
  }
  return(ent_react_human)
}

#' Creates and writes csv of the all species of analysis' entities and reactions counts for all Human pathways in Reactome pulling data from the Reactome API.
#'
#' @description Run Reactome APIs to extract entity and reaction counts per given pathway for Humans
#' The API will pull the data of the live version of Reactome and so this function only needs to be run when HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
#' This function will run Reactomes Analysis Service APIs "https://reactome.org/AnalysisService/species/homoSapiens/speciesCode?sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1" per species to get entities and reactions counts.
#' Where no counts exist for any of the species for a given pathway, API's "https://reactome.org/ContentService/data/participants/path_id/referenceEntities" and "https://reactome.org/ContentService/data/pathway/path_id/containedEvents" to get the Human entities and reactions counts for the missing pathways.
#' The function will write a matrix of entity and reaction counts for all pathways and species in csv to the dataLocation/Reactome_entities_reactions/ directory named ("Reactome_entities_reactions_v", version, ".csv").
#' @param version Version number to name output file with
#' @inheritParams filterInterProResults
#' @include global_functions.R
#' @return File path to written file
#' @export
createEntitiesReactionsCsv <- function(version, dataLocation = ".") {
  checkClass(parameter = "version", value = version, classType = "integer", functionName = "createEntitiesReactionsCsv")
  species_name_code <- R.utils::getOption(localVars, "species_list")
  rel_path_from_root <- R.utils::getOption(localVars, "rel_path_from_root")
  reactome_entreact_csv_path <- file.path(dataLocation, "g2pScanData", "Reactome_entities_reactions")
  dir.create(reactome_entreact_csv_path, recursive = T)
  species <- names(species_name_code)
  all_pathways_in_reactome <- callApi(query = "https://reactome.org/ContentService/data/schema/Pathway/min?species=9606&page=1&offset=20000", simplifyVector = TRUE)
  
  ent_react_data <- data.frame(matrix(nrow = nrow(all_pathways_in_reactome), ncol = 3 + length(species) * 2))
  colnames(ent_react_data) <- c("pathway", "Human_entities" , paste0(species, "_entities"), "Human_reactions", paste0(species, "_reactions"))
  ent_react_data$pathway <- sort(as.character(all_pathways_in_reactome$stId))
  for(sp in species){
    logger::log_info("Getting {sp} data")
    sp_code <- species_name_code[[sp]]$reactomeID
    species_data <- getHumanSpeciesEntitiesReactionsCounts(sp, sp_code)
    for(column in c("Human_entities", "Human_reactions", paste0(sp, "_entities"), paste0(sp, "_reactions"))) {
      values <- species_data[order(as.character(species_data$pathway)), column]
      ent_react_data[ent_react_data$pathway %in% species_data$pathway , column] <- values
    }
  }
  
  # get missing human entity counts where species show no data
  no_ent_data_paths <- ent_react_data[is.na(ent_react_data$Human_entities), "pathway"]
  logger::log_info("Getting Human Entities and Reactions")
  missing_human_data <- getHumanEntitiesReactionsCounts(no_ent_data_paths)
  # add missing counts to final matrix
  for(column in c("Human_entities", "Human_reactions")) {
    values <- missing_human_data[order(as.character(missing_human_data$pathway)), column]
    ent_react_data[ent_react_data$pathway %in% missing_human_data$pathway , column] <- values
  }
  ent_react_data[is.na(ent_react_data)] <- 0 # Where no data is available count is 0.
  
  # Write matrix to file
  file_name <- file.path(reactome_entreact_csv_path, paste0("Reactome_entities_reactions_v", version, ".csv"))
  write.csv(x = ent_react_data, file = file_name, row.names = F)
  logger::log_info("Entities and Reactions file written to {file_name}.  This is of Reactome's most up to date entities and reactions information using version {version}.")
  return(file_name)
}
