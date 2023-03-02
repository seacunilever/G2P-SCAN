# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Filter Entities' and Reaction's counts for a list of pathways
#'
#' @description Reads in matrix of all pathway and species (created by updateEntitiesReactionCounts function) and filter for the pathways of interest.
#' @param pathways List of Reactome pathway identifiers to filter by. When NULL no filter is applied and all pathways are returned
#' @param dataToFilter Data.frame to filter. Rownames should be pathways
#' 
#' @importFrom R.utils getOption
#' @importFrom logger log_fatal log_warn
#' @include global_functions.R
#' @return Filtered entity and reaction matrix
filterEntitiesReactionsPathways <- function(dataToFilter, pathways = NULL) {
  reactome_v <- R.utils::getOption(localVars, "reactome_version")
  checkClass(parameter = "dataToFilter", value = dataToFilter, classType = "data.frame", functionName = "filterEntitiesReactionsPathways")
  if(!is.null(pathways)) {
    checkClass(parameter = "pathways", value = pathways, classType = "character", functionName = "filterEntitiesReactionsPathways")
  } else {
    pathways <- rownames(dataToFilter)
  }
  # Filter Pathways
  filtered_data <- dataToFilter[rownames(dataToFilter) %in% pathways, ]
  not_found <- pathways[!pathways %in% rownames(dataToFilter)]
  if(length(not_found) > 0 ) {
    if(nrow(filtered_data) == 0) {
      stop(paste0("No variables given as pathways are identified as a pathway in Reactome version ", reactome_v, "."))
      logger::log_fatal("No variables given as pathways are identified as a pathway in Reactome version {reactome_v}.")
    } else {
      logger::log_warn("Some variables given as pathways are not identified as a pathway in Reactome version {reactome_v}. These pathways will have no entity and reaction counts because of this. Pathways not found: {not_found}")
      for(p in not_found) {
        filtered_data[p,] <- rep(NA, ncol(filtered_data))
      }
    }
  }
  return(filtered_data)
}

#' Filter Entities and Reaction counts for a list of species
#'
#' @description Reads in matrix of all pathway and species (created by updateEntitiesReactionCounts function) and filters for the species of interest.
#' @inheritParams filterOrthologueSpecies
#' @param dataToFilter data.frame to filter. Column names must include species with either '_entities' or '_reactions' suffixes.
#' 
#' @importFrom R.utils getOption
#' @importFrom logger log_fatal
#' @include global_functions.R
#' @return Filtered entity and reaction matrix
filterEntitiesReactionsSpecies <- function(dataToFilter, species = NULL) {
  reactome_v <- R.utils::getOption(localVars, "reactome_version")
  checkClass(parameter = "dataToFilter", value = dataToFilter, classType = "data.frame", functionName = "filterEntitiesReactionsPathways")
  if(!any(c("Human_entities", "Human_reactions") %in% colnames(dataToFilter))) {
    stop("Human data not found in entity and reaction matrix.")
    logger:log_fatal("Human data not found in entity and reaction matrix. Please check matrix.")
    
  }
  species <- checkSpecies(species = species, functionName = "filterEntitiesReactionsSpecies")
  
  #Filter Species
  species_filter_list <- c("Human_entities", "Human_reactions", paste0(species, "_entities"), paste0(species, "_reactions"))
  filtered_data <- dataToFilter[, colnames(dataToFilter) %in% species_filter_list]
  
  return(filtered_data)
}


#' Get Entities and Reaction counts for a list of pathways and species
#'
#' @description Reads in matrix of all pathway and species (created by updateEntitiesReactionCounts function) and filter for the pathways and species of interest.
#' @inheritParams filterEntitiesReactionsSpecies
#' @inheritParams filterEntitiesReactionsPathways 
#' @inheritParams filterInterProResults
#' @importFrom R.utils getOption
#' @include global_functions.R
#' @return Filtered entity and reaction matrix
#' @export
getEntitiesReactionsCounts <- function(pathways = NULL, species = NULL, onMissedUpdate = "update", dataLocation = ".") {
  updateEntitiesReactionsCounts(onMissedUpdate, dataLocation = dataLocation) # set localVars reactome_ent_react_file to get path to right entities and reactions file - update if needed.
  ent_react_path <- R.utils::getOption(localVars, "reactome_ent_react_file")
  ent_react_data <- read.csv(ent_react_path, row.names = 1)
  colnames(ent_react_data) <- gsub("[.]{2}", ". ", colnames(ent_react_data))
  
  ent_react_data_path_filtered <- filterEntitiesReactionsPathways(dataToFilter = ent_react_data, pathways = pathways)
  ent_react_data_path_species_filtered <- filterEntitiesReactionsSpecies(dataToFilter = ent_react_data_path_filtered, species = species)
  return(ent_react_data_path_species_filtered)
}
