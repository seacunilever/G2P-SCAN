# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get version number of HumanMine
#' @importFrom InterMineR getVersion
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @return Version number
getHumanMineVersion <- function() {
  hm <- R.utils::getOption(localVars, "hm")
  version <- InterMineR::getVersion(hm)
  return(version)
}

#' Get version number of the Reactome database currently used in HumanMine
#'
#' @include global_functions.R
#' @return Version number
getHumanMineReactomeVersion <- function() {
  #get current Human Reactome version
  reactome_version_result <- runHumanMineQuery(templateName = "Gene_Pathway",
                                               constraintPath = c("Gene.pathways.dataSets.dataSource.dataSets.name", "Gene.symbol"),
                                               constraintOperators = c("=", "="),
                                               constraintValues = list("Reactome", "BCHE"),
                                               select = c("Gene.symbol", "Gene.pathways.dataSets.dataSource.dataSets.version")
  )
  humanMine_reactome_version <- unique(reactome_version_result[, "Gene.pathways.dataSets.dataSource.dataSets.version"])
  humanMine_reactome_version_number <- as.numeric(stringr::str_replace_all(humanMine_reactome_version, c(";.*" = "")))
  return(humanMine_reactome_version_number)
}

#' Get version number of Reactome resources used in analysis
#' @param var "reactome_datatree" or "reactome_ent_react_file" to pull version of
#' @importFrom R.utils getOption 
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @return Version number
getVersionOfReactomeDataUsed <- function(var) {
  if(!var %in% c("reactome_datatree", "reactome_ent_react_file")) {
    logger::log_fatal("Variable var must be either 'reactome_datatree' or 'reactome_ent_react_file'.")
  }
  tryCatch({
  used_version<- R.utils::getOption(localVars, var)
  version <- as.numeric(stringr::str_replace_all(used_version, c(".*_v" = "", "[.].*$" = "")))
  
  return(version)
  }, error = function(e){
    if(var == "reactome_datatree") {
      function_name <- 'updateReactomeHierarchy()'
    } else {
      function_name <- 'updateEntitiesReactionsCounts()'
    }
    logger::log_info("No {var} variable is set in localVars yet, please run {function_name}.")
  })
}


#' Get version number of Panther currently used in HumanMine
#'
#' @include global_functions.R
#' @return Version number
getHumanMineOrthologueVersion <- function() {
  #get current Human Reactome version
  orthologue_version_result <- runHumanMineQuery(templateName = "Gene_Orth",
                                                      constraintPath = c("Gene.symbol", "Gene.homologues.dataSets.dataSource.dataSets.name"),
                                                      constraintOperators = c("=", "="),
                                                      constraintValues = list("BCHE", "panther"),
                                                      constraintm.index = 1,
                                                      select = c("Gene.primaryIdentifier","Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.dataSets.dataSource.dataSets.version")
  )
  orthologue_version <- orthologue_version_result$Gene.homologues.dataSets.dataSource.dataSets.version[1]
  orthologue_version_number <- as.numeric(stringr::str_replace_all(orthologue_version, c(";.*" = "")))
  return(orthologue_version)
}

#' Get version number of UniProt used by the UniProt API
#' @importFrom httr RETRY
#' @include global_functions.R
#' @return UniProt release date
getUniProtVersion <- function() {
  queryResult <- httr::HEAD("https://rest.uniprot.org/uniprotkb/P12345.rdf?include=yes") # random query
  uniprotVersion <- queryResult$headers$`x-uniprot-release`
  return(uniprotVersion)
}

#' Get version number of InterPro used by the InterPro API and sub databses used by InterPro
#' @include global_functions.R
#' @return List of all InterPro versions, including InterPro itself, the protein instances and member databases with name and version and description.
getInterProVersions <- function() {
  query <-'https://www.ebi.ac.uk/interpro/api/'
  query_result <- callApi(query)
  databases <- query_result$databases
  version_list <- list()
  dont_include_dbs <- c("CATH-Gene3D", "CDD", "PROSITE profiles", "PROSITE patterns", "SMART", "SUPERFAMILY")
  for(db in databases) {
    if(db$type == "entry" | db$type == "protein") {
      if(!db$name %in% dont_include_dbs) {
      version_list[[db$name]]["version"] <- db$version
      version_list[[db$name]]["description"] <- db$description
      }
    }
  }
  return(version_list)
}


#' Write openxlsx worksheet of analysis database versions
#' @description Write openxlsx worksheet sheet giving all versions of databases used in the G2P-SCAN analysis
#' @param wb A openxlsx workbook object
#' @importFrom openxlsx addWorksheet writeDataTable
#' @return wb with addition worksheet appended
#' @export
writeVersionWb <- function(wb) {
  hm <-getHumanMineVersion()
  hm_reactome <- getHumanMineReactomeVersion()
  hm_reactome_hierachy_used <- getVersionOfReactomeDataUsed("reactome_datatree")
  hm_reactome_ent_react_used <- getVersionOfReactomeDataUsed("reactome_ent_react_file")
  hm_ortho <- getHumanMineOrthologueVersion()
  uniprot <- getUniProtVersion()
  
  interpro <- t(as.data.frame(getInterProVersions()))
  interpro <- interpro[c("InterPro", rownames(interpro)[rownames(interpro) != "InterPro"]), ]
  all_vals <- list(hm,
                   hm_reactome,
                   hm_reactome_hierachy_used,
                   hm_reactome_ent_react_used,
                   hm_ortho,
                   uniprot,
                   interpro[, "version"]) 
  any_null <- unlist(lapply(all_vals, is.null))
  if(any(any_null)) {
    all_vals[any_null] <- "Unknown"
  }
  version_dt <- data.frame("Type" = c("InterMine",
                        "InterMine",
                        "Static stored data",
                        "Static stored data",
                        "InterMine",
                        "Uniprot API",
                        "InterPro API",
                        rep("InterPro member databases", nrow(interpro)-1)
                        ), 
             "Data source" = c("HumanMine",
                               "HumanMine installment of Reactome",
                               "Reactome Hierarchy",
                               "Reactome Entities & Reactions",
                               "HumanMine installment of Panther",
                               "UniProt",
                               rownames(interpro)),
             "Version" = unlist(all_vals))
  
  version_sheet <- openxlsx::addWorksheet(wb = wb, sheetName = "DB versions")  
  openxlsx::writeDataTable(wb = wb, version_sheet, x = version_dt, rowNames = F, na.string = "")
  return(wb)
}

