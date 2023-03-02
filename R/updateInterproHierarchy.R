# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get InterPro family branches and hierarchy levels of InterPro families
#'
#' @description Build data.frame representation of the InterPro Hierarchy using 'ParentChildTreeText.txt' (from InterPro ftp site) data.frame passed to function. Each family will be given as a row in the data.frame with InterPro ID, name, hierarchy level and family branch/group.
#' The order of the hierarchy is maintained in the data.frame. A family branch number represents one branch of the family tree and every family with the same family branch number belong to the same most parental family. The level indicates which generation the family is of within the family tree of that branch.
#' Level 1 means the family is the most parental family. Level 2 means the family is a child of the level 1 family (1st generation child). Level 3 means that the family is the child of a level 2 family (2nd generation child). Levels are defined as number of parents of an individual family plus 1. The number of parents of a family is determined by the number of "--" prefixing the family identifier.
#'
#' @param dt A data.frame of the InterPro's 'ParentChildTreeText.txt' (from InterPro ftp site). Quotes should be removed from data in this data.frame. InterPro IDs should be given in column V1 and order of these ids along with "--" indicating hierarchical level must remain.
#' @include global_functions.R
#' @return data.frame of family identifier, family name, hierarchy level and hierarchy branch number
getInterproLevels <- function(dt) {
  checkClass(parameter = "dt", value = dt, classType = "data.frame", functionName = "getInterproLevels")
  family_levels <- data.frame(data=NA,stringsAsFactors = FALSE ) 
  familyCount <- 0
  for(i in 1:nrow(dt))
  {
    fam <- dt[i,1] # get identifier
    family_levels[i,1] <- as.character(fam)
    family_levels[i,2] <- dt[i,2]
    level <- grepl("--",fam) # Has the identifier got '--'
    if(!level){ 
      family_levels[i,3] <-1 # When there is no '--' the level is 1
      familyCount <- familyCount+1
      family_levels[i,4] <- familyCount
    }
    else{
      family_levels[i,3] <- length(gregexpr("--",fam)[[1]])+1 # When '--' are found the level is the count +1
      family_levels[i,4]<- familyCount
    }
  }
  family_levels$data <- gsub("--", "", family_levels$data)
  colnames(family_levels) <- c("interpro_id", "family_names", "hierarchy_level", "family_branch")
  return(family_levels)
}

#' Update InterPro Hierarchy csv
#' 
#'
#' @description Looks for a InterPro hierarchy csv of the current InterPro API database version at directory dataLocation/interpro_hierarchy_with_levels.
#' If correct version is found then this csv version is set as the InterPro hierarchy global variable for use during analysis, if it is not found then a new hierarchy csv is created by downloading hierarchy from "https://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt" and transforming data to included hierarchy information columns.
#' @details The InterPro API database version is determined by the API query "https://www.ebi.ac.uk/interpro/api/protein/uniprot/" and this version is searched for in the directory dataLocation/interpro_hierarchy_with_levels. Versions of the csv files are identified by the characters after the last underscore of the file names. If version number of the database match any of the file names suffixes then that file is identified as the most up to data InterPro hierarchy csv.
#' When there is no match between InterPro database require version number and files found, a new hierarchy file will be written to the directory dataLocation/interpro_hierarchy_with_levels with file name 'Interpro_hierarchy_\emph{versionNumber}.csv'.
#' A new hierarchy file is created from downloaded file from the InterPro  ftp site at "https://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt" and then each branch of the hierarchy in this file is processed to identify it's level. Levels are defined as number of parents of an individual family plus 1. The number of parents of a family in the downloaded file is determined by the number of "--" prefixing the family identifier.
#' A new table of family identifiers, the branch number it belongs to and the hierarchy level of the family is created and it is this file saved as the new InterPro hierarchy file for use in the pipeline. 
#' The hierarchy file determined as the current InterPro version has it's file path saved to R.utils' option 'interpro_levels' for later use.
#' @inheritParams filterInterProResults
#' @importFrom R.utils downloadFile setOption getOption
#' @importFrom logger log_info log_fatal
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @export
updateInterproHierarchy <- function(dataLocation = ".") {
  interpro_version_query <- query <- "https://www.ebi.ac.uk/interpro/api/protein/uniprot/"
  interpro_version_result <- query_result <- httr::RETRY("GET", query, encode = "application/json" )
  interpro_version <- query_result[["all_headers"]][[1]][["headers"]][["interpro-version"]]
  
  family_rData_path <- file.path(dataLocation, "g2pScanData", "interpro_hierarchy_with_levels")
  data_available <- list.files(family_rData_path)
  update = FALSE
  if(length(data_available) != 0) {
    data_available_version <- stringr::str_replace_all(data_available, c("^.*_" = "", ".csv" = ""))
    if(any(grepl(interpro_version, data_available_version))) {
      logger::log_info("InterPro hierarchy version {interpro_version} already exists - keeping version.")
      interpro_hierarchy_path <- file.path(family_rData_path,
                                           data_available[grepl(interpro_version, data_available_version)])
      
    } else {
      update = TRUE
    }
  } else {
    update = TRUE
  }
  
  if(update) {
    logger::log_info("InterPro hierarchy being updated with version {interpro_version}.")
    return_val <- R.utils::downloadFile("https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt",
                                        method = "auto")
    if(return_val == "ParentChildTreeFile.txt") {
      tree_text <- readLines(return_val) # '::' is the delimiter between columns but ':' is found in family names. So to make sep one byte change "::" to ";"
      tree_text <- gsub("::",";",tree_text) # Make ; the delimiter as : is found in family names
      tree_text <- gsub('"',"", tree_text) # Remove string quotes
      file.remove(file.path(return_val)) # delete downloaded file.
      family_tree_formatted_text <-  as.data.frame(do.call( rbind, strsplit(tree_text,';' )), stringAsFactors=F) 
      family_levels <- getInterproLevels(family_tree_formatted_text)
      interpro_hierarchy_path <- file.path(family_rData_path, paste0("Interpro_hierarchy_", interpro_version, ".csv"))
      if(!dir.exists(family_rData_path)) {
        dir.create(family_rData_path, recursive = T)
      }
      write.csv(family_levels, file = interpro_hierarchy_path, row.names = F)
      logger::log_info("InterPro hierarchy update complete. Current version: {interpro_version}.")
      
    } else {
      logger::log_fatal("Unable to get InterPro ParentChildTreeFile.txt from 'https://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt'. Update of hierarchy can not proceed.")
      stop("Unable to get InterPro ParentChildTreeFile.txt from 'https://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt'. Update of hierarchy can not proceed.")
    }
  }
  
  R.utils::setOption(localVars, "interpro_levels", 
                     interpro_levels <- interpro_hierarchy_path)
  return(interpro_hierarchy_path)
}



