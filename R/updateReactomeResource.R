# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Get file path of specified version of Reactome pathway hierarchy tree if available
#'
#' @description This function looks in the dataLocation/g2pScanData folder and in the directory specified from the 'folder' parameter and searches for data.tree with same version as the 'version' parameter specifies. Versions numbers of the files are determined by the characters after the last underscore in the file name. If a file of matching version is available that file path is returned.
#' If matching version is not found then NULL is returned.
#' When parameter version is set to string "latest", then the files are searched and the one with the highest version number has it's file path returned.
#' @param dataName The string of what data is being searched for. Used for logger messages.
#' @param folder The folder location to search for files in the package 'data' folder. 
#' @param version Version number or 'latest' to get specific versions
#' @inheritParams filterInterProResults
#' @importFrom logger log_info
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @include global_functions.R
#' @return File path of Rdata file of Reactome pathway hierarchy if available, NULL if not.
getReactomeFileOfVersionPath <- function(dataName, folder, version, dataLocation = ".") {
  if(!(version %in% "latest" & is.character (version))) {
    checkClass(parameter = "version", value = version, classType = "numeric", functionName = "getReatomeHierarchyOfVersionPath")
  }
  reactome_hierarchy_rData_path <- file.path(dataLocation, "g2pScanData", folder)
  if(!dir.exists(reactome_hierarchy_rData_path)) {
    dir.create(path = reactome_hierarchy_rData_path, recursive = T)
  }
  data_available <- list.files(reactome_hierarchy_rData_path)
  update = FALSE
  if(length(data_available) != 0) {
    data_available_version <- as.numeric(stringr::str_replace_all(data_available, c(".*_v" = "", "[.].*$" = "")))
    if(version != "latest") {
      if(any(grepl(version, data_available_version))) {
        logger::log_info("{dataName} version {version} already exists - using this version.")
        reactome_hierarchy_path <- file.path(reactome_hierarchy_rData_path,
                                             data_available[grepl(version, data_available_version)])

      } else {
        reactome_hierarchy_path = NULL
      }
      
    } else {
      version <- max(data_available_version)
      reactome_hierarchy_path <- file.path(reactome_hierarchy_rData_path,
                                           data_available[grepl(version, data_available_version)])
      logger::log_info("{dataName} is using version {version} of Reactome. This is the latest hierarchy saved in the package.")
      
    }
    
  } else {
    reactome_hierarchy_path = NULL
  }
  if(is.null(reactome_hierarchy_path)) {
    logger::log_info("{dataName} version {version} is not found.")
  }
  return(reactome_hierarchy_path)
}

#' Update and write Reactome Pathway Hierarchy data.tree
#' 
#'
#' @description Create and save Reactome hierarchy data.tree using latest Reactome version data from 'https://reactome.org/download/current/ReactomePathwaysRelation.txt'. File is download and read by function and transformed in corresponding data.tree hierarchy.
#' The data.tree has first node 'topNode' which then the most parental Reactome pathways branch from. The data.tree is written as Rdata object file to 'dataLocation/g2pScanData/pathway_hierarchy_datatrees' location. The full file path of this file is returned by the function.
#' @param version Version number of Reactome used for file name output
#' @inheritParams filterInterProResults
#' @importFrom data.tree Node 
#' @importFrom R.utils downloadFile setOption getOption
#' @importFrom logger log_info log_fatal
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @include global_functions.R
#' @return Full file path of newly written data.tree file
createReactomeHierarchyDataTree <- function(version, dataLocation = ".") {
  checkClass(parameter = "version", value = version, classType = "integer", functionName = "createReactomeHierarchyDataTree")
  reactome_hierarchy_rData_path <- file.path(dataLocation, "g2pScanData", "pathway_hierarchy_datatrees")
  logger::log_info("Downloading pathway hierarchy txt file from Reactome...")
  return_val <- R.utils::downloadFile("https://reactome.org/download/current/ReactomePathwaysRelation.txt",
                                      method = "auto")
  if(return_val == "ReactomePathwaysRelation.txt") {
    tree_dt<- read.delim(return_val, header = F)
    file.remove(return_val) # delete downloaded file.
    pathTree_Human <- tree_dt[grepl("R-HSA*",tree_dt[,1]), c(1, 2)]
    parent_not_child <- as.character(unique(pathTree_Human[!(pathTree_Human[, 1] %in% pathTree_Human[, 2]), 1]))
    colnames(pathTree_Human) <- c("parent", "child")
    
    #create data.tree
    topNode <- data.tree::Node$new("top")
    for(p in parent_not_child){
      
      assign(p, topNode$AddChild(p))
      subChildTable <- pathTree_Human[pathTree_Human[, 1] == p, 2]
      for(child in subChildTable){
        assign(child, get(p)$AddChild(child))
      }
    }
    
    pathTree_Human <- pathTree_Human[!(pathTree_Human[, 1] %in% parent_not_child),]
    parent_not_child <- as.character(unique(pathTree_Human[!(pathTree_Human[,1] %in% pathTree_Human[, 2]), 1]))
    while(nrow(pathTree_Human) > 0){
      for(p in parent_not_child){
        subChildTable <- pathTree_Human[pathTree_Human[, 1] == p, 2]
        for(child in subChildTable){
          assign(child, get(p)$AddChild(child))
        }
      }
      pathTree_Human <- pathTree_Human[!(pathTree_Human[, 1] %in% parent_not_child),]
      parent_not_child <- as.character(unique(pathTree_Human[!(pathTree_Human[, 1] %in% pathTree_Human[, 2]), 1]))
    }
    
    #write data.tree object
    file_name <- file.path(reactome_hierarchy_rData_path, paste0("pathway_hierarchy_reactome_v", version, ".Rdata"))
    save(topNode, file = file_name)
    logger::log_info("Rdata file for reactome hierarcy data.tree is saved to '{file_name}'. This is of Reactome's most up to date pathway information using version {version}.")
  } else {
    logger::log_fatal("Unable to get Reactome ReactomePathwaysRelation.txt from 'https://reactome.org/download/current/ReactomePathwaysRelation.txt'. Update of hierarchy can not proceed.")
    stop("Unable to get InterPro ReactomePathwaysRelation.txt from 'https://reactome.org/download/current/ReactomePathwaysRelation.txt'. Update of hierarchy can not proceed.")
  }
  return(file_name)
}

#' Update Reactome Pathway Hierarchy data.tree for package use
#' 
#'
#' @description Update and save Reactome Hierarchy data.tree if the version of Reactome in HumanMine has changed since last data.tree was saved. This function looks in the package's dataLocation/g2pScanData/pathway_hierarchy_datatrees location and searches for data.tree with same version of current Reactome in HumanMine (obtained by InterMineR call). If a data.tree of matching version is available no update is needed and global R.utils option variable reactome_datatree is set as the matching data.tree full file path.Versions numbers of the files are determined by the characters after the last underscore in the file name.
#' If matching version file is not found, hierarchy text file is downloaded from "https://reactome.org/download/current/ReactomePathwaysRelation.txt" and transformed into a data.tree, matching the current version of the current HumanMine Reactome instance. If Reactome version and HumanMine Reactome instance versions do not match, because an update has been missed, there's an option to use the most up to date data.tree found in the package's data/pathway_hierarchy_datatrees directory or the create a new data.tree using the most up to data version from Reactome.
#' Any new data.trees created are saved as a Rdata file at dataLocation/g2pScanData/pathway_hierarchy_datatrees of the package directory and this file path is set as the R.utils option variable reactome_datatree of localVars for use throughout the package.
#' @param onMissedUpdate Character string to indicate what to do when Reactome file cannot be updated to appropriate version.
#' Either "update" to use the Reactome API most up to date version. "use latest" to use the most up to date version of the file saved in the dataLocation appropriate repository already or given an integer to state what version of Reactome file to use - this must be already present in the dataLocation appropriate repository. 
#' @inheritParams filterInterProResults
#' @importFrom data.tree Node 
#' @importFrom R.utils setOption
#' @importFrom logger log_warn log_fatal
#' @importFrom stringr str_replace_all
#' @include global_variables.R
#' @include global_functions.R
#' @include getVersions.R
#' @export
updateReactomeHierarchy <- function(onMissedUpdate = "update", dataLocation = ".") {
  run_update <- FALSE
  if(!onMissedUpdate %in% c("update", "use latest")) {
    if(class(onMissedUpdate) != "numeric") {
      e <- "onMissedUpdate parameter for function updateReactomeHierarchy includes invalid input: '{onMissedUpdate}'. Must be only 'update' or 'use latest' or a given version number."
      logger::log_fatal(e)
      stop(e)
      }
  }
  logger::log_info("Checking to see if the Reactome pathway hierarchy needs updating...")
  
  humanMine_reactome_version_number <- getHumanMineReactomeVersion() 
  datatree_up_to_date_path <- getReactomeFileOfVersionPath(dataName = "Reactome pathway hierarchy",
                                                           folder = "pathway_hierarchy_datatrees",
                                                           version = humanMine_reactome_version_number,
                                                           dataLocation = dataLocation)
  if(is.null(datatree_up_to_date_path)) {
    #check current live Reactome version 
    reactome_version <-  callApi("https://reactome.org/AnalysisService/database/version")
    if(reactome_version == humanMine_reactome_version_number) {
      ##update
      logger::log_info("Reactome API available is {reactome_version} so correct update can proceed.")
      run_update <- TRUE
    } else {
      # cannot update correctly as current reactome api version if not of the same version we need
      if(onMissedUpdate == "update") {
        logger::log_info("Looking for Reactome entities and reactions file using Reactome API current version ({reactome_version}). ")
        
        current_aval <- getReactomeFileOfVersionPath(dataName =  "Reactome pathway hierarchy",
                                                     folder = "pathway_hierarchy_datatrees",
                                                     version = as.numeric(reactome_version),
                                                     dataLocation = dataLocation)

        if(!is.null(current_aval)) {
          hierarchy_use_path <- current_aval
        } else {
          #save new data.tree based on reactome API version.
          logger::log_info("Reactome API available is {reactome_version} and this version will be used to update the Reactome pathway hierarcy.")
          run_update <- TRUE
          logger::log_warn("The version of the Reactome pathway hierarchy is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Pathway Hierarchy being used: {reactome_version}, HumanMine version : {humanMine_reactome_version_number}.")
        } 
      }
      else if(class(onMissedUpdate) == "numeric") {
        #get file path wit highest version.
        logger::log_info("Reactome API available is {reactome_version} which does not match required version. Instead version {onMissedUpdate} has been selected to be used for the pathway hierarchy.")
        
        hierarchy_use_path <- getReactomeFileOfVersionPath(dataName = "Reactome pathway hierarchy",
                                                           folder = "pathway_hierarchy_datatrees",
                                                           version = onMissedUpdate,
                                                           dataLocation = dataLocation)
        if(!is.null(hierarchy_use_path)) {
          logger::log_info("The Reactome pathway hierarchy file of version {onMissedUpdate} is found and will be used for analysis.")
          logger::log_warn("The version of the Reactome pathway hierarchy is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Pathway Hierarchy being used: {onMissedUpdate}, HumanMine version : {humanMine_reactome_version_number}.")
        } else{
          logger::log_info("The Reactome pathway hierarchy file of version {onMissedUpdate} is NOT found. Reactome API version file will be looked for instead.")
          api_version <- getReactomeFileOfVersionPath(dataName =  "Reactome pathway hierarchy",
                                                       folder = "pathway_hierarchy_datatrees",
                                                       version = as.numeric(reactome_version),
                                                       dataLocation = dataLocation)
          if(is.null(api_version)){
            run_update <- TRUE
          } else {
            logger::log_warn("The version of the Reactome pathway hierarchy is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Pathway Hierarchy being used: {reactome_version}, HumanMine version : {humanMine_reactome_version_number}.")
            
            hierarchy_use_path <- api_version
          }
        }
      } else {
        #get file path wit highest version.
        logger::log_info("Reactome API available is {reactome_version} which does not match required version. Using the most up to data already available pathway hierarchy will be used instead.")
        
        hierarchy_use_path <- getReactomeFileOfVersionPath(dataName =  "Reactome pathway hierarchy",
                                                           folder = "pathway_hierarchy_datatrees",
                                                           version = "latest",
                                                           dataLocation = dataLocation)
        version <- as.numeric(stringr::str_replace_all(hierarchy_use_path, c("^.*_v" = "", "[.].*$" = "")))
        logger::log_warn("The version of the Reactome pathway hierarchy is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Pathway Hierarchy being used: {version}, HumanMine version : {humanMine_reactome_version_number}.")
        
      }
    }
  } else {
    hierarchy_use_path <- datatree_up_to_date_path
  }
  if(run_update) {
    logger::log_info("Updating pathway hierarchy file using Reactome API with data version {reactome_version}.")
    hierarchy_use_path <- createReactomeHierarchyDataTree(reactome_version, dataLocation = dataLocation)
  } 
  R.utils::setOption(localVars, "reactome_datatree", 
                     "reactome_datatree" <- hierarchy_use_path)
  return(hierarchy_use_path)
}

#' Get Entities and Reaction count matrix for all Species vs Human
#'
#' @description Update and save Reactome entity and reaction counts if the version of Reactome in HumanMine has changed since last data.tree was saved. This function looks in the package's data/Reactome_entities_reactions location and searches for csv with same version of current Reactome in HumanMine (obtained by InterMineR call). If a csv of matching version is available no update is needed and global R.utils option variable reactome_ent_react_file is set as the matching csv full file path.
#' If matching versioned file is not found, Reactome API is used to get entities and reactions and transformed into the appropriate format. If Reactome version and HumanMine Reactome instance versions do not match, because an update has been missed, there's an option to use the most up to date csv found in the package's data/Reactome_entities_reactions directory or the create a file using the most up to data version from Reactome.
#' Any new csv created are saved at data/Reactome_entities_reactions of the package directory and this file path is set as the R.utils option variable reactome_ent_react_file of localVars for use throughout the package.
#' Update of the entities and reaction is performed by running Reactome APIs to extract entity and Reaction counts per pathway for Human and all 6 species to be used in the Genes to pathways tool -
#' C. elegans, D. rerio, D. melanogaster, M. musculus, R. norvegicus and S. cerevisiae. 
#' 
#' @param outputDir file path to directory where file should be saved
#' @inheritParams filterInterProResults
#' @importFrom R.utils getOption
#' @importFrom logger log_info
#' @include global_variables.R
#' @include getReactomeEntitiesReactions.R
#' @return Entity Reaction matrix
#' @export
updateEntitiesReactionsCounts <- function(onMissedUpdate = "update", dataLocation = ".") {
  run_update <- FALSE
  if(!onMissedUpdate %in% c("update", "use latest")) {
    if(class(onMissedUpdate) != "numeric") {
      e <- "onMissedUpdate parameter for function updateEntitiesReactionsCounts includes invalid input: '{onMissedUpdate}'. Must be only 'update' or 'use latest' or a given version number."
      logger::log_fatal(e)
      stop(e)
      }
  }
  logger::log_info("Checking to see if the Reactome entities and reactions file needs updating...")
  humanMine_reactome_version_number <- getHumanMineReactomeVersion() 
  csv_up_to_date_path <- getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file",
                                                      folder = "Reactome_entities_reactions",
                                                      version = humanMine_reactome_version_number,
                                                      dataLocation = dataLocation)
  
  
  reactome_version <-  callApi("https://reactome.org/AnalysisService/database/version")
  
  if(is.null(csv_up_to_date_path)) {
    #check current live Reactome version 
    reactome_version <-  callApi("https://reactome.org/AnalysisService/database/version")
    if(reactome_version == humanMine_reactome_version_number) {
      ##update
      logger::log_info("Reactome API available is {reactome_version} so correct update can proceed.")
      run_update <- TRUE
    } else {
      logger::log_info("Reactome API available is {reactome_version} and does not match required HumanMine reactome version {humanMine_reactome_version_number}.")
      
      # cannot update correctly as current reactome api version if not of the same version we need
      if(onMissedUpdate == "update") {
        logger::log_info("Looking for Reactome entities and reactions file using Reactome API current version ({reactome_version}). ")
        
        current_aval <- getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file",
                                                     folder = "Reactome_entities_reactions",
                                                     version = as.numeric(reactome_version),
                                                     dataLocation = dataLocation)
        if(!is.null(current_aval)) {
          enreact_use_path <- current_aval
        } else {
          #save new data.tree based on reactome API version.
          logger::log_info("Reactome API available is {reactome_version} and this version will be used to update the Reactome entities and reactions file.")
          run_update <- TRUE 
        }
      }else if(class(onMissedUpdate) == "numeric") {
        #get file path wit highest version.
        logger::log_info("Reactome API available is {reactome_version} which does not match required version. Instead version {onMissedUpdate} has been selected to be used for the entities and reactions file.")
        logger::log_info("Looking to see entities and reactions file of version {onMissedUpdate} already exists.")
        
        enreact_use_path <- getReactomeFileOfVersionPath(dataName =  "Reactome entities and reactions file",
                                                         folder = "Reactome_entities_reactions",
                                                         version = onMissedUpdate,
                                                         dataLocation = dataLocation)
        if(!is.null(enreact_use_path)) {
          logger::log_info("Entities and reactions file of version {onMissedUpdate} is found and will be used for analysis.")
          logger::log_warn("The version of the Reactome entities and reactions file is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Entities and reactions file being used: {onMissedUpdate}, HumanMine version : {humanMine_reactome_version_number}.")
        } else{
          logger::log_info("Entities and reactions file of version {onMissedUpdate} is NOT found. Reactome API version file will be looked for instead.")
          api_version <- getReactomeFileOfVersionPath(dataName =  "Reactome entities and reactions file",
                                                      folder = "Reactome_entities_reactions",
                                                      version = as.numeric(reactome_version),
                                                      dataLocation = dataLocation)
          if(is.null(api_version)){
            run_update <- TRUE
          } else {
            logger::log_warn("The version of the Reactome entities and reactions file is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Entities and reactions file being used: {reactome_version}, HumanMine version : {humanMine_reactome_version_number}.")
            enreact_use_path <- api_version
          }
      
        }
      }else {
        #get file path wit highest version.
        logger::log_info("Reactome API available is {reactome_version} which does no match required version. Most up to data already available file of entities and reactions will be used instead.")
        enreact_use_path <- getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file",
                                                         folder = "Reactome_entities_reactions",
                                                         version = "latest",
                                                         dataLocation = dataLocation)
        version <- as.numeric(stringr::str_replace_all(enreact_use_path, c("^.*_v" = "", "[.].*$" = "")))
        logger::log_warn("The version of the Reactome pathway entities and reactions data is different to the version used by HumanMine and therefore some mismatching of pathway identifiers is possible. Entities and reactions version being used: {version}, HumanMine version : {humanMine_reactome_version_number}.")
        
      }
    }
  } else {
    enreact_use_path <- csv_up_to_date_path
  }
  if(run_update) {
    logger::log_info("Update Entities and Reactions file using Reactome API with data version {reactome_version}.")
    enreact_use_path <- createEntitiesReactionsCsv(reactome_version, dataLocation = dataLocation)
  }
  R.utils::setOption(localVars, "reactome_ent_react_file", 
                     enreact_use_path)

  
  return(enreact_use_path)
}



