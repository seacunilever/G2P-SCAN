# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Map pathways to Reactome pathway hierarchy
#'
#' @description Map Reactome pathways to the Reactome hierarchy to get a hierarchy level, 'parental', 'intermediate' or 'terminal', for each pathways passed to the function.
#'
#' @details
#' Pathways are mapped to a data.tree of all Reactome pathways using saved hierarchy data.tree Rdata file of Reactome hierarchy structure of all pathways. The data.tree is subset using the pathway identifiers in the pathwayData parameter and the hierarchy of the subset of pathways is maintained.
#' Each pathway is then labelled with parental, intermediate and terminal levels based on their position in the subset hierarchy. Parental level pathways are defined as those which have children in the hierarchy of selected pathways but do not have any parent themselves.
#' Terminal level pathways are defined as those which have parent pathways, but no children and intermediate pathways have both parent and children pathways. Hierarchy level is determined based on the filter hierarchy of pathways.
#' 
#' The Rdata file can be found in the pathway_hierarchy_datatrees directory of the dataLocation directory. The file with the matching version to HumanMine version of Reactome is used if exists and if not use or creation of the file is based on theMissedUpdate parameter and updateReactomeHierarchy
#' @param pathways A dataframe including  Reactome pathway "Pathway Identifiers" of which hierarchy levels should returned
#' @inheritParams updateReactomeHierarchy 
#' @inheritParams filterInterProResults
#' @importFrom R.utils getOption
#' @importFrom data.tree Prune FindNode
#' @include global_functions.R
#' @include updateReactomeResource.R
#' 
#' @return Pathways dataframe with an appended "Pathway Level" column with 'parental', 'intermediate' or 'terminal' labels for each pathway
getPathwayLevels <- function(pathways, onMissedUpdate = "update", dataLocation = ".") {
  checkClass(parameter = "pathways", value = pathways, classType = "data.frame", functionName = "getPathwayLevels")
  checkColumns(dt = pathways, columns = "Pathway Identifier", dtParameter = "pathways", functionName = "getPathwayLevels")
  # Build Pathway Hierarchy for the list of pathways and filter based on hierarchy level.
  updateReactomeHierarchy(onMissedUpdate = onMissedUpdate, dataLocation = dataLocation) # set localVars reactome_data.tree to get path to right hierarchy datatree - update if needed.
  hierarchy_tree_path <- R.utils::getOption(localVars, "reactome_datatree")
  pathTree <-
    get(load(hierarchy_tree_path)) # load result of pathway hierarchy conversion to data.tree.
  data.tree::Prune(pathTree, function(x) {
    x$name %in% pathways$"Pathway Identifier"
  }) # Prune hierarchy to only pathways from pathwaysList
  
  for(p in pathways$"Pathway Identifier") {
    if (length(data.tree::FindNode(pathTree, p)$children) == 0 ||
        is.null(data.tree::FindNode(pathTree, p)$children)) {
      # terminal pathways have no children
      level <- "terminal" 
    } else if (data.tree::FindNode(pathTree, p)$level == 2) {
      # as the root node is "top", parental pathways have a level of 2 (not 1)
      level  <- "parental"
    } else {
      level  <- "intermediate"
    }
    pathways[pathways$"Pathway Identifier" == p, "Pathway Level"] <-
      level
  }
  
  return(pathways)
}


#' Filter Pathways
#'
#' @description Filter InterMine pathway results for Human, Reactome pathways with optional filter based on Reactome hierarchy level
#'
#' @details
#' This function filters the results of HumanMine 'Gene_Pathway' query to get pathways for a gene input, for only Human, Reactome defined pathways. 
#' The remaining pathways are mapped to the data.tree of all Reactome pathways saved currently as a Rdata object found at the pathway_hierarchy_datatrees directory of the dataLocation directory (and version determined by the updateReactomeHierarchy function). The hierarchy data.tree is filtered using the pathway identifiers in the pathwayData parameter and the hierarchy of the subset of pathways is maintained.
#' Each pathway is then labelled with parental, intermediate and terminal levels based on their position in the subset hierarchy. The pathways can then be filtered based on this hierarchy level.
#'
#' @param pathwayData A data.frame of HumanMine 'Gene_Pathway' query output
#' @param pathwayLevels Character vector indicating what pathways levels should be returned, 'parental','intermediate' and/or 'terminal'. The default value = c('parental','intermediate','terminal') to include all pathways. Any combination of pathway levels can be used.
#' \itemize{
#' \item {parental = Pathways which have no parent pathway - the top of Reactome's pathway hierarchy}
#' \item {intermediate = Pathways which are not parental or terminal and have both parent and children pathways}
#' \item {terminal = Pathways which have no children pathways - the end of Reactome's pathway hierarchy}
#' }
#' @inheritParams filterInterProResults
#' @include global_functions.R
#' @return A data.frame of HumanMine 'Gene_Pathway' query including pathway name and identifier for each gene, link to pathway Reactome site, number of geneInput in pathway and the Reactome hierarchy level (for that subset of pathways).
filterReactomePathways <-
  function(pathwayData,
           pathwayLevels = c('parental', 'intermediate', 'terminal'),
           onMissedUpdate = "update",
           dataLocation = ".") {
    checkClass(parameter = "pathwayData", value = pathwayData, classType = "data.frame", functionName = "filterReactomePathways")
    checkClass(parameter = "pathwayLevels", value = pathwayLevels, classType = "character", functionName = "filterReactomePathways")
    
    # Filter the results to only contain human pathways and columns "Gene Name", "Pathway Name", "Data Sets Name" and "Pathways Identifier"
    pathwaysList <-
      unique(pathwayData[pathwayData$"Gene.organism.shortName" == "H. sapiens", !colnames(pathwayData) %in% c("Gene.organism.shortName")])
    if(!is.null(pathwaysList)) {
      # Rename pathwaysList column names
      colnames(pathwaysList) <-
        c("Derived Gene",
          "Pathways",
          "Derived Database",
          "Pathway Identifier")
      
      pathwaysList <-
        pathwaysList[pathwaysList$"Derived Database" == "reactome", ] # Filter to Reactome data only
      pathways_hierarchy <- getPathwayLevels(pathwaysList, onMissedUpdate = onMissedUpdate, dataLocation = dataLocation)
      final_pathways <- pathways_hierarchy[pathways_hierarchy$"Pathway Level" %in% pathwayLevels, ]
    } else {
      final_pathways <- pathwaysList # setting in the case pathwayList is null or pathwayLevels is c('parental', 'intermediate', 'terminal')
    }
    
    return(final_pathways)
  }


#' Get Pathways
#'
#' @description Get Reactome pathways which gene input are found in (using HumanMine via InterMineR package) and filter human Reactome pathways based on hierarchy level.
#'
#' @details
#' This function takes a vector of gene symbols and calls HumanMine's template query 'Gene_Pathway' which returns all pathways genes in the vector belong to.
#' Pathways results from HumanMine are  filtered for only Human Reactome defined pathways. (Reactome verison used is dependent on verison loaded in HumanMine - see https://www.humanmine.org/humanmine/dataCategories.do).
#' Pathways are further filtered for only pathways which are the same hierarchy level as set in the pathwayLevels parameter.
#' The hierarchy level is determined by a subset of the hierarchy data.tree created where the hierarchy levels are as follows:
#' \itemize{
#' \item {parental = Pathways which have no parent pathway - the top of Reactome's pathway hierarchy}
#' \item {intermediate = Pathways which are not parental or terminal and have both parent and children pathways}
#' \item {terminal = Pathways which have no children pathways - the end of Reactome's pathway hierarchy}
#' }
#' The data.tree of all Reactome pathways is read from or saved to a Rdata object found at the pathway_hierarchy_datatrees directory of the dataLocation directory and version of the Reactome hierarchy used is determined by the updateReactomeHierarchy function and onMissedUpdate parameter. The Reactome hierarchy of the same version of Reactome that HumanMine uses will be used when possible.
#' @param inputGenes Character vector of gene symbols to be searched
#' @inheritParams updateReactomeHierarchy 
#' @inheritParams filterReactomePathways 
#' @inheritParams filterInterProResults
#' @importFrom logger log_fatal log_warn
#' @include global_functions.R
#' @include checkGeneExists.R
#' @return A data.frame of HumanMine 'Gene_Pathway' query including pathway name and identifier for each gene, link to pathway Reactome site, number of geneInput in pathway and the Reactome hierarchy level (for that subset of pathways).
#' @export
getPathways <-
  function(inputGenes,
           pathwayLevels = c('parental', 'intermediate', 'terminal'),
           onMissedUpdate = "update",
           dataLocation = ".") {
    # Error Handling
    checkClass(parameter = "inputGenes", value = inputGenes, classType = "character", functionName = "getPathways")
    checkClass(parameter = "pathwayLevels", value = pathwayLevels, classType = "character", functionName = "getPathways")
    if(any(!pathwayLevels %in% c('parental', 'intermediate', 'terminal'))) {
      oddOne <- pathwayLevels[!pathwayLevels %in% c('parental', 'intermediate', 'terminal')]
      logger::log_fatal("pathwayLevels parameter for function getPathways includes invalid input: '{oddOne}'. Must be only 'parental', 'intermediate' or 'terminal'")
      stop("pathwayLevels parameter for function getPathways includes invalid input. Must be only 'parental', 'intermediate' or 'terminal'")
    }
    gene_to_path_result <- runHumanMineQuery(templateName = "Gene_Pathway",
                                             constraintPath = "Gene.symbol",
                                             constraintOperators = "=",
                                             constraintValues = list(c(toupper(inputGenes),
                                                                       tolower(inputGenes))),
                                             select = c("Gene.symbol",
                                                        "Gene.pathways.name",
                                                        "Gene.pathways.dataSets.dataSource.dataSets.name",
                                                        "Gene.pathways.identifier",
                                                        "Gene.organism.shortName"
                                             )
    )
    
    if(is.null(gene_to_path_result)) {
      gene_not_found_reason <- checkGenesExists(inputGenes)
      stop("No Human Reactome pathways returned for gene input.")
    } else {
      final_pathways <- filterReactomePathways(pathwayData = gene_to_path_result, pathwayLevels = pathwayLevels, dataLocation = dataLocation, onMissedUpdate = onMissedUpdate)
      
      # Warn when a gene is not found in any pathways
      
      not_found <- inputGenes[! toupper(inputGenes) %in% gene_to_path_result$Gene.symbol]# Not found in search of pathways
      if(length(not_found) > 0) {
        logger::log_info("Gene symbol(s) {paste0(not_found,collapse = ', ' )} have not be found to be mapped to any pathway. Checking if these genes are found as recognised gene symbols in HumanMine...")
        gene_not_found_reason <- checkGenesExists(not_found)
        logger::log_info("If gene symbol(s) are found then they just do not have any pathway mapped to them in the Reactome database.")
      }
      
      filtered_out <- inputGenes[(toupper(inputGenes) %in% gene_to_path_result$Gene.symbol) & (!toupper(inputGenes) %in% final_pathways$`Derived Gene`)]# Not found in search of pathways
      if(length(filtered_out) > 0) {
        logger::log_info("Gene symbol(s) {paste0(filtered_out, collapse = ', ')} were mapped to pathways but after hierarchical filter these pathways are lost. In the remianing pathway list these genes are not found.")
      }
      if(nrow(final_pathways)>0) {
        # Create Reactome link
        final_pathways$"Pathway Link" <-
          paste0("https://reactome.org/content/detail/",
                 final_pathways$"Pathway Identifier")
        
        # Count the number of times the pathways is matched with the given input genes (count the number of appearances of that pathways in the table)
        for (i in 1:nrow(final_pathways)) {
          pathway <- final_pathways[i, "Pathway Identifier"]
          final_pathways[i, "Input Gene Count"] <-
            nrow(final_pathways[final_pathways$"Pathway Identifier" == pathway, ])
        }
      } else { # end analysis is there not pathways
        logger::log_fatal("No pathways returned using gene input ({paste0(inputGenes, collapse = ',')}) and pathway filtering ({paste0(pathwayLevels, collapse = ',')}). Analysis ending.")
        stop("No pathways returned")
        
      }
    }
    
    return(final_pathways)
  }
