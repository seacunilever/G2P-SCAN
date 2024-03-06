# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Check parameter is given class
#'
#' @description Check a given value is expected class and throw fatal error if not 
#'
#' @param parameter The name of the parameter to be checked for use in the error message
#' @param value The value to have class checked - the parameter value
#' @param classType The class which the function should check for
#' @param functionName The name of the function which this function is called from for use in the error message
#' @importFrom logger log_fatal
#' @importFrom glue glue
checkClass <- function(parameter, value, classType, functionName) {
  if(class(value) != classType) {
    e <- glue::glue("'{parameter}' parameter must be of class {classType}. '{parameter}' parameter value passed to the '{functionName}' function is class {class(value)}.")
    logger::log_fatal(e)
    stop(e)
  }
}

#' Call API query and return results as a list
#' @description Function to an API query, get content and convert to a list or atomic vector of results.
#'
#' @param query API query URL
#' @param fail_on_404 TRUE/FALSE whether an 404 API error causes a log_fatal and stop().
#' @param simplifyVector Parameter of fromJSON function controlling if fromJSON function should return atomic vector instead of list. TRUE = returns atomic vector, FALSE = returns list
#' @param failonErrors TRUE/FALSE whether errors (but not 404) should cause a log_fatal and stop()
#' @importFrom jsonlite fromJSON
#' @importFrom httr RETRY content
#' @importFrom logger log_warn
#' @importFrom glue glue
#' 
#' @return List of results
callApi <- function(query, failOn404 = TRUE, failonErrors = TRUE, simplifyVector = FALSE) {
  query_result <- httr::RETRY("GET", query, encode = "application/json", )
  length_result <- length(query_result$content)
  if(query_result$status_code == 200) {
    if(length(query_result$content) != 0) { # Reactome's participants API doesnt correctly present a 404 error so this check is when the query is a bad request.
      data_output <- jsonlite::fromJSON(httr::content(query_result, "text", encoding = "UTF-8"), simplifyVector = simplifyVector)
      return(data_output)
    } else {
      message <- glue::glue("API unsuccessful with error 'No result returned'. Likely to be an unrecognised database identifier queried. Query: {query}.")
      if(failOn404) { # different logging as different function handle error different. Some functions can handle no response from a query, others it is fatal.
        
        logger::log_fatal(message)
        stop(message)
      } else {
        logger::log_warn("{message} Result value will be set as NULL.")
      }
    }
  } else if(query_result$status_code == 404){
    message <- glue::glue("API unsuccessful with status code '{query_result$status_code}'. Likely to be an unrecognised database identifier queried. Query: {query}.")
    if(failOn404) { # different logging as different function handle error different. Some functions can handle no response from a query, others it is fatal.
      
      logger::log_fatal(message)
      stop(message)
    } else {
      logger::log_warn("{message} Result value will be set as NULL.")
    }
  } else {
    message <- glue::glue("API unsuccessful with status code '{query_result$status_code}'. Failed query: {query}.")
    
    if(failonErrors) { # different logging as different function handle error different. Some functions can handle no response from a query, others it is fatal.
      logger::log_fatal(message)
      stop(message)
    } else {
      logger::log_warn("{message} Result value will be set as NULL.")
    }
  }
}

#' Build and call HumanMine query
#' @description Build and call HumanMine query by using and editing HumanMine prebuilt templates. Return results as data.frame of selected data columns.
#'
#' @details This function first calls InterMineR::getTemplateQuery using the HumanMine database and parameter templateName to get an initial template query.
#' Secondly, the template query is edited by adding custom constraint using InterMineR::setConstraints and InterMineR::setQuery and parameters constraintPath (the variable to be constrained), constraintOperators (the criteria to filter), 
#' the constraintValue (the value to filter on) and constraintm.index (the position of the constraint).
#' The HumanMine query is then run via InterMineR using edited template returning data selected by the parameter select.
#' @param templateName Name passed to getTemplateQuery
#' @param constraintPath Path passed to setConstraints - The variables which should be constrained
#' @param constraintOperators Operators passed to setConstraints - The relationship between constraint path and constraint value, eg '='
#' @param constraintValues Values passed to setConstraints - Values to be filtered on constraintPath
#' @param constraintm.index m.index passed to setConstraints - The constraint index
#' @param select Select passed to setQuery - Which values to return
#' 
#' @importFrom InterMineR getTemplateQuery setConstraints setQuery runQuery
#' @include global_variables.R
#' @return Data.frame of results from query with columns defined in select
runHumanMineQuery <- function(templateName,
                              constraintPath,
                              constraintOperators,
                              constraintValues,
                              constraintm.index,
                              select) {
  hm <- R.utils::getOption(localVars, "hm")
  template_query <- InterMineR::getTemplateQuery(im = hm, name = templateName)
  
  query_constraints <- InterMineR::setConstraints(paths = constraintPath,
                                     operators = constraintOperators,
                                     values = constraintValues,
                                     m.index = constraintm.index)
  # Create a query based on template query using geneToOrthoLDOConstraint 
  query <- InterMineR::setQuery(inheritQuery = template_query,
                    where = query_constraints,
                    select = select
  )
  results <- InterMineR::runQuery(hm, query)
  return(results)
}


#' Check passed species 
#'
#' @description Check species passed to function are within the accepted list of species.
#'
#' @param species Character vector of passed species
#' @inheritParams checkClass
#' @param includeHuman Boolean determining whether the function should allow "Human" as an input or not and whether Human should be added to output if species is NULL.
#' @importFrom logger log_fatal log_warn
#' @importFrom glue glue
#' @importFrom R.utils getOption
#' @include global_variables.R
#' @return Vector of species, if NULL is passed to the function all available species are set. 
checkSpecies <- function(species, functionName, includeHuman = FALSE) {
  
  species_name_code <- R.utils::getOption(localVars, "species_list")
  species_name <- names(species_name_code)
  if(includeHuman) {
    species_name <- c("Human", species_name)
  }
  not_species <- species[!species %in% species_name]
  if(!is.null(species)) {
    checkClass(parameter = "species", value = species, classType = "character", functionName = functionName)
    species <- species_name[species_name %in% species ] # order species
  } else {
    species <- species_name
  }
  
  if(length(not_species) > 0) {
    
    if(length(not_species) == length(species)) {
      logger::log_fatal("All variables given as species are not identified as one of {paste0(species_name, collapse = ', ')}. The following is not found: {paste0(not_species, collapse = ", ")}")
    } else {
      logger::log_warn("Some variables given as species are not identified as one of {paste0(species_name, collapse = ', ')}. The following is not found: {paste0(not_species, collapse = ", ")}")
    }
    stop(glue::glue("Unrecognised species passed to function {functionName}."))
  }
  return(species)
}

#' Check if columns are missing from data frames 
#'
#' @description Throw error if columns are missing from input data.frame
#'
#' @param dt A data.frame to check
#' @param columns The required columns to be in the data.frame
#' @param dtParameter Name of the parameter of the data.frame of the function 
#' @inheritParams checkClass
checkColumns <- function(dt, columns, dtParameter, functionName) {
  missing <- columns[!columns %in% colnames(dt)]
  if(length(missing) > 0){
    logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter {dtParameter} passed to {functionName}")
    stop(glue::glue("Missing columns in dataframe passed to {functionName}"))
  }
}

#' Initiate cluster for parallel compute
#' 
#' @description Set up PSOCK clusters and their environment readyy for use in the G2P-SCAN package
#' 
#' @param cores Number of cores to parallelise
#' @importFrom doParallel registerDoParallel 
#' @importFrom parallel makeCluster clusterExport clusterEvalQ
#' @importFrom doSNOW registerDoSNOW
#' @importFrom logger log_info
#' @importFrom R.utils getOption setOption
#' @include global_variables.R
#' @return Built cluster
#' @export
initiateCluster <- function(cores) {
  logger::log_info("Setting up clusters for parallel compute.")
  cluster <- parallel::makeCluster(cores, type = "PSOCK")
  
  rel_path_from_root_R <- file.path(R.utils::getOption(localVars, "rel_path_from_root"), "R")
  parallel::clusterExport(cl = cluster, varlist = c("rel_path_from_root_R"), envir = environment())
  parallel::clusterEvalQ(cl = cluster, expr = {
    library(InterMineR)
    library(logger)
    library(jsonlite)
    library(R.utils)
    library(glue)
    library(httr)
  })
  R.utils::setOption(localVars, "cluster", cluster)
  doSNOW::registerDoSNOW(cl = cluster)
  return(cluster)
}


#' Stop cluster and clear cluster pkgglobalenv variable to NULL
#' 
#' @importFrom parallel stopCluster 
#' @importFrom logger log_info
#' @importFrom R.utils getOption
#' @include global_variables.R
#' @export
endCluster <- function() {
  cluster <- R.utils::getOption(localVars, "cluster")
  if(!is.null(cluster)) {
    parallel::stopCluster(cl = cluster)
    logger::log_info("Cluster for parallel computing ended.")
  }
}
