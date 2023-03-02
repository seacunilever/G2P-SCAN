# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Filter UniProt Protein API results
#'
#' @description Function to take a single gene UniProt protein API result (as a data.frame) and filter results to obtain a singular protein accession based on the following priority choosing order:
#' \itemize{
#' \item{1. Reviewed entries (if no reviewed entries are present, the following steps will occur on unreviewed data)}
#' \item{2. Strongest evidence of existence. High-low: Experimental evidence at protein level, Experimental evidence at transcript level, Protein inferred from homology, Protein predicted, Protein uncertain.}
#' \item{3. Longest sequence}
#' \item{4. Version 1 of sequence}
#' \item{5. The first result of the remaining data}
#' }
#' After the filtering down of results, the function will output a singular protein accession number.
#' @include global_functions.R
#' @param proteinData UniProt protein API output for single gene input as a data.frame
#' @return Single UniProt protein accession 
filterUniprotProteins <- function(proteinData) {
  if((class(proteinData) == "list" & length(proteinData) == 0) | is.null(proteinData)) { # this is what is return when no results from the API are given.
    final_accession <- NA
  } else {
    checkClass("proteinData", proteinData, "data.frame", "filterUniprotProteins")
    columns_needed <- c("info.type",
                        "proteinExistence",
                        "sequence.length",
                        "sequence.version",
                        "accession")
    missing <- columns_needed[!columns_needed %in% colnames(proteinData)]
    if(length(missing) > 0){
      logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter proteinData passed to filterUniprotProteins")
      stop("Missing columns in dataframe passed to filterUniprotProteins.")
    }
    
    reviewed <- proteinData[proteinData$info.type == "Swiss-Prot", ]
    if (nrow(reviewed) >= 1){
      proteins_to_filter <- reviewed
    } else {
      proteins_to_filter <- proteinData
    }
    proteins_evidence_sort <- proteins_to_filter[order(proteins_to_filter$proteinExistence), ] #relying on alphabetical order of protein exist category remain the priority order. c('Experimental evidence at protein level','Experimental evidence at transcript level','Protein inferred by homology','Protein predicted','Protein uncertain')
    proteins_top_evidence <- proteins_evidence_sort[proteins_evidence_sort$proteinExistence == proteins_evidence_sort$proteinExistence[1],]
    max_seq_length <- max(proteins_top_evidence$sequence.length)
    max_seq_proteins <- proteins_top_evidence[proteins_top_evidence$sequence.length == max_seq_length, ]
    final_protein <- max_seq_proteins[max_seq_proteins$sequence.version == min(max_seq_proteins$sequence.version), ]
    
    final_accession <- final_protein$accession[1]
  }
  return(final_accession)
}

#' Run UniProt API
#'
#' @description Run one of 2 types of UniProt protein APIs for a singular search term. 
#' \itemize{
#' \item{1. Database search - Get UniProt entries by using UniProt cross reference(search term's database origin) and ID search term. Takes the searchID along with the species database identifier which relates to the databases Intermine uses. If searchID starts with "^ENSMUSG", the Emsembl database will be used to search.
#' https://www.ebi.ac.uk/proteins/api/proteins/database:searchID?offset=0&size=-1}
#' \item{2. Taxonomy search - Get UniProt entries by searching using a exact_gene search term and taxonomy ID. The searchID inputted in to the function will be edited to remove "Dmel_" and "CELE_" from identifier before searching. 
#' https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&exact_gene=seachID&taxid=taxid}
#' }
#' The function will filter any results to match the taxonomy of the species searched, this is because UniProt can return multiple species for one searchID.
#' @param searchID Gene identifier, symbol or primary ID, to search.
#' @param species species to search
#' @param callType whether to search using the database type or by taxonomy. 'Database', 'Taxonomy' or "Accession"
#' @importFrom logger log_fatal log_info
#' @importFrom R.utils getOption
#' @importFrom jsonlite flatten
#' @include global_functions.R
#' @return Result of the API call as data.frame
#' @export
callUniprotProteinAPI <- function(searchID, species, callType) {
  logger::log_info("searchID: {searchID} ")
  checkClass(parameter = "searchID", value = searchID, classType = "character", functionName = "callUniprotProteinAPI")
  if(is.null(species)) {
    logger::log_fatal("NULL is passed to function callUniprotProteinAPI where species is needed.")
    stop("NULL is passed to function callUniprotProteinAPI where species is needed.")
  }
  species <- checkSpecies(species = species, functionName = "callUniprotProteinAPI", includeHuman = TRUE)
  if(!any(tolower(callType) %in% c("database", "taxonomy"))) {
    logger::log_fatal("callType passed to callUniprotProteinAPI is neither 'database' or 'taxonomy'. Value passed: {callType}.")
    stop("callType passed to getOrthologueGenes is neither 'database' or 'taxonomy'.")
  }
  if(species == "Human") {
    species_name_code <- R.utils::getOption(localVars, "human_codes")
  } else {
    species_name_code <- R.utils::getOption(localVars, "species_list")
  }
  tax <- species_name_code[[species]]$taxonomy
  if(tolower(callType) == "database") {
    if(grepl("^ENSMUSG", searchID)) { # For Human gene identifier from HumanMine, id very from HumanMine ids and Ensembl ID. When "ENSMUSG" is present in the term then Ensembl database needs to be used.
      db <- "Ensembl"
    } else {
      db <- species_name_code[[species]]$interMineDB
    }
    requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/proteins/", db,":",gsub(":", "%3A", searchID),"?offset=0&size=-1")
  } else if(tolower(callType) == "taxonomy") {
    id_edit <- gsub("Dmel_", "", searchID)
    id_edit <- gsub("CELE_", "", id_edit)
    
    requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&exact_gene=", id_edit,"&taxid=", tax)
  } else if(tolower(callType) == "accession") {
    requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=", searchID)
  }
  api_result <- callApi(query = requestURL, failOn404 = FALSE, simplifyVector = TRUE)
  if(length(api_result) != 0) {
    result_tax <- jsonlite::flatten(api_result[api_result$organism$taxonomy == tax, ])
  } else {
    result_tax <- api_result
  }
  return(result_tax)
}

#' Run a batch of search terms through the UniProt protein APIs and filter protein results to one protein accession per search term. 
#' 
#' @description This function can run one of two types of UniProt protein APIs for multiple search terms in batch and parallel. Results will be filtered on a priority criteria to result in one protein accession per search term.
#' UniProt protein API calls (callType)
#' \itemize{
#' \item{1. Database search - Get UniProt entries by using UniProt cross reference(search term's database origin) and ID search term. Takes the searchID along with the species database identifier which relates to the databases Intermine uses. If searchID starts with "^ENSMUSG", the Emsembl database will be used to search.
#' https://www.ebi.ac.uk/proteins/api/proteins/database:searchID?offset=0&size=-1}
#' \item{2. Taxonomy search - Get UniProt entries by searching using a exact_gene search term and taxonomy ID. The searchID inputted in to the function will be edited to remove "Dmel_" and "CELE_" from identifier before searching. 
#' https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&exact_gene=seachID&taxid=taxid}
#' }
#' Results of the API call are initially filtered to match the taxonomy ID of the species searched before further filtering per gene result using the priority criteria for protein accession selection below:
#' \itemize{
#' \item{1. Reviewed entries (if no reviewed entries are present, the following steps will occur on unreviewed data)}
#' \item{2. Strongest evidence of existence. High-low: Experimental evidence at protein level, Experimental evidence at transcript level, Protein inferred from homology, Protein predicted, Protein uncertain.}
#' \item{3. Longest sequence}
#' \item{4. Version 1 of sequence}
#' \item{5. The first result of the remaining data}
#' }
#' The input data.frame to the function will be appended with the column "protein.accession" where the final selected protein accession number per row is given. The "search_term" columns is removed in the output.
#' 
#' Running of the API calls and filtering of results are run sequentially but each search terms run in parallel.
#' 
#' @param dataToSearch Data.frame of data to search, must include column "search_term".
#' @param cluster Cluster to run parallel operation on
#' @param stopCluster If TRUE 'parallel::stopCluster(cl = cluster)' is run and localVar variable for cluster is set to NULL. If FALSE the cluster is not stopped
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel makeCluster stopCluster
#' @importFrom logger log_fatal
#' @importFrom rprojroot find_root_file
#' @include global_functions.R
#' @return return dataframe including protein accession results
runAllProteinSearch <- function(dataToSearch, callType, cluster = NULL, stopCluster = TRUE) {
  checkClass(parameter = "dataToSearch", value = dataToSearch, classType = "data.frame", functionName = "runAllProteinSearch")
  columns_needed <- c("search_term",
                      "Gene.homologues.homologue.organism.shortName"
  )
  missing <- columns_needed[!columns_needed %in% colnames(dataToSearch)]
  if(length(missing) > 0){
    logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter dataToSearch passed to runAllProteinSearch")
    stop("Missing columns in dataframe passed to runAllProteinSearch")
  }
  checkSpecies(species = unique(dataToSearch$"Gene.homologues.homologue.organism.shortName"), functionName = "runAllProteinSearch", includeHuman = TRUE)
  database_search_results <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(database_search_results) <- c("Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.homologue.symbol" , "Gene.homologues.homologue.organism.shortName", "protein.accession")
  search_terms <- dataToSearch$search_term
  
  if(!is.null(cluster)){
    all_api_results <- foreach::foreach(X = search_terms) %dopar% {
      callUniprotProteinAPI(searchID = X, species = dataToSearch[dataToSearch$search_term == X, "Gene.homologues.homologue.organism.shortName"], callType = callType)
    }
    protein_accessions <- foreach::foreach(X = all_api_results, .combine=c ) %dopar% {
      filterUniprotProteins(proteinData = X)
    }
    if(stopCluster) {
      endCluster()
    }
  } else {
    all_api_results <- lapply(X = search_terms,
                              FUN = function(X) {callUniprotProteinAPI(searchID = X,
                                                                       species = dataToSearch[dataToSearch$search_term == X, "Gene.homologues.homologue.organism.shortName"],
                                                                       callType = callType)})
    protein_accessions <- unlist(lapply(X = all_api_results, FUN = function(X){
      filterUniprotProteins(proteinData = X) }
    ))
  }
  database_search_results <- dataToSearch
  database_search_results$'protein.accession' <- protein_accessions
  
  return(database_search_results)
}

#' Get Protein Accessions for gene inputs
#'
#' @description For every gene row passed to this function in the orthologueData and allHumanGenes parameters, UniProt APIs will be called to determine a protein accession for each gene. 
#' Steps: \itemize{
#' \item {1. Assign primary IDs of input data as the search term. When no primary ID is provided but gene symbol is, gene symbol is assigned as the search term. }
#' \item {2. Run UniProt's protein API using database identifiers. UniProt entries are retrieved using UniProt cross reference and its ID. This step takes the searchID along with the species database identifier which relates to the databases Intermine uses. If searchID starts with "^ENSMUSG", the Emsembl database will be used to search.
#' API URL = https://www.ebi.ac.uk/proteins/api/proteins/database:searchID?offset=0&size=-1}
#' \item {3. Filter results to retrieve a singular protein accession for the gene search term using prioritisation criteria below
#' \itemize{
#' \item{1. Reviewed entries (if no reviewed entries are present, the following steps will occur on unreviewed data)}
#' \item{2. Strongest evidence of existence. High-low: Experimental evidence at protein level, Experimental evidence at transcript level, Protein inferred from homology, Protein predicted, Protein uncertain.}
#' \item{3. Longest sequence}
#' \item{4. Version 1 of sequence}
#' \item{5. The first result of the remaining data}
#' }}
#' \item {4. For entries which no protein accession were retrieved at step 3, run UniProt's protein API using taxonomy. When a gene symbol is available the search term is assigned to the gene symbol, when not available primary ID is used. UniProt entries are retrieved by searching using a exact_gene search term and  the taxonomy ID of the species. The search term will be edited to remove "Dmel_" and "CELE_" from identifier before searching. 
#' API URL = https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&exact_gene=searchID&taxid=taxid#'}
#' \item {5. Combine results, returning data.frame of following columns: \itemize{
#' \item{Gene.homologues.homologue.primaryIdentifier}
#' \item{Gene.homologues.homologue.symbol}
#' \item{Gene.homologues.homologue.organism.shortName}
#' \item{protein.accession}
#' }
#' }
#' }
#' 
#' If a cluster is passed to the function the running of the API calls and filtering of results are run sequentially but each search terms run in parallel.
#' 
#' @param orthologueData Data from output of function getOrthologueGenes
#' @param allHumanGenes Data.frame of "Gene.primaryIdentifier" and "Gene.symbol" as columns in that order. Function is not sensitive to column names but expects primary IDs in the first columns and gene symbols in the second
#' @inheritParams runAllProteinSearch
#' @importFrom logger log_info 
#' @importFrom R.utils getOption
#' @importFrom parallel mclapply
#' @include global_functions.R
#' @return Data.frame of input data appended with protein.accession results
#' @export
getProteins <- function(orthologueData, allHumanGenes, cluster = NULL, stopCluster = TRUE) {
  checkOrthologueDataFrame(orthologueData, functionName = "getProteins")
  checkClass(parameter = "allHumanGenes", value = allHumanGenes, classType = "data.frame", functionName = "getProteins")
  allHumanGenes$Gene.homologues.homologue.organism.shortName <- "Human"
  allHumanGenes$Gene.homologues.homologue.symbol <- allHumanGenes[, 2]
  colnames(allHumanGenes) <- c("Gene.homologues.homologue.primaryIdentifier", "Gene.symbol", "Gene.homologues.homologue.organism.shortName", "Gene.homologues.homologue.symbol" )
  all_data_to_search <- rbind(orthologueData[, c("Gene.symbol", "Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.homologue.organism.shortName", "Gene.homologues.homologue.symbol")],
                              allHumanGenes[,c("Gene.symbol", "Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.homologue.organism.shortName", "Gene.homologues.homologue.symbol")])
  all_data_to_search$search_term <- all_data_to_search$Gene.homologues.homologue.primaryIdentifier
  all_data_to_search[is.na(all_data_to_search$search_term), "search_term"] <- all_data_to_search[is.na(all_data_to_search$search_term), "Gene.homologues.homologue.symbol"]
  all_data_to_search <- all_data_to_search[rownames(unique(all_data_to_search[, c("Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.homologue.organism.shortName", "Gene.homologues.homologue.symbol")])),]
  
  logger::log_info("Running UniProt protein API using database identifiers for {nrow(all_data_to_search)} genes identifiers.")
  database_search_results <- runAllProteinSearch(dataToSearch = all_data_to_search, callType = "database", cluster = cluster, stopCluster = stopCluster)
  
  not_found_to_search <- database_search_results[is.na(database_search_results$protein.accession), c("Gene.symbol", "Gene.homologues.homologue.primaryIdentifier", "Gene.homologues.homologue.symbol", "Gene.homologues.homologue.organism.shortName")]
  not_found_to_search[not_found_to_search$Gene.homologues.homologue.symbol != "Gene symbol not found", "search_term"] <- not_found_to_search[not_found_to_search$Gene.homologues.homologue.symbol != "Gene symbol not found", "Gene.homologues.homologue.symbol"]
  not_found_to_search[not_found_to_search$Gene.homologues.homologue.symbol == "Gene symbol not found", "search_term"] <- not_found_to_search[not_found_to_search$Gene.homologues.homologue.symbol == "Gene symbol not found", "Gene.homologues.homologue.primaryIdentifier"]
  
  logger::log_info("Running UniProt protein API using taxonomy IDs for {nrow(not_found_to_search)} genes identifiers.")
  taxonomy_search_results <- runAllProteinSearch(not_found_to_search, callType = "taxonomy", cluster = cluster, stopCluster = stopCluster)
  
  database_search_results_fitlered <- database_search_results[!database_search_results$Gene.homologues.homologue.primaryIdentifier %in% not_found_to_search$Gene.homologues.homologue.primaryIdentifier, ]
  final_protein_accessions <- rbind(database_search_results_fitlered, taxonomy_search_results)
  na_num <- nrow(final_protein_accessions[is.na(final_protein_accessions$"protein.accession"), ])
  assigned_num <- nrow(final_protein_accessions[!is.na(final_protein_accessions$"protein.accession"), ])
  logger::log_info("Protein accessions assigned for {assigned_num} gene identifier. No protein accessions assigned for {na_num} gene identifiers.")
  return(final_protein_accessions)
}

#' Get matrix of proteins
#'
#' @description Create data.frame matrix of human genes as rows and columns as species, including human, where protein accessions are listed as values.
#' @param proteinResult Data from output of function getProteins
#' @inheritParams getProteins
#' @inheritParams runAllProteinSearch
#' @inheritParams getOrthologueGenes
#' @importFrom logger log_fatal
#' @importFrom R.utils getOption
#' @include global_functions.R
#' @return data.frame matrix of human genes as rows and columns as species,including Human, where proteins are listed as values.
#' @export
getProteinMatrix <- function(proteinResult, orthologueData, species) {
  checkOrthologueDataFrame(orthologueData, functionName = "getProteinMatrix")
  checkClass(parameter = "proteinResult", value = proteinResult, classType = "data.frame", functionName = "getProteinMatrix")
  columns_needed <- c("protein.accession",
                      "Gene.symbol",
                      "Gene.homologues.homologue.symbol",
                      "Gene.homologues.homologue.organism.shortName"
  )
  if(!"Human" %in% proteinResult$Gene.homologues.homologue.organism.shortName) {
    logger::log_fatal("No Human data found in proteinResult dataframe passed to getProteinMatrix")
    stop("No Human data found in proteinResult dataframe passed to getProteinMatrix")
  }
  missing <- columns_needed[!columns_needed %in% colnames(proteinResult)]
  if(length(missing) > 0){
    logger::log_fatal("Columns {paste0(missing, collapse = ', ')} are missing from parameter proteinResult passed to getProteinMatrix")
    stop("Missing columns in dataframe passed to getProteinMatrix")
  }
  
  protein_ortho_merged <- merge(x = orthologueData, y = proteinResult[, c("Gene.homologues.homologue.primaryIdentifier",
                                                                          "Gene.homologues.homologue.organism.shortName",
                                                                          "Gene.homologues.homologue.symbol",
                                                                          "protein.accession")],
                                by = c("Gene.homologues.homologue.primaryIdentifier",
                                       "Gene.homologues.homologue.organism.shortName",
                                       "Gene.homologues.homologue.symbol"),
                                all.x = TRUE)
  human_genes <- unique(protein_ortho_merged$Gene.symbol)
  
  proteins_sorted <- protein_ortho_merged[order(proteinResult$Gene.symbol), ] 
  # get  input species and maintain order of species.
  species_ordered <- c("Human", checkSpecies(species = species, functionName = "getProteinMatrix", includeHuman = F))
  # create empty matrix
  protein_matrix <- data.frame(matrix(nrow = length(human_genes), ncol = length(species_ordered)))
  rownames(protein_matrix) <- sort(human_genes)
  colnames(protein_matrix) <- species_ordered
  # fill matrix
  for(sp in species_ordered) {
    if(sp == "Human"){
      species_results <- proteinResult[proteinResult$Gene.homologues.homologue.organism.shortName == sp, ]
    } else {
      species_results <- protein_ortho_merged[protein_ortho_merged$Gene.homologues.homologue.organism.shortName == sp, ]
    }
    # concat when there are more than one protein for 1 human gene
    species_results[,"protein.accession.string"] <- sapply(X= species_results$Gene.symbol,
                                                           FUN = function(p, species_results){
                                                             proteins <- species_results[species_results$Gene.symbol == p, "protein.accession"]
                                                             if(length(proteins)>0) {
                                                               output <- paste0(proteins, collapse = ", ")
                                                             } else {
                                                               output <- NA
                                                             }
                                                             output
                                                           }, species_results)
    species_human_genes <- species_results$Gene.symbol
    protein_matrix[species_human_genes, sp] <- species_results$protein.accession.string
  }
  return(protein_matrix)
}