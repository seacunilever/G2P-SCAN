# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Run Genes2Pathways pipeline
#'
#' @description The main wrapper function of the G2P-SCAN pipeline. This function runs the G2P-SCAN pipeline from gene input to final counts tables across model species of genes, proteins, families, entities and reactions.
#' @details The pipeline is run using the following steps: \enumerate{
#' \item{Get Pathways \itemize{\item{ 
#' This is the first step of the pipeline where input genes are mapped to pathways. Using HumanMine via InterMineR to get all Reactome pathways which any of the input genes passed to the function are are found in.
#' Pathways which have any of the input genes in are filter based on the user pathwayLevel parameter input. The Reactome pathway hierarchy is used to determine hierarchy level of the identified pathways.
#' The Reactome pathway hierarchy is read in or written to (in case of an update) to the g2pScanData/pathway_hierarchy_datatrees directory and the Reactome version matching HumanMines instance of Reactome will be used when available. }}
#' }
#' \item{Get Genes in Pathways \itemize{\item{Taking the identified pathways from step 1, all genes found in those pathways are retrieved using HumanMine's Reactome data via InterMineR queries.}}
#' }
#' \item{Get entity and reaction counts for each pathway and each species \itemize{\item{Taking the identified pathways from step 1, the number of entities and reactions for each of the pathways are identified for each chosen model species and Human.
#' Entities and Reactions are pulled using Reactome's Content and Analysis service APIs and are read in from or written to (in case of an update) the g2pScanData/Reactome_entities_reactions directory and the Reactome version matching HumanMine's instance of Reactome will be used when available.}}
#' }
#' \item{Get orthologues of pathway genes for each species \itemize{\item{For all human genes identified in step 2, orthologue genes are pulled using Panther via HumanMine/InterMineR queries.}}
#' }
#' \item{Get proteins for identified orthologues \itemize{\item{All species genes identified in step 4, are queried in UniProt to map a protein accession to each gene.
#' Only a single protein is accepted in the pipeline per gene queried and so, when multiple proteins are returned for one gene query, the following criteria is used:
#' \enumerate{
#' \item{Reviewed entries (if no reviewed entries are present, the following steps will occur on unreviewed data)}
#' \item{Strongest evidence of existence. High-low: Experimental evidence at protein level, Experimental evidence at transcript level, Protein inferred from homology, Protein predicted, Protein uncertain.}
#' \item{Longest sequence}
#' \item{Version 1 of sequence}
#' \item{The first result of the remaining data}
#' }}}
#' }
#' \item{Get Families \itemize{\item{All proteins from step 5 are queried using the InterPro API service to identify associated families to each of the proteins.
#' Only InterPro entries which are classified as a 'Family' are considered for the pipeline analysis (not 'domain' or 'superfamily'). Families returned for a single proteins are filtered to only include those which are a level 2 (first child) family for a branch in the hierarchy of InterPro families.
#' The InterPro hierarchy is read in or written to (in case of an update) to the g2pScanData/interpro_hierarchy_with_levels directory.}}
#' }s
#' \item{Count and output data \itemize{\item{Pathways - Count the number of pathways which are mapped to the input genes and meet the pathway level criteria}
#'       \item{Genes - Count the number of unique genes/orthologues in each identified pathway and per species and the number of Human genes which have orthologues mapped to them for each pathway and model species}
#'       \item{Entities & Reactions - Output entity and reation counts as returned by step 3 individually}
#'       \item{Proteins - Count the number of unique selected proteins for each identified pathway and species}
#'       \item{Families - Count the number of unique selected InterPro families or each identified pathway and species and the number of proteins per pathway and species which have no families assigned to them (unassigned)}
#'       \item{Output - Counts are combined per pathway and all pathways together and excels can be created of counts and background data (for pathways genes, protein and families)}}
#' }
#' }
#' 
#' Both the protein and Family steps (steps 5 an 6) can have queries run in parallel if a number of cores is defined. When a numeric value is passed to the cores parameter, PSOCK clusters are created to divide queries and run in parallel. If no cores are defined then queries are run sequentially. 
#' 
#' @param cores The number of cores to run in parallel. If NULL analysis will be run sequentially. 
#' @param g2pScanData The directory location of the folder "g2pScanData". If this folder does not exist state where this folder of data should be written
#' @inheritParams getPathways
#' @param outputPrefix Character vector to prefix all output files
#' @inheritParams filterReactomePathways 
#' @inheritParams getEntitiesReactionsCounts
 
#' @inheritParams getOrthologueMatrix
#' @inheritParams writeOutput
#' @param reactomeOnMissingUpdate Character string to indicate what to do when Reactome files (for pathway hierarchy and entity and reaction counts) cannot be updated to appropriate version.
#' Either "update" to use the Reactome API most up to date version. "use latest" to use the most up to date version of the file saved in the g2pScanData appropriate repository already or given an integer to state what version of Reactome file to use - this must be already present in the g2pScanData appropriate repository.  
#' @param useSynonyms Boolean value indicating when an inputted gene string is not identified as predominant identifier for a gene in HumanMine whether the identified main synonym should be used instead. TRUE will use synonyms when dominant gene symbol is found instead.
#' @include global_functions.R
#' @include global_variables.R
#' @importFrom parallel detectCores
#' @export
runGenes2Pathways <- function(inputGenes,
                              pathwayLevels,
                              species = NULL,
                              orthologueFilter = "all",
                              orthologueOutput = "both",
                              cores = (parallel::detectCores() - 1),
                              g2pScanData = ".",
                              outputDir = ".",
                              outputPrefix  = "g2pRun",
                              pathwaysTabs = TRUE,
                              countSummary = TRUE,
                              versions = TRUE,
                              inputSummary = TRUE,
                              reactomeOnMissingUpdate = "update",
                              useSynonyms = FALSE) {

  
  if(is.null(cores)){
    cl <- NULL
  } else {
    cl <- initiateCluster(cores = cores)
  }
  inputGenes <- data.frame("input" = unique(inputGenes), "synonym_replaced" = unique(inputGenes))
  foundDT <- checkGenesExists(inputGenes$input)
  if(useSynonyms & length(foundDT[foundDT$Found == 'Synonym', 'Gene']) != 0) {
    logger::log_info("useSynonym is set as TRUE therefore {foundDT[foundDT$Found == 'Synonym', 'Gene']} will be replaced with synonym {foundDT[foundDT$Found == 'Synonym', 'Synonym']} for the following analysis.")
    synonyms <- foundDT[foundDT$Found == "Synonym", ]
    for(s in synonyms$Gene){
      inputGenes[inputGenes$input == s, "synonym_replaced"] <- synonyms[synonyms$Gene == s, "Synonym"]
    }
  }
  pathways <- getPathways(inputGenes = inputGenes$synonym_replaced, pathwayLevels = pathwayLevels, onMissedUpdate = reactomeOnMissingUpdate, dataLocation = g2pScanData)
  genes_in_path <- getGenes(pathways$'Pathway Identifier')
  all_genes_symbols <- genes_in_path$`Gene Symbol`
  entities_reactions <- getEntitiesReactionsCounts(pathways = pathways$'Pathway Identifier', species = species, onMissedUpdate = reactomeOnMissingUpdate, dataLocation = g2pScanData)
  orthologues_genes <- getOrthologueGenes(genes = all_genes_symbols, species = species, orthologueFilter = orthologueFilter)
  orthologue_matrix <- getOrthologueMatrix(orthologueData = orthologues_genes, species = species, allHumanGenes = all_genes_symbols, orthologueOutput = orthologueOutput)
  orthologue_matrix_priid <- getOrthologueMatrix(orthologueData = orthologues_genes, species= species, allHumanGenes = all_genes_symbols, orthologueOutput = "primaryIdentifier")
  save_list <- list("pathways" = pathways, 
                    "genes_in_path" = genes_in_path,
                    "entities_reactions" = entities_reactions,
                    "orthologues_genes" = orthologues_genes,
                    "orthologue_matrix" = orthologue_matrix)
  save(save_list, file = "runall_save.Rdata")
  proteins <- getProteins(orthologueData = orthologues_genes, allHumanGenes = genes_in_path[, c("Gene Identifier", "Gene Symbol")], cluster = cl, stopCluster = FALSE)
  protein_matrix <- getProteinMatrix(proteinResult = proteins, orthologueData = orthologues_genes, species = species)
  families <- getFamilies(proteinResult = proteins, cluster = cl, stopCluster = TRUE, dataLocation = g2pScanData)
  family_matrix <- getFamilyMatrix(familyResult = families, proteinResult = proteins, species = species)
  pathway_counts <- getAllCounts(inputGenes = inputGenes$synonym_replaced,
                                 pathwayDT = pathways,
                                 genesInPath = genes_in_path,
                                 orthologueMatrix = orthologue_matrix_priid, 
                                 proteinDT = proteins,
                                 proteinMatrix = protein_matrix,
                                 familyDT = families,
                                 entitiesReactions = entities_reactions,
                                 species = species
  )
  writeOutput(pathwayDT = pathways,
              inputGenes = inputGenes,
              foundDT = foundDT,
              species = species,
              orthologueFilter = orthologueFilter,
              pathwayLevels = pathwayLevels,
              pathwaysCountsList = pathway_counts,
              pathwaysTabs = pathwaysTabs,
              countSummary = countSummary,
              outputDir =outputDir,
              outputFileName = paste0(outputPrefix, "_counts.xlsx"),
              versions = versions,
              inputSummary = inputSummary)
  
  writeSupplementaryData(pathwayDT = pathways,
                         inputGenes = inputGenes,
                         foundDT = foundDT,
                         species = species,
                         orthologueFilter = orthologueFilter,
                         pathwayLevels = pathwayLevels,
                         genesInPath = genes_in_path,
                         orthologueMatrix = orthologue_matrix_priid,
                         proteinMatrix = protein_matrix,
                         familyMatrix = family_matrix,
                         outputDir = outputDir,
                         outputFileName = paste0(outputPrefix, "_data.xlsx"),
                         versions = versions,
                         inputSummary = inputSummary)
  return(list("pathways" = pathways, 
              "genes_in_path" = genes_in_path,
              "entities_reactions" = entities_reactions,
              "orthologues_genes" = orthologues_genes,
              "orthologue_matrix" = orthologue_matrix,
              "proteins" = proteins,
              "protein_matrix" = protein_matrix,
              "families" = families,
              "family_matrix" = family_matrix,
              "allCounts" = pathway_counts))

}
