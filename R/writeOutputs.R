# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


#' Create singular table for all pathways and species in a workbook sheet
#'
#' @description Add worksheet to workbook object of all pathways and species counts in one table.
#' @param wb The openxlsx workbook object to write to
#' @param pathwaysCountsList A list of data.frames of the output of getAllCounts
#' @importFrom openxlsx writeDataTable addWorksheet
#' @return openxlsx wb with new sheet appended onto it
#' @export
getCountSummaryWb <- function(wb,
                              pathwaysCountsList) {
  
  shared_columns <- c("Input.Found",
                      "Input.Found.Count",
                      "Mapped.Human.Genes",
                      "Total.Gene.Count",
                      "Protein.Count",
                      "Family.Count",
                      "Proteins.Unmapped.to.families.Count",
                      "Entity.Count",
                      "Reaction.Count")
  
  all_pathways <- lapply(names(pathwaysCountsList), FUN = function(path_id) {
    path_data <- as.data.frame(pathwaysCountsList[[path_id]])
    path_name <- path_data$Human$`Pathway Name`
    path_coverage <-  path_data$Human$`Coverage %`
    pathway_info_dt <- data.frame("Pathway ID" = path_id,
                                  "Pathway Name" = path_name,
                                  "Human_pathway_coverage_percent" = path_coverage)
    rownames(pathway_info_dt) <- path_id
    all_species <- lapply(colnames(path_data), FUN= function(sp) {
      sp_data <- path_data[, sp]
      sp_data <- lapply(sp_data, function(x) if(identical(x, character(0))| is.null(x)) NA_character_ else x)
      sp_data$`Input Found` <- paste0(sp_data$`Input Found`, collapse = ", ")
      sp_data <- as.data.frame(sp_data)
      sp_data_filtered <- sp_data[, shared_columns] 
      colnames(sp_data_filtered) <- paste0(sp, "_", colnames(sp_data_filtered))
      rownames(sp_data_filtered) <- path_id
      if(sp == "Human") {
        sp_data_filtered <- cbind(pathway_info_dt, sp_data_filtered)
      }
      return(sp_data_filtered)
    })
    
    all_species_dt <- do.call("cbind", all_species)
    
    return(all_species_dt)
  })
  all_pathways_dt <- do.call("rbind", all_pathways)
  counts_sheet <- openxlsx::addWorksheet(wb = wb, sheetName = "count summary")
  colnames(all_pathways_dt) <- gsub("[.]", " ", colnames(all_pathways_dt))
  openxlsx::writeDataTable(wb = wb, counts_sheet, x = all_pathways_dt, rowNames = F, na.string = "" )
  return(wb)
}



#' Write count workbook excel file
#'
#' @description Create and write workbook with optional tabs of input summary, database versions, count summary and tabs of counts per Human pathway.
#' @inheritParams getCountSummaryWb
#' @param foundDT The output of function checkGenesExists of data.frame of gene, found status and any synonyms found
#' @param species A character vector of species run in the analysis
#' @param pathwayLevels A character vector of the pathway levels used for the analysis
#' @param orthologueFilter A character vector of the orthologue filter used for the analysis
#' @param pathwaysTabs Boolean value. TRUE to include a tab per pathway count matrix in counts output
#' @param countSummary Boolean value. TRUE to include a table for all counts (all pathways and species) in a singular tab in the counts output
#' @param outputDir Output directory to write output file to. If directory does not exist, directory will be created at the current working directory
#' @param outputFileName What to name the output file, with '.xlsx' suffix.
#' @param versions  Boolean value. TRUE to include version tab
#' @param inputSummary Boolean value. TRUE to include a summary of the inputs of the analysis
#' @importFrom openxlsx createWorkbook saveWorkbook
#' @importFrom logger log_fatal
#' @export
writeOutput <- function(inputGenes = NULL,
                        pathwayDT = NULL,
                        foundDT = NULL,
                        species = NULL,
                        orthologueFilter = NULL,
                        pathwayLevels = NULL,
                        pathwaysCountsList,
                        pathwaysTabs = TRUE,
                        countSummary = TRUE,
                        outputDir = ".",
                        outputFileName  = "counts.xlsx",
                        versions = TRUE,
                        inputSummary = TRUE) {
  wb <- openxlsx::createWorkbook()
  if((pathwaysTabs | countSummary) & is.null(pathwaysCountsList)) {
    logger::log_fatal("No data passed to parameter pathwaysCountsList")
  }
  if(inputSummary){
    writeInputSummary(wb = wb,
                      inputGenes = inputGenes,
                      pathwayDT = pathwayDT,
                      foundDT = foundDT,
                      species = species,
                      orthologueFilter = orthologueFilter,
                      pathwayLevels =  pathwayLevels)
  }
  if(versions) {
    writeVersionWb(wb) 
  }
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }
  if(countSummary){
    getCountSummaryWb(wb = wb, pathwaysCountsList = pathwaysCountsList)
  }
  if(pathwaysTabs) {
    getAllPathwayCountsWb(wb = wb, pathwaysCountsList = pathwaysCountsList)
  }
  
  openxlsx::saveWorkbook(wb = wb, file = file.path(outputDir, outputFileName), overwrite = T)
}




#' Create supplementary worksheet
#'
#' @description Create workbook object where data behind counts are outputted, including optional genes in pathway, orthologue, proteins and families data. 
#' @param wb The openxlsx workbook object to write to
#' @param data Data to write to the workbook
#' @param sheetName The name to name the worksheet 
#' @importFrom openxlsx addWorksheet writeDataTable
#' @return openxlsx wb with new sheet appended onto it
createSupplementarySheet <- function(wb, data, sheetName) {
  sheet <- openxlsx::addWorksheet(wb = wb, sheetName = sheetName)  
  openxlsx::writeDataTable(wb = wb, sheet, x = data, rowNames = F, na.string = "" )
  return(wb)
}


#' Write supplementary data workbook excel file
#'
#' @description Create and write workbook object where data behind counts are outputted, including optional genes in pathway, orthologue, proteins and families data. Optional tabs also include input summary and database versions.
#' @inheritParams getCountsForPathway
#' @inheritParams writeOutput
#' @param familyMatrix Data.frame of output of getFamilyMatrix where human genes are as rows and columns as species,including Human, where families are listed as values.
#' @param outputDir Where to save output file. If directory does not exist, directory will be create at the current working directory.
#' @param outputFileName What to name the output file, with '.xlsx' suffix.
#' @param versions Include version tab
#' @importFrom openxlsx createWorkbook saveWorkbook
#' @export
writeSupplementaryData <- function(inputGenes = NULL,
                                   pathwayDT = NULL,
                                   foundDT = NULL,
                                   species = NULL,
                                   orthologueFilter = NULL,
                                   pathwayLevels = NULL,
                                   genesInPath = NULL,
                                   orthologueMatrix = NULL,
                                   proteinMatrix = NULL,
                                   familyMatrix = NULL,
                                   outputDir = ".",
                                   outputFileName  = "suppData.xlsx",
                                   inputSummary = TRUE,
                                   versions = TRUE) {
  wb <- openxlsx::createWorkbook()
  if(inputSummary){
    writeInputSummary(wb = wb,
                      inputGenes = inputGenes,
                      pathwayDT = pathwayDT,
                      foundDT = foundDT,
                      species = species,
                      orthologueFilter = orthologueFilter,
                      pathwayLevels =  pathwayLevels)
  }
  if(versions) {
    writeVersionWb(wb) 
  }
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }
  for(data_type in c("genesInPath", "orthologueMatrix", "proteinMatrix", "familyMatrix")) {
    
    if(!is.null(data_type)) {
      data_dt <- data.frame(get(data_type))
      if(data_type == "orthologueMatrix" | data_type == "proteinMatrix" | data_type == "familyMatrix") {
        columns <- colnames(data_dt)
        data_dt$'Human Gene' <- rownames(data_dt)
        data_dt <- data_dt[, c("Human Gene", columns)]
      }
      wb <- createSupplementarySheet(wb = wb, data = data_dt, sheetName = data_type)
    }
  }
  openxlsx::saveWorkbook(wb = wb, file = file.path(outputDir, outputFileName), overwrite = T)
}

#' Write input summary worksheet
#'
#' @description Create workbook object for record of data inputted into the analysis and info on how they mapped to pathways.
#' @inheritParams writeOutput
#' @inheritParams getPathways
#' @inheritParams getCountsForPathway
#' @inheritParams getCountSummaryWb
#' @inheritParams getOrthologueGenes
#' @importFrom openxlsx addWorksheet createStyle writeData addStyle setColWidths
#' @include global_functions.R
#' @return openxlsx wb with new sheet appended onto it
writeInputSummary <- function(wb, inputGenes, pathwayDT, foundDT, species, orthologueFilter, pathwayLevels) {
  checkClass(parameter = "inputGenes", value = inputGenes, classType = "data.frame", functionName = "writeInputSummary")
  checkColumns(dt = pathwayDT, columns = c("Pathway Identifier", "Pathways", "Derived Gene"), dtParameter = "pathwayDT", functionName = "writeInputSummary")
  checkColumns(dt = foundDT, columns = c("Gene", "Found", "Synonym"), dtParameter = "foundDT", functionName = "writeInputSummary")
  checkColumns(dt = inputGenes, columns = c("input", "synonym_replaced"), dtParameter = "inputGenes", functionName = "writeInputSummary")
  species <- checkSpecies(species = species, functionName = "writeInputSummary", includeHuman = F)
  checkClass(parameter = "pathwayLevels", value = pathwayLevels, classType = "character", functionName = "writeInputSummary")
  
  input <-  openxlsx::addWorksheet(wb = wb, sheetName = "Input Summary")
  header1 <- openxlsx::createStyle(fontSize = 14, textDecoration = c("bold"))
  header2 <- openxlsx::createStyle(fontSize = 11, textDecoration = c("bold", "underline"))
  italics <- openxlsx::createStyle(fontSize = 11, textDecoration = "italic")
  #Gene Summary
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = 1, x = "Gene Summary")
  openxlsx::addStyle(wb = wb, sheet = input, style = header1, cols = 1, rows = 1)
  ##input genes
  n_input <- length(unique(inputGenes$synonym_replaced))
  input_table <- inputGenes[, c("synonym_replaced", "input")]
  input_table[input_table$synonym_replaced == input_table$input, "input"] <- NA
  
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = 2, x = paste0(n_input, " Input Gene(s) Searched"))
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 1, rows = 2)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 2, startRow = 2, x =  "Input Gene(s) Synonyms")
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 2, rows = 2)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = 3, x = input_table, col.names = F, rowNames = F)
  
  ##mapped to pathways
  mapped_row <- n_input + 4
  mapped <- inputGenes[inputGenes$synonym_replaced %in% unique(pathwayDT$'Derived Gene'), c("synonym_replaced", "input")]
  mapped[mapped$synonym_replaced == mapped$input, "input"] <- NA
  n_mapped <- nrow(mapped)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = mapped_row, x = paste0(n_mapped, " Gene(s) Mapped To Pathways"))
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 1, rows = mapped_row)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = mapped_row +1 , x = mapped, col.names = F, rowNames = F)
  
  #not mapped
  not_mapped_row <- mapped_row + n_mapped + 2 
  not_mapped <- inputGenes[inputGenes$input %in% foundDT[foundDT$Found!=FALSE, "Gene"] & !inputGenes$synonym_replaced %in% pathwayDT$`Derived Gene`, c("synonym_replaced", "input")]
  not_mapped[not_mapped$synonym_replaced == not_mapped$input, "input"] <- NA
  n_not_mapped <- nrow(not_mapped)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = not_mapped_row, x = paste0(n_not_mapped, " Genes(s) Not Mapped To Pathways"))
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 1, rows = not_mapped_row)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = not_mapped_row +1 , x = not_mapped, col.names = F, rowNames = F)
  
  # Not Found
  not_found_row <- not_mapped_row + n_not_mapped + 2
  not_found <- foundDT[foundDT$Found == FALSE, "Gene"]
  n_not_found <- length(not_found)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = not_found_row, x = paste0(n_not_found, " Genes(s) Not Found"))
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 1, rows = not_found_row)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 1, startRow = not_found_row +1 , x = not_found, array = T)
  
  #Run selection
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = 1, x = "Run Selection")
  openxlsx::addStyle(wb = wb, sheet = input, style = header1, cols = 4, rows = 1)
  
  #Species Run
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = 2, x =  "Species Run")
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 4, rows = 2)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = 3, x = species, array = T)
  
  #Orthologues
  ortho_row <- length(species) + 4
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = ortho_row, x =  "Orthologue filter")
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 4, rows = ortho_row)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = ortho_row + 1 , x = orthologueFilter)
  
  #Pathways
  path_row <- ortho_row + 3 
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = path_row, x =  "Pathway Levels")
  openxlsx::addStyle(wb = wb, sheet = input, style = header2, cols = 4, rows = path_row)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 4, startRow = path_row + 1 , x = pathwayLevels)
  
  
  ##time
  openxlsx::writeData(wb = wb, sheet = input, startCol = 6, startRow = 1, x = "Run Date")
  openxlsx::addStyle(wb = wb, sheet = input, style = header1, cols = 6, rows = 1)
  openxlsx::writeData(wb = wb, sheet = input, startCol = 6, startRow = 2 , x = date())
  
  openxlsx::setColWidths(wb, input, cols = 1:6, widths = "auto")
  return(wb)
}
