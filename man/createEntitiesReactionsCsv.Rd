% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeEntitiesReactions.R
\name{createEntitiesReactionsCsv}
\alias{createEntitiesReactionsCsv}
\title{Creates and writes csv of the all species of analysis' entities and reactions counts for all Human pathways in Reactome pulling data from the Reactome API.}
\usage{
createEntitiesReactionsCsv(version, dataLocation = ".")
}
\arguments{
\item{version}{Version number to name output file with}

\item{dataLocation}{The directory path of where the 'g2pScanData' folder should be found}
}
\value{
File path to written file
}
\description{
Run Reactome APIs to extract entity and reaction counts per given pathway for Humans
The API will pull the data of the live version of Reactome and so this function only needs to be run when HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
This function will run Reactomes Analysis Service APIs "https://reactome.org/AnalysisService/species/homoSapiens/speciesCode?sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1" per species to get entities and reactions counts.
Where no counts exist for any of the species for a given pathway, API's "https://reactome.org/ContentService/data/participants/path_id/referenceEntities" and "https://reactome.org/ContentService/data/pathway/path_id/containedEvents" to get the Human entities and reactions counts for the missing pathways.
The function will write a matrix of entity and reaction counts for all pathways and species in csv to the dataLocation/Reactome_entities_reactions/ directory named ("Reactome_entities_reactions_v", version, ".csv").
}
