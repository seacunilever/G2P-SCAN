% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeEntitiesReactions.R
\name{getHumanEntitiesReactionsCounts}
\alias{getHumanEntitiesReactionsCounts}
\title{Get Human only Entities and Reactions counts and produce data.frame of results}
\usage{
getHumanEntitiesReactionsCounts(pathways)
}
\arguments{
\item{pathways}{List of Reactome pathway identifier of interested}
}
\value{
Human Entity Reaction data.frame matrix for given pathways only
}
\description{
Run Reactome analysis service API to extract entity and reaction counts for given Human pathways.
The API will pull the data of the live version of Reactome and so this function only needs to be run when HumanMine instance of Reactome is updated to maintain equivalent versions of Reactome across the tool.
This function is intended to be run to retrieve Human pathway entities and reaction counts for pathways which do not have data in any for the species comparison run previously. Pathways are passed to the function and for each pathway functions getHumanEntitiesCounts and getHumanReactionsCounts are run and results collated.
}
