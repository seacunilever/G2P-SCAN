% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPathways.R
\name{filterReactomePathways}
\alias{filterReactomePathways}
\title{Filter Pathways}
\usage{
filterReactomePathways(
  pathwayData,
  pathwayLevels = c("parental", "intermediate", "terminal"),
  onMissedUpdate = "update",
  dataLocation = "."
)
}
\arguments{
\item{pathwayData}{A data.frame of HumanMine 'Gene_Pathway' query output}

\item{pathwayLevels}{Character vector indicating what pathways levels should be returned, 'parental','intermediate' and/or 'terminal'. The default value = c('parental','intermediate','terminal') to include all pathways. Any combination of pathway levels can be used.
\itemize{
\item {parental = Pathways which have no parent pathway - the top of Reactome's pathway hierarchy}
\item {intermediate = Pathways which are not parental or terminal and have both parent and children pathways}
\item {terminal = Pathways which have no children pathways - the end of Reactome's pathway hierarchy}
}}

\item{dataLocation}{The directory path of where the 'g2pScanData' folder should be found}
}
\value{
A data.frame of HumanMine 'Gene_Pathway' query including pathway name and identifier for each gene, link to pathway Reactome site, number of geneInput in pathway and the Reactome hierarchy level (for that subset of pathways).
}
\description{
Filter InterMine pathway results for Human, Reactome pathways with optional filter based on Reactome hierarchy level
}
\details{
This function filters the results of HumanMine 'Gene_Pathway' query to get pathways for a gene input, for only Human, Reactome defined pathways. 
The remaining pathways are mapped to the data.tree of all Reactome pathways saved currently as a Rdata object found at the pathway_hierarchy_datatrees directory of the dataLocation directory (and version determined by the updateReactomeHierarchy function). The hierarchy data.tree is filtered using the pathway identifiers in the pathwayData parameter and the hierarchy of the subset of pathways is maintained.
Each pathway is then labelled with parental, intermediate and terminal levels based on their position in the subset hierarchy. The pathways can then be filtered based on this hierarchy level.
}
