# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.


# These are variables which are used throughout the package functions. Everything within the Options are treated as global variables. 
# Create local variable where InterMineR's HumanMine instance can be 'saved'. This because hm is used to define the species database used for Intermine queries in most functions.
localVars <- R.utils::Options()
hm <-
  R.utils::setOption(localVars, "hm", InterMineR::initInterMine(mine = listMines()["HumanMine"]))


# Add more species here.
species_list <-
  R.utils::setOption(localVars, "species_list", species_list <- list(
    "R. norvegicus" = list(reactomeID = 48895, taxonomy = 10116, interMineDB = "RGD"),
    "M. musculus" = list(reactomeID = 48892,taxonomy = 10090, interMineDB = "MGI"),
    "D. rerio" = list(reactomeID = 68323, taxonomy = 7955, interMineDB = "ZFIN"),
    "C. elegans" = list(reactomeID = 68320, taxonomy = 6239, interMineDB = "WormBase"),
    "D. melanogaster" = list(reactomeID = 56210, taxonomy = 7227, interMineDB = "FlyBase"),
    "S. cerevisiae" = list(reactomeID = 68322, taxonomy = 559292, interMineDB = "SGD") # code is reactome species codes
  ))

human_codes <- R.utils::setOption(localVars, "human_codes", human_codes <- list(
  "Human" = list(taxonomy = 9606, interMineDB = "GeneID")
))

rel_path_from_root <- R.utils::setOption(localVars, "rel_path_from_root",
                                         rel_path_from_root <- rprojroot::find_root_file(criterion = rprojroot::has_file("DESCRIPTION")))


