fake_proteins1 <- data.frame("info.type" = c("TrEMBL", "TrEMBL", "TrEMBL", "TrEMBL", "TrEMBL"),
                                            "sequence.version" = c(1,1,1,1,2),
                                            "sequence.length" = c(505, 353, 222, 547, 547),
                                            "proteinExistence" = c("Evidence at protein level", "Evidence at transcript level", "Evidence at protein level", "Evidence at transcript level", "Evidence at protein level"),
                                      "accession" = c("B1PS59", "Q9I9P4", "B3DHM7", "B1PS58", "F1QS48"))

fake_proteins2 <- data.frame("info.type" = "Swiss-Prot",
                             "sequence.version" = 1,
                             "sequence.length" = 298,
                             "proteinExistence" = "Evidence at protein level",
                             "accession" = "P00546")

fake_proteins3 <- data.frame("info.type" = c("Swiss-Prot", "TrEMBL", "TrEMBL", "TrEMBL", "TrEMBL", "TrEMBL", "TrEMBL"),
                             "sequence.version" = c(4, 1, 1, 1, 1, 1, 1),
                             "sequence.length" = c(393, 158, 393, 354, 234, 182, 187),
                             "proteinExistence" = c("Evidence at protein level", "Evidence at transcript level", "Evidence at transcript level", "Evidence at transcript level", "Evidence at protein level", "Evidence at protein level", "Evidence at protein level"),
                             "accession" = c("P04637", "Q53GA5", "K7PPA8", "H2EHT1", "A0A087X1Q1", "A0A087WXZ1", "A0A087WT22"))

fake_proteins4 <- data.frame("info.type" = c("Swiss-Prot", "TrEMBL", "Swiss-Prot", "TrEMBL"),
                            "sequence.version" = c(1, 1, 1, 1),
                            "sequence.length" = c(291, 158, 247, 50),
                            "proteinExistence" = c("Evidence at protein level", "Evidence at protein level", "Evidence at protein level", "Evidence at protein level"),
                            "accession" = c("P46737", "A3KGA8", "P46737-2", "E9Q0P6"))
fake_filter_input1 <- data.frame("wrong" = c("Swiss-Prot", "TrEMBL", "Swiss-Prot", "TrEMBL"),
                             "sequence.version" = c(1, 1, 1, 1),
                             "sequence.length" = c(291, 158, 247, 50),
                             "proteinExistence" = c("Evidence at protein level", "Evidence at protein level", "Evidence at protein level", "Evidence at protein level"),
                             "accession" = c("P46737", "A3KGA8", "P46737-2", "E9Q0P6"))
fake_filter_input2 <- list("sequence.length" = c(291, 158, 247, 50))
test_that("Test protein filter results", {
  #needs all filters
  accession1 <- filterUniprotProteins(fake_proteins1)
  # only 1 result
  accession2 <- filterUniprotProteins(fake_proteins2)
  # only 1 reviewed
  accession3 <- filterUniprotProteins(fake_proteins3)
  # 2 reviewed choice between sequence length
  accession4 <- filterUniprotProteins(fake_proteins4)
  accession5 <- filterUniprotProteins(list())
  expect_equal(accession1, "F1QS48")
  expect_equal(accession2, "P00546")
  expect_equal(accession3, "P04637")
  expect_equal(accession4, "P46737")
  expect_equal(accession5, NA)
})

test_that("Test protein filter input", {
  #missing column
  expect_error(error1 <- filterUniprotProteins(fake_filter_input1),
               "Missing columns in dataframe passed to filterUniprotProteins.")
  expect_error(error2 <- filterUniprotProteins("wrong"),
               "'proteinData' parameter must be of class data.frame. 'proteinData' parameter value passed to the 'filterUniprotProteins' function is class character.")
  })

test_that("Test callUniprotProteinAPI input", {
  expect_error(error1 <- callUniprotProteinAPI(searchID = 1, species = "Human", callType = "database"),
               "'searchID' parameter must be of class character. 'searchID' parameter value passed to the 'callUniprotProteinAPI' function is class numeric.")
  expect_error(error2 <- callUniprotProteinAPI(searchID = "ENSMUSG00000031201", species = "wrong", callType = "Database"),
               "Unrecognised species passed to function callUniprotProteinAPI.")
  expect_error(error3 <- callUniprotProteinAPI(searchID = "ENSMUSG00000031201", species = "Human", callType = "wrong"),
               "callType passed to getOrthologueGenes is neither 'database' or 'taxonomy'")
  })

test_that("Test runAllProteinSearch input", {
  expect_error(error1 <- runAllProteinSearch(dataToSearch = "wrong", callType = "database"),
               "'dataToSearch' parameter must be of class data.frame. 'dataToSearch' parameter value passed to the 'runAllProteinSearch' function is class character.")
  expect_error(error2 <- runAllProteinSearch(dataToSearch = data.frame("wrong" = "wrong", "Gene.homologues.homologue.organism.shortName" = "Human"), callType = "database"),
               "Missing columns in dataframe passed to runAllProteinSearch")
  expect_error(error3 <- runAllProteinSearch(dataToSearch = data.frame("search_term" = "test", "Gene.homologues.homologue.organism.shortName" = "test"), callType = "database"),
               "Unrecognised species passed to function runAllProteinSearch.")
  expect_error(error4 <- runAllProteinSearch(dataToSearch = data.frame("search_term" = "test", "Gene.homologues.homologue.organism.shortName" = "Human"), callType = "1"),
               "callType passed to getOrthologueGenes is neither 'database' or 'taxonomy'")
})

fake_ortho_data1 <- data.frame("Gene.primaryIdentifier" = c("10000", "10000", "10000", "1017", "1017", "1017", "1017", "1017", "1017", "1019"),
                               "Gene.symbol" = c("AKT3", "AKT3", "AKT3", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK4"),
                               "Gene.homologues.homologue.primaryIdentifier" = c("MGI:1345147", "RGD:62390", "ZDB-GENE-110309-3", "FBgn0004107", "MGI:104772", "RGD:70486", "S000000364", "WBGene00019362", "ZDB-GENE-040426-2741", "MGI:88357"),
                               "Gene.homologues.homologue.symbol" = c("Akt3", "Akt3", "Gene symbol not found", "Gene symbol not found", "Cdk2", "Cdk2", "Gene symbol not found", "Gene symbol not found", "Gene symbol not found", "Cdk4"),
                               "Gene.homologues.homologue.organism.shortName" = c("M. musculus", "R. norvegicus", "D. rerio", "D. melanogaster", "M. musculus", "R. norvegicus", "S. cerevisiae", "C. elegans", "D. rerio", "M. musculus"),
                               "Gene.homologues.type" = c("least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue", "least diverged orthologue"))

test_that("Test getProteins input", {
  expect_error(error1 <- getProteins(orthologueData = unlist(fake_ortho_data1), allHumanGenes = fake_ortho_data1[,c("Gene.primaryIdentifier", "Gene.symbol")]),
               "'orthologueData' parameter must be of class data.frame. 'orthologueData' parameter value passed to the 'getProteins' function is class character.")
  expect_error(error1 <- getProteins(orthologueData = fake_ortho_data1[, !colnames(fake_ortho_data1) %in% "Gene.primaryIdentifier"], allHumanGenes = fake_ortho_data1[, c("Gene.primaryIdentifier", "Gene.symbol")]),
               "Missing columns in dataframe passed to getProteins")
  expect_error(error1 <- getProteins(orthologueData = fake_ortho_data1, allHumanGenes = "wrong"),
               "'allHumanGenes' parameter must be of class data.frame. 'allHumanGenes' parameter value passed to the 'getProteins' function is class character.")
  
})
protein_results <- getProteins(orthologueData = fake_ortho_data1, allHumanGenes = fake_ortho_data1[, c("Gene.primaryIdentifier", "Gene.symbol")])

test_that("Test correct run of getProteins", {
  expect_equal(colnames(protein_results), c("Gene.symbol",
                                            "Gene.homologues.homologue.primaryIdentifier",
                                            "Gene.homologues.homologue.organism.shortName",
                                            "Gene.homologues.homologue.symbol",
                                            "search_term",
                                            "protein.accession"))
  expect_equal(nrow(protein_results), 13)
  expect_equal(all(fake_ortho_data1$Gene.homologues.homologue.primaryIdentifier %in% protein_results$Gene.homologues.homologue.primaryIdentifier), TRUE)
  # Cannot test for correct protein outputs as these can change with updates of UniProt.
})

test_that("Test getProteinMatrix input", {
  expect_error(error1 <- getProteinMatrix(proteinResult = 123, orthologueData = fake_ortho_data1), 
               "'proteinResult' parameter must be of class data.frame. 'proteinResult' parameter value passed to the 'getProteinMatrix' function is class numeric.")
  expect_error(error2 <- getProteinMatrix(proteinResult = protein_results[, !colnames(protein_results) %in% c("Gene.symbol", "protein.accession")], orthologueData = fake_ortho_data1),
               "Missing columns in dataframe passed to getProteinMatrix")
  expect_error(error3 <- getProteinMatrix(proteinResult = protein_results, orthologueData = "fake_ortho_data1"),
               "'orthologueData' parameter must be of class data.frame. 'orthologueData' parameter value passed to the 'getProteinMatrix' function is class character.")
  expect_error(error4 <- getProteinMatrix(proteinResult = protein_results, orthologueData = fake_ortho_data1[, !colnames(fake_ortho_data1) %in% c("Gene.symbol", "Gene.primaryIdentifier")]),
               "Missing columns in dataframe passed to getProteinMatrix")
               })

fake_protein_result <- data.frame("Gene.symbol" = c("AKT3", "AKT3", "AKT3", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK4", "AKT3", "CDK2", "CDK4"),
                                  "Gene.homologues.homologue.primaryIdentifier" = c("MGI:1345147", "RGD:62390", "ZDB-GENE-110309-3", "FBgn0004107", "MGI:104772", "RGD:70486", "S000000364", "WBGene00019362", "ZDB-GENE-040426-2741", "MGI:88357", "10000", "1017", "1019"),
                                  "Gene.homologues.homologue.organism.shortName" = c("M. musculus", "R. norvegicus", "D. rerio", "D. melanogaster", "M. musculus", "R. norvegicus", "S. cerevisiae", "C. elegans", "D. rerio", "M. musculus", "Human", "Human", "Human"),
                                  "Gene.homologues.homologue.symbol" = c("Akt3", "Akt3", "Gene symbol not found", "Gene symbol not found", "Cdk2", "Cdk2", "Gene symbol not found", "Gene symbol not found", "Gene symbol not found", "Cdk4", "AKT3", "CDK2", "CDK4"),
                                  "search_term" = c("MGI:1345147", "RGD:62390", "ZDB-GENE-110309-3", "FBgn0004107", "MGI:104772", "RGD:70486", "S000000364", "WBGene00019362", "ZDB-GENE-040426-2741", "MGI:88357", "10000", "1017", "1019"),
                                  "protein.accession" = c("Q9WUA6", "Q63484", "E7EXT6", "P23573", "P97377", "Q63699", "P00546", "O61847", "Q7ZWB1", "P30285", "Q9Y243", "P24941", "P11802"))
                                  

test_that("Test correct run of getProteinMatrix", {
  protein_matrix_result <- getProteinMatrix(proteinResult = fake_protein_result, orthologueData = fake_ortho_data1)
  expect_equal(rownames(protein_matrix_result), c("AKT3", "CDK2", "CDK4"))
  expect_equal(protein_matrix_result$Human, c("Q9Y243", "P24941", "P11802"))
  expect_equal(protein_matrix_result$`R. norvegicus`, c("Q63484", "Q63699", NA))
  expect_equal(protein_matrix_result$Human, c("Q9Y243", "P24941", "P11802"))
  expect_equal(protein_matrix_result$`M. musculus`, c("Q9WUA6", "P97377", "P30285"))
  expect_equal(protein_matrix_result$`D. rerio`, c("E7EXT6", "Q7ZWB1", NA ))
  expect_equal(protein_matrix_result$`C. elegans`, c( NA, "O61847", NA))
  expect_equal(protein_matrix_result$`D. melanogaster`, c(NA, "P23573", NA))
  expect_equal(protein_matrix_result$`S. cerevisiae`, c(NA, "P00546", NA))
})

