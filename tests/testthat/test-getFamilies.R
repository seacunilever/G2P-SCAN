fake_protein_result <- data.frame("Gene.symbol" = c("AKT3", "AKT3", "AKT3", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK2", "CDK4", "AKT3", "CDK2", "CDK4"),
                                  "Gene.homologues.homologue.primaryIdentifier" = c("MGI:1345147", "RGD:62390", "ZDB-GENE-110309-3", "FBgn0004107", "MGI:104772", "RGD:70486", "S000000364", "WBGene00019362", "ZDB-GENE-040426-2741", "MGI:88357", "10000", "1017", "1019"),
                                  "Gene.homologues.homologue.organism.shortName" = c("M. musculus", "R. norvegicus", "D. rerio", "D. melanogaster", "M. musculus", "R. norvegicus", "S. cerevisiae", "C. elegans", "D. rerio", "M. musculus", "Human", "Human", "Human"),
                                  "Gene.homologues.homologue.symbol" = c("Akt3", "Akt3", "Gene symbol not found", "Gene symbol not found", "Cdk2", "Cdk2", "Gene symbol not found", "Gene symbol not found", "Gene symbol not found", "Cdk4", "AKT3", "CDK2", "CDK4"),
                                  "search_term" = c("MGI:1345147", "RGD:62390", "ZDB-GENE-110309-3", "FBgn0004107", "MGI:104772", "RGD:70486", "S000000364", "WBGene00019362", "ZDB-GENE-040426-2741", "MGI:88357", "10000", "1017", "1019"),
                                  "protein.accession" = c("Q9WUA6", "Q63484", "E7EXT6", "P23573", "P97377", "Q63699", "P00546", "O61847", "Q7ZWB1", "P30285", "Q9Y243", "P24941", "P11802"))


test_that("Test getFamilies input", {
  expect_error(error1 <- getFamilies(proteinResult = 123),
               "'proteinResult' parameter must be of class data.frame. 'proteinResult' parameter value passed to the 'getFamilies' function is class numeric.")
  expect_error(error2 <- getFamilies(proteinResult = fake_protein_result[, !colnames(fake_protein_result) %in% c("Gene.symbol", "protein.accession")]),
               "Missing columns in dataframe passed to getFamilies")
})

fam_result <- getFamilies(proteinResult = fake_protein_result)

test_that("Test correct run of getFamilies", {
  expect_equal(colnames(fam_result), c("protein.accession",
                                       "Gene.homologues.homologue.organism.shortName",
                                       "filtered_family"))
  expect_equal(nrow(fam_result), 13)
  expect_equal(all(fake_protein_result$protein.accession %in% fam_result$protein.accession), TRUE)
  # Cannot test for correct protein outputs as these can change with updates of UniProt.
})