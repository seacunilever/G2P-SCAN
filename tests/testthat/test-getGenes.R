httptest::with_mock_dir("httpbin-getGenes", {
 test_that("Correct run of getGenes, with one pathway", {
    genesInPath <- getGenes("R-HSA-112311")
    expect_equal(colnames(genesInPath), c("Gene Identifier", "Gene Symbol", "Pathway Identifier"))
    expect_equal(nrow(genesInPath), 10)
    expect_equal(genesInPath$`Gene Identifier`, c("1312", "217", "220074", "4128", "43", "590", "6531", "6532", "6580", "6582"))
    expect_equal(genesInPath$`Gene Symbol`, c("COMT", "ALDH2", "LRTOMT", "MAOA", "ACHE", "BCHE", "SLC6A3", "SLC6A4", "SLC22A1", "SLC22A2"))
    expect_equal(unique(genesInPath$`Pathway Identifier`), "R-HSA-112311")
    
})
  
  test_that("Correct run of getGenes, with two pathway", {
    genesInPath <- getGenes(c("R-HSA-112311", "R-HSA-1430728"))
    expect_equal(colnames(genesInPath), c("Gene Identifier", "Gene Symbol", "Pathway Identifier"))
    expect_equal(nrow(genesInPath), 2073)
    expect_equal(genesInPath$`Gene Identifier`[c(10,100,1000,2000)], c("6582", "10894", "4719", "94235"))
    expect_equal(genesInPath$`Gene Symbol`[c(10,100,1000,2000)], c("SLC22A2", "LYVE1", "NDUFS1", "GNG8"))
    expect_equal(unique(genesInPath$`Pathway Identifier`), c("R-HSA-112311", "R-HSA-1430728"))
    
  })  
})
