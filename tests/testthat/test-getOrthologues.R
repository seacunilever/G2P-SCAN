httptest::with_mock_dir("httpbin-getOrthologues", {
  # tests for getOrthologuesGenes
  test_that("Incorrect gene input", {
    expect_error(getOrthologueGenes(genes = c(1,2,3), orthologueFilter = "ldo"),
                 "'genes' parameter must be of class character. 'genes' parameter value passed to the 'getOrthologueGenes' function is class numeric.")
    
  })
  test_that("Incorrect orthologue input", {
    expect_error(getOrthologueGenes(genes = "BCHE", orthologueFilter = "1"),
                 "OrthologueFilter passed to getOrthologueGenes is neither 'LDO' or 'All'")
    
    expect_error(getOrthologueGenes(genes = "BCHE", orthologueFilter = 1),
                 "'orthologueFilter' parameter must be of class character. 'orthologueFilter' parameter value passed to the 'getOrthologueGenes' function is class numeric.")
  })
  
  test_that("Incorrect species input", {
    expect_error(getOrthologueGenes(genes = "BCHE", orthologueFilter = "ldo", species = "fake"),
                 "Unrecognised species passed to function filterOrthologueSpecies")
    
    expect_error(getOrthologueGenes(genes = "BCHE", orthologueFilter = "ldo", species = c("fake", "R. norvegicus")),
                 "Unrecognised species passed to function filterOrthologueSpecies")
  })
  
  test_that("Test when no results are returned", {
    expect_error(getOrthologueGenes("afkadkf", "ldo", species = NULL),
                 "No orthologue results returned.") # no orthologue results at all (before filtering)
    expect_error(getOrthologueGenes(genes = "BCHE", orthologueFilter = "ldo", species = "S. cerevisiae"),
                 "No orthologue results returned.") # no orthologue results after species filter
    expect_error(getOrthologueGenes(genes = "ZNF683", orthologueFilter = "ldo", species = "S. cerevisiae"),
                 "No orthologue results returned.") # no orthologue results after species and orthologue filter
    
  })
  
  test_that("Test when species = NULL, correct run", {
    null_species <- getOrthologueGenes(genes = "BCHE", orthologueFilter = "ldo")
    expect_equal(dim(null_species), c(4, 6))
    expect_equal(colnames(null_species), c("Gene.primaryIdentifier",
                                           "Gene.symbol",
                                           "Gene.homologues.homologue.primaryIdentifier",
                                           "Gene.homologues.homologue.symbol",
                                           "Gene.homologues.homologue.organism.shortName",
                                           "Gene.homologues.type"))
    expect_equal(unique(null_species$Gene.symbol), "BCHE")
    expect_equal(null_species$Gene.homologues.homologue.symbol, c("Gene symbol not found",
                                                                  "Bche",
                                                                  "Bche",
                                                                  "Gene symbol not found"
    ))
    expect_equal(null_species$Gene.homologues.homologue.organism.shortName, c("D. melanogaster",
                                                                              "M. musculus",
                                                                              "R. norvegicus",
                                                                              "C. elegans"
    ))
  })
  test_that("Test correct species filter", {
    mouse <- getOrthologueGenes(genes = "BCHE", orthologueFilter = "ldo", species = "M. musculus")
    expect_equal(dim(mouse), c(1, 6))
    expect_equal(unique(mouse$Gene.homologues.homologue.organism.shortName), "M. musculus")
  })
  
  test_that("Test correct orthologue filter", {
    ldos <- getOrthologueGenes(genes = c("BCHE", "AADACL4"), orthologueFilter = "ldo", species = "M. musculus")
    all <- getOrthologueGenes(genes = c("BCHE", "AADACL4"), orthologueFilter = "ALL", species = "M. musculus")
    expect_equal(dim(ldos), c(2, 6))
    expect_equal(dim(all), c(7, 6))
    expect_equal(unique(ldos$Gene.homologues.type), "least diverged orthologue")
    expect_equal(unique(all$Gene.homologues.type), c("least diverged orthologue", "orthologue"))
  })
  
  # tests for getOrthologueMatrix
  ortho_data_all <- getOrthologueGenes(genes = c("ACHE", "BCHE", "AADACL4"), orthologueFilter = "all", species = NULL)
  
  test_that("Test incorrect orthologueData input", {
    expect_error(getOrthologueMatrix(orthologueData = data.frame(c(1,2,3), c("a","b", "c")), allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "both"),
                 "Missing columns in dataframe passed to getOrthologueMatrix")
    expect_error(getOrthologueMatrix(orthologueData = ortho_data_all[, 1:5], allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "both"),
                 "Missing columns in dataframe passed to getOrthologueMatrix")
    
    expect_error(getOrthologueMatrix(orthologueData = "not_data_frame", allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "both"),
                 "orthologueData' parameter must be of class data.frame. 'orthologueData' parameter value passed to the 'getOrthologueMatrix' function is class character.")
  })
  test_that("Test incorrect allHumanGenes input", {
    expect_error(getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = 123),
                 "'allHumanGenes' parameter must be of class character. 'allHumanGenes' parameter value passed to the 'getOrthologueMatrix' function is class numeric.")
  })
  test_that("Test incorrect orthologueOutput input", {
    expect_error(getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = 123),
                 "'orthologueOutput' parameter must be of class character. 'orthologueOutput' parameter value passed to the 'getFormatOrthologue' function is class numeric.")
    expect_error(getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "wrong"),
                 "orthologueOutput passed to getFormatOrthologue is neither 'primaryIdentifier', 'geneSymbol' or 'both'")
  })
  
  ortho_matrix_both <- getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "both")
  ortho_matrix_priID <- getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "primaryIdentifier")
  ortho_matrix_symbol <- getOrthologueMatrix(orthologueData = ortho_data_all, allHumanGenes = c("ACHE", "BCHE", "AADACL4"), orthologueOutput = "genesymbol")
  
  test_that("Test species order is maintained", {
    all_species_name_code <- R.utils::getOption(localVars, "species_list")
    all_species_name <- names(all_species_name_code)
    expect_equal(colnames(ortho_matrix_both), all_species_name)
  })
  
  test_that("Test correct run and test they use correct orthologue output", {
    
    expect_equal(dim(ortho_matrix_both), c(3,6))
    expect_equal(dim(ortho_matrix_priID), c(3,6))
    expect_equal(dim(ortho_matrix_symbol), c(3,6))
    
    expect_equal(ortho_matrix_both[1,1],
                 c("RGD:1308878(RGD1308878), RGD:1559644(RGD1559644), RGD:1563334(RGD1563334), RGD:1565761(Aadacl4), RGD:1587396(LOC691196)"))
    expect_equal(ortho_matrix_priID[1,1],
                 c("RGD:1308878, RGD:1559644, RGD:1563334, RGD:1565761, RGD:1587396"))
    expect_equal(ortho_matrix_symbol[1,1],
                 c("RGD1308878, RGD1559644, RGD1563334, Aadacl4, LOC691196"))
  })
  
})





