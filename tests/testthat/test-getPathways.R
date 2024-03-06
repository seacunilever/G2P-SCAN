library(httptest)

with_mock_dir("httpbin-getPathways", {
  
  test_that("Correct run of getPathways, 2 genes keeping intermediate and terminal", {
    ACHE_BCHE_I_T <- getPathways(c("ACHE", "BCHE"),
                                 c("intermediate", "terminal"))
    expect_equal(colnames(ACHE_BCHE_I_T),
                 c("Derived Gene", "Pathways", "Derived Database", "Pathway Identifier", "Pathway Level", "Pathway Link", "Input Gene Count"))
    expect_equal(all(!"parental" %in% ACHE_BCHE_I_T$`Pathway Level`), TRUE)
    expect_equal(nrow(ACHE_BCHE_I_T), 16)
    expect_equal(unique(ACHE_BCHE_I_T$`Pathway Identifier`),
                 c("R-HSA-1483206", "R-HSA-556833", "R-HSA-112311", "R-HSA-2980736", "R-HSA-1483257", "R-HSA-1483191", "R-HSA-422085", "R-HSA-112315"))
    expect_equal(as.list(ACHE_BCHE_I_T[4, ]), list("Derived Gene" = "ACHE",
                                                   "Pathways" = "Peptide hormone metabolism",
                                                   "Derived Database" = "reactome",
                                                   "Pathway Identifier" = "R-HSA-2980736",
                                                   "Pathway Level"  = "intermediate",
                                                   "Pathway Link" = "https://reactome.org/content/detail/R-HSA-2980736",
                                                   "Input Gene Count" = 2
    ))
    
  })
  
  
  test_that("Check handling when all genes are incorrect", {
    expect_error(getPathways(c("FakeGene", "FakeGene2"),
                             c("parental", "intermediate", "terminal")),
                 "All input genes - FakeGene, FakeGene2 - are NOT found gene symbols.")
  })
  
  test_that("Check Handling when no genes are found in HumanMine", {
    expect_error(getPathways(c("CDk3"),
                             c("parental", "intermediate", "terminal")),
                 "No Human Reactome pathways returned for gene input.")
  })  
  
  test_that("Correct run of getPathways, 1 gene keeping parental, intermediate and terminal", {
    ESR1_P_I_T <- getPathways(c("ESR1"),
                              c("parental","intermediate", "terminal"))
    expect_equal(colnames(ESR1_P_I_T),
                 c("Derived Gene", "Pathways", "Derived Database", "Pathway Identifier", "Pathway Level", "Pathway Link", "Input Gene Count"))
    expect_equal(nrow(ESR1_P_I_T), 34)
    
    expect_equal(unique(ESR1_P_I_T$`Pathway Identifier`),
                 c("R-HSA-2219530", "R-HSA-5688426", "R-HSA-1643685", "R-HSA-5663202", "R-HSA-8939211", "R-HSA-9018519", "R-HSA-9009391", "R-HSA-74160", "R-HSA-212436", "R-HSA-9006925",
                   "R-HSA-392499", "R-HSA-199418", "R-HSA-383280", "R-HSA-1251985", "R-HSA-5689896", "R-HSA-2219528", "R-HSA-6811558", "R-HSA-1257604", "R-HSA-597592", "R-HSA-73857",
                   "R-HSA-8931987", "R-HSA-8939256", "R-HSA-8939902", "R-HSA-3108232", "R-HSA-2990846", "R-HSA-4090294", "R-HSA-162582", "R-HSA-1236394", "R-HSA-9006931", "R-HSA-9006934",
                   "R-HSA-8866910", "R-HSA-8878171", "R-HSA-8878166", "R-HSA-8864260"))
    expect_equal(as.list(ESR1_P_I_T[4, ]), list("Derived Gene" = "ESR1",
                                                "Pathways" = "Diseases of signal transduction by growth factor receptors and second messengers",
                                                "Derived Database" = "reactome",
                                                "Pathway Identifier" = "R-HSA-5663202",
                                                "Pathway Level"  = "intermediate",
                                                "Pathway Link" = "https://reactome.org/content/detail/R-HSA-5663202",
                                                "Input Gene Count" = 1
    ))
    
  })
})
