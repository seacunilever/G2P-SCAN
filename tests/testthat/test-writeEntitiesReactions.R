httptest::with_mock_dir("httpbin-writeEandR", {
  test_that("test update of E&R", {
    createEntitiesReactionsCsv(version = as.integer(1), dataLocation = "Resource_file_test_directory_empty")
    after_dirs <- list.dirs(path = "Resource_file_test_directory_empty", recursive = T)
    expect_equal(after_dirs, c("Resource_file_test_directory_empty", "Resource_file_test_directory_empty/g2pScanData", "Resource_file_test_directory_empty/g2pScanData/Reactome_entities_reactions"))
    file_output <- file.path("Resource_file_test_directory_empty", "g2pScanData", "Reactome_entities_reactions", "Reactome_entities_reactions_v1.csv")
    expect_equal(file.exists(file_output), TRUE)
    actual_output <- read.csv(file_output)
    for(i in after_dirs){
      R.utils::removeDirectory(i, recursive = T, mustExist = F)
    }
    expect_equal(dim(actual_output), c(2546, 15))
    expect_equal(as.list(actual_output[1, ]), list("pathway" = "R-HSA-1059683",
                                                   "Human_entities" = 17,
                                                   "R..norvegicus_entities" = 8,
                                                   "M..musculus_entities" = 8,
                                                   "D..rerio_entities" = 6,
                                                   "C..elegans_entities" = 0 ,
                                                   "D..melanogaster_entities" = 0,
                                                   "S..cerevisiae_entities" = 0,
                                                   "Human_reactions" = 20,
                                                   "R..norvegicus_reactions" = 18,
                                                   "M..musculus_reactions" = 18,
                                                   "D..rerio_reactions" = 18,
                                                   "C..elegans_reactions" = 0,
                                                   "D..melanogaster_reactions" = 0,
                                                   "S..cerevisiae_reactions" = 0))
    
    expect_equal(as.list(actual_output[100, ]), list("pathway" = "R-HSA-1237112",
                                                     "Human_entities" =  24,
                                                     "R..norvegicus_entities" = 6,
                                                     "M..musculus_entities" = 6,
                                                     "D..rerio_entities" = 6,
                                                     "C..elegans_entities" = 5 ,
                                                     "D..melanogaster_entities" = 5,
                                                     "S..cerevisiae_entities" = 5,
                                                     "Human_reactions" = 7,
                                                     "R..norvegicus_reactions" = 7,
                                                     "M..musculus_reactions" = 7,
                                                     "D..rerio_reactions" = 7,
                                                     "C..elegans_reactions" = 6,
                                                     "D..melanogaster_reactions" = 6,
                                                     "S..cerevisiae_reactions" = 6
    ))
    
    expect_equal(as.list(actual_output[1000, ]), list("pathway" = "R-HSA-392499",
                                                      "Human_entities" =  2205,
                                                      "R..norvegicus_entities" = 1562,
                                                      "M..musculus_entities" = 1616,
                                                      "D..rerio_entities" = 1376,
                                                      "C..elegans_entities" = 769 ,
                                                      "D..melanogaster_entities" = 870,
                                                      "S..cerevisiae_entities" = 431,
                                                      "Human_reactions" = 798,
                                                      "R..norvegicus_reactions" = 687,
                                                      "M..musculus_reactions" = 695,
                                                      "D..rerio_reactions" = 672,
                                                      "C..elegans_reactions" = 525,
                                                      "D..melanogaster_reactions" = 557,
                                                      "S..cerevisiae_reactions" = 419))
    
  })
})