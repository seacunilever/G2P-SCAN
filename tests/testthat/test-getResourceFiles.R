# test getReactomeFileOfVersionPath

test_that("Test when directories/files do not exist", {
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "reactome_hierarchy", folder = "pathway_hierarchy_datatrees", version = 76, dataLocation = "Resource_file_test_directory_empty"),
    NULL)

  after_dirs <- list.dirs(path = "Resource_file_test_directory_empty", recursive = T)
  for(i in after_dirs){
    R.utils::removeDirectory(i, recursive = T, mustExist = F)
  }
  
  expect_equal(after_dirs, c("Resource_file_test_directory_empty", "Resource_file_test_directory_empty/g2pScanData", "Resource_file_test_directory_empty/g2pScanData/pathway_hierarchy_datatrees"))
})

test_that("Test if correct file is being returned", {
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "reactome_hierarchy", folder = "pathway_hierarchy_datatrees", version = 76, dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/pathway_hierarchy_datatrees/test_v76.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "reactome_hierarchy", folder = "pathway_hierarchy_datatrees", version = 77, dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/pathway_hierarchy_datatrees/test_v77.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "reactome_hierarchy", folder = "pathway_hierarchy_datatrees", version = "latest", dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/pathway_hierarchy_datatrees/test_v78.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "reactome_hierarchy", folder = "pathway_hierarchy_datatrees", version = 79, dataLocation = "Resource_file_test_directory"),
    NULL)
  
})

test_that("Test when directories/files do not exist", {
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file", folder = "Reactome_entities_reactions", version = 76, dataLocation = "Resource_file_test_directory_empty"),
    NULL)

  after_dirs <- list.dirs(path = "Resource_file_test_directory_empty", recursive = T)
  for(i in after_dirs){
    R.utils::removeDirectory(i, recursive = T, mustExist = F)
  }
  
  expect_equal(after_dirs, c("Resource_file_test_directory_empty", "Resource_file_test_directory_empty/g2pScanData", "Resource_file_test_directory_empty/g2pScanData/Reactome_entities_reactions"))
})

test_that("Test if correct file is being returned", {
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file", folder = "Reactome_entities_reactions", version = 76, dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/Reactome_entities_reactions/test_v76.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file", folder = "Reactome_entities_reactions", version = 77, dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/Reactome_entities_reactions/test_v77.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file", folder = "Reactome_entities_reactions", version = "latest", dataLocation = "Resource_file_test_directory"),
    "Resource_file_test_directory/g2pScanData/Reactome_entities_reactions/test_v78.Rdata")
  expect_equal(
    getReactomeFileOfVersionPath(dataName = "Reactome entities and reactions file", folder = "Reactome_entities_reactions", version = 79, dataLocation = "Resource_file_test_directory"),
    NULL)
  
})

