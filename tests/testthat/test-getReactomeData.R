reactome_v <- R.utils::getOption(localVars, "reactome_version")

test_that("Test error when all pathways input is wrong", {
  expect_error(getEntitiesReactionsCounts(pathways = c("wrong1", "wrong2", "wrong3"), species = "R. norvegicus"),
               paste0("No variables given as pathways are identified as a pathway in Reactome version ", reactome_v, "."))
})

test_that("Test error when all species input is wrong", {
  expect_error(getEntitiesReactionsCounts(pathways = NULL, species = c("wrong1", "wrong2", "wrong3")),
               "Unrecognised species passed to function filterEntitiesReactionsSpecies.")
})

#If this were to fail in test it would be caused by function failing OR names pathway no longer existed in updated version of Reactome.
test_that("Check Filtering correct", {
  output <- getEntitiesReactionsCounts(pathways = c("R-HSA-391903", "R-HSA-391906", "R-HSA-391908", "R-HSA-392023", "R-HSA-392154", "R-HSA-392170", "R-HSA-392451", "R-HSA-392499"),
                                      species = c("R. norvegicus", "M. musculus"))
  expect_equal(ncol(output), 6)
  expect_equal(nrow(output), 8)
})

