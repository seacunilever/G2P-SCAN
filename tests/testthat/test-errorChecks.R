## Test checkClass - class types tested are those which the function is used for in the package.
test_that("Check handling of list when character vector is expected", {
  expect_error(checkClass(parameter = "test_run_parameter", value = list("R-HSA-112311"), classType = "character", functionName = "test_run_function"),
               "'test_run_parameter' parameter must be of class character. 'test_run_parameter' parameter value passed to the 'test_run_function' function is class list.")
})
test_that("Check handling of number when character vector is expected", {
  expect_error(checkClass(parameter = "test_run_parameter", value = 2, classType = "character", functionName = "test_run_function"),
               "'test_run_parameter' parameter must be of class character. 'test_run_parameter' parameter value passed to the 'test_run_function' function is class numeric.")
})
test_that("Check handling of dataframe when character vector is expected", {
  expect_error(checkClass(parameter = "test_run_parameter", value = data.frame("test" = 1), classType = "character", functionName = "test_run_function"),
               "'test_run_parameter' parameter must be of class character. 'test_run_parameter' parameter value passed to the 'test_run_function' function is class data.frame.")
})
test_that("Check handling of character when character vector is expected", {
  expect_error(checkClass(parameter = "test_run_parameter", value = "correct", classType = "character", functionName = "test_run_function"),
               NA)
})

## Test checkSpecies 
test_that("Check output of checkSpecies", {
  ##correct inputs
  correct_species <- checkSpecies(species = c("R. norvegicus", "M. musculus"), functionName = "test_run_function", includeHuman = FALSE)
  expect_equal(correct_species, c("R. norvegicus", "M. musculus"))
  correct_species2 <- checkSpecies(species =  c("R. norvegicus", "M. musculus", "D. rerio", "C. elegans"), functionName = "test_run_function", includeHuman = FALSE)
  expect_equal(correct_species2, c("R. norvegicus", "M. musculus", "D. rerio", "C. elegans"))
  correct_species3 <- checkSpecies(species =  c("R. norvegicus", "M. musculus", "D. rerio", "C. elegans"), functionName = "test_run_function", includeHuman = TRUE)
  expect_equal(correct_species3, c("R. norvegicus", "M. musculus", "D. rerio", "C. elegans"))
  correct_species4 <- checkSpecies(species = c("R. norvegicus", "M. musculus", "Human"), functionName = "test_run_function", includeHuman = TRUE)
  expect_equal(correct_species4, c("Human", "R. norvegicus", "M. musculus"))
  correct_species5 <- checkSpecies(species = NULL, functionName = "test_run_function", includeHuman = TRUE)
  expect_equal(correct_species5, c( "Human", "R. norvegicus", "M. musculus", "D. rerio", "C. elegans", "D. melanogaster", "S. cerevisiae"))
})
test_that("Check handling when species input is incorrect", {
  ##correct inputs
  expect_error(checkSpecies(species = c("R. norvegicus", "M. musculus", "Human"), functionName = "test_run_function", includeHuman = FALSE),
               "Unrecognised species passed to function test_run_function.")
  expect_error(checkSpecies(species = c("R. norvegicus", "M. musculus", "wrongSpecies"), functionName = "test_run_function", includeHuman = FALSE),
               "Unrecognised species passed to function test_run_function.")

})

## Test checkColums
test_that("Check output of checkSpecies", {
  ##correct inputs
  dt_test <- data.frame("col1" = c(1,1,1), "col2" = c(1,1,1), "col3" = c(1,1,1), "col4" = c(1,1,1), "col5" = c(1,1,1), "col6" = c(1,1,1), "col7" = c(1,1,1), "col8" = c(1,1,1))
  expect_error(checkColumns(dt= dt_test, columns = c("col1", "col2", "col7"), dtParameter = "test_parameter", functionName = "test_function"), 
               NA)
  expect_error(checkColumns(dt= dt_test, columns = c("col1", "col2", "WRONG"), dtParameter = "test_parameter", functionName = "test_function"), 
               "Missing columns in dataframe passed to test_function")
  })
