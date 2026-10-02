# Same values and same R-level layout: identical(), attribute order, and the
# compact/expanded row.names form. Raw serialize() bytes are deliberately not
# compared: they also encode internal string flags that depend on the R
# version and on how a vector was built (e.g. rbind vs c()), not on its value.
expect_same_representation <- function(object, expected) {
  testthat::expect_identical(object, expected)
  testthat::expect_identical(names(attributes(object)), names(attributes(expected)))
  if (is.data.frame(object) && is.data.frame(expected)) {
    testthat::expect_identical(.row_names_info(object, 0L), .row_names_info(expected, 0L))
    testthat::expect_identical(lapply(object, function(x) names(attributes(x))),
                               lapply(expected, function(x) names(attributes(x))))
  }
}
