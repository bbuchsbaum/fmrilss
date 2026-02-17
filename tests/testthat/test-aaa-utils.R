# Tests for internal utilities in aaa_utils.R
# These are internal functions but we test them through their usage

test_that(".as_base_matrix handles NULL", {
  expect_null(fmrilss:::.as_base_matrix(NULL))
})

test_that(".as_base_matrix handles base matrix", {
  m <- matrix(1:6, 2, 3)
  result <- fmrilss:::.as_base_matrix(m)
  expect_equal(result, m)
  expect_true(is.matrix(result))
})

test_that(".as_base_matrix handles data.frame", {
  df <- data.frame(a = 1:3, b = 4:6)
  result <- fmrilss:::.as_base_matrix(df)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 2))
})

test_that(".as_base_matrix handles Matrix objects", {
  skip_if_not_installed("Matrix")
  m <- Matrix::Matrix(1:6, 2, 3)
  result <- fmrilss:::.as_base_matrix(m)
  expect_true(is.matrix(result))
  expect_false(inherits(result, "Matrix"))
})

test_that("%||% returns first value if not NULL", {
  `%||%` <- fmrilss:::`%||%`
  expect_equal(5 %||% 10, 5)
  expect_equal("a" %||% "b", "a")
  expect_equal(0 %||% 10, 0)  # 0 is not NULL
  expect_equal(FALSE %||% TRUE, FALSE)  # FALSE is not NULL
})

test_that("%||% returns second value if first is NULL", {
  `%||%` <- fmrilss:::`%||%`
  expect_equal(NULL %||% 10, 10)
  expect_equal(NULL %||% "default", "default")
})

test_that(".default_trial_names generates correct names", {
  names1 <- fmrilss:::.default_trial_names(1)
  expect_equal(names1, "Trial_1")

  names5 <- fmrilss:::.default_trial_names(5)
  expect_equal(names5, c("Trial_1", "Trial_2", "Trial_3", "Trial_4", "Trial_5"))
})

test_that(".trial_names_from_cols returns colnames if present", {
  X <- matrix(1:12, 4, 3)
  colnames(X) <- c("A", "B", "C")
  result <- fmrilss:::.trial_names_from_cols(X)
  expect_equal(result, c("A", "B", "C"))
})

test_that(".trial_names_from_cols falls back to default names", {
  X <- matrix(1:12, 4, 3)
  result <- fmrilss:::.trial_names_from_cols(X)
  expect_equal(result, c("Trial_1", "Trial_2", "Trial_3"))
})

test_that(".trial_names_from_cols uses n parameter", {
  result <- fmrilss:::.trial_names_from_cols(NULL, n = 4)
  expect_equal(result, c("Trial_1", "Trial_2", "Trial_3", "Trial_4"))
})

test_that(".set_beta_dimnames sets row names", {
  beta <- matrix(1:6, 2, 3)
  result <- fmrilss:::.set_beta_dimnames(beta, trial_names = c("T1", "T2"))
  expect_equal(rownames(result), c("T1", "T2"))
})

test_that(".set_beta_dimnames sets column names", {
  beta <- matrix(1:6, 2, 3)
  result <- fmrilss:::.set_beta_dimnames(beta, voxel_names = c("V1", "V2", "V3"))
  expect_equal(colnames(result), c("V1", "V2", "V3"))
})

test_that(".set_beta_dimnames generates default names when NULL", {
  beta <- matrix(1:6, 2, 3)
  result <- fmrilss:::.set_beta_dimnames(beta)
  expect_equal(rownames(result), c("Trial_1", "Trial_2"))
})

test_that(".set_beta_dimnames preserves existing names", {
  beta <- matrix(1:6, 2, 3)
  rownames(beta) <- c("Existing1", "Existing2")
  result <- fmrilss:::.set_beta_dimnames(beta, trial_names = c("New1", "New2"))
  # Should preserve existing since rownames already set
  expect_equal(rownames(result), c("Existing1", "Existing2"))
})
