if (!exists("skip_if_installed")) {
  skip_if_installed <- function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      testthat::skip(sprintf("%s installed", package))
    }
    invisible(TRUE)
  }
}

