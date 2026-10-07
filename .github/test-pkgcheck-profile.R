# Run in a fresh R process: checks both upstream behavior and our narrow hook.
source(".github/pkgcheck.Rprofile")
invisible(loadNamespace("pkgcheck"))
check <- get("pkgchk_has_vignette", asNamespace("pkgcheck"))
root <- tempfile("pkgcheck-format-")
dir.create(file.path(root, "vignettes"), recursive = TRUE)
fixture <- function(format) {
  writeLines(c("---", 'title: "Format check"', paste0("output: ", format),
               "---", "Body."), file.path(root, "vignettes", "test.Rmd"))
  check(list(pkg = list(path = root)))
}
stopifnot(fixture("html_document"), fixture("albersdown::albers_vignette"),
          !fixture("pdf_document"), !fixture("word_document"))
unlink(root, recursive = TRUE)
cat("pkgcheck format compatibility: HTML accepted; PDF/Word rejected\n")
