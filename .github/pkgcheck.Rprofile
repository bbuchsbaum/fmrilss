options(
  repos = c(
    CRAN = "https://cloud.r-project.org",
    bbuchsbaum = "https://bbuchsbaum.r-universe.dev"
  )
)

# pkgcheck recognizes HTML vignettes by output-format name. albers_vignette
# produces HTML but its name does not contain "html". Preserve upstream checks
# and inspect the actual Pandoc output type for this one supported format.
setHook(packageEvent("pkgcheck", "onLoad"), function(...) {
  ns <- asNamespace("pkgcheck")
  original <- get("pkgchk_has_vignette", envir = ns)
  compatible <- function(checks) {
    if (isTRUE(original(checks))) return(TRUE)
    paths <- list.files(file.path(checks$pkg$path, "vignettes"),
                        pattern = "\\.[rR]md$", recursive = TRUE, full.names = TRUE)
    any(vapply(paths, function(path) {
      tryCatch({
        output <- rmarkdown::yaml_front_matter(path)$output
        name <- "albersdown::albers_vignette"
        if (is.character(output)) {
          if (!identical(output, name)) return(FALSE)
          args <- list()
        } else {
          if (!name %in% names(output)) return(FALSE)
          args <- output[[name]]
          if (is.null(args) || identical(args, "default")) args <- list()
        }
        format <- do.call(albersdown::albers_vignette, args)
        inherits(format, "rmarkdown_output_format") &&
          format$pandoc$to %in% c("html", "html4", "html5")
      }, error = function(e) FALSE)
    }, logical(1)))
  }
  assignInNamespace("pkgchk_has_vignette", compatible, ns = "pkgcheck")
})
