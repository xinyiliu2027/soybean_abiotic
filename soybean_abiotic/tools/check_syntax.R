# Run from the repository root. Parses scripts without loading analysis packages.
files <- unlist(lapply(c("scripts", "examples", "legacy"), function(d) {
  list.files(d, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
}))
failed <- character()
for (f in files) {
  tryCatch({parse(file = f, encoding = "UTF-8"); cat("PASS", f, "\n")},
           error = function(e) {failed <<- c(failed, f); cat("FAIL", f, conditionMessage(e), "\n")})
}
if (length(failed)) quit(status = 1L)
