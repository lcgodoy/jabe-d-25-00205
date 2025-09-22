create_debug <- function(file) {
  aux <- readLines(file)
  lines_to_modify <-
    aux |>
    trimws() |>
    grep(pattern = "^target \\+=", x = _)
  for (k in lines_to_modify)
    aux[k] <- paste(aux[k], 'print("target = ", target());')
  writeLines(text = aux, con = gsub("\\.stan$", "-debug\\.stan", file))
}
