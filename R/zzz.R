.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("palimpsestr")
  packageStartupMessage(sprintf("palimpsestr %s — Probabilistic Palimpsest Analysis", v))
}
