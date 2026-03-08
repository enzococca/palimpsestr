.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("palimpsestr")
  packageStartupMessage(sprintf("palimpsestr %s \u2014 Probabilistic Palimpsest Analysis", v))
}
