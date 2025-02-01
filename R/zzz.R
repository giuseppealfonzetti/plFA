.onAttach <- function(libname, pkgname) {

  # Load lavaan
  startup_msg <- paste0(
    cli::cli_fmt(cli::cli_alert_info("Loading required package: {.pkg lavaan}"))
  )
  rlang::inform(startup_msg, class = "packageStartupMessage")
  attach_lavaan(without = "cfa")

  # Inform conflicts
  conflict_msg <- paste0(
    cli::rule(
      left = cli::style_bold("Conflicts"),
      right = paste("plFA", utils::packageVersion("plFA"))
    ),
    "\n",
    bullets_mask("plFA", "lavaan", "cfa")
  )
  rlang::inform(conflict_msg, class = "packageStartupMessage")
}

attach_lavaan <- function(without) {
  do.call("library", list(
    package = "lavaan",
    exclude = without,
    character.only = TRUE
  ))
}

format_pkgfun <- function(pkg, fun) {
  paste0(cli::col_blue(pkg), "::", cli::col_green(fun, "()"))
}

bullets_mask <- function(winner, loser, fun) {
  paste0(
    cli::col_red(cli::symbol$cross), " ",
    format_pkgfun(winner, fun),
    " masks ",
    format_pkgfun(loser, fun)
  )
}
