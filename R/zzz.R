.onAttach <- function(libname, pkgname) {

  # Load lavaan
  # cli::cli_alert_info(
  #   "Loading required package: {.pkg lavaan}",
  #   class = "packageStartupMessage"
  # )
  # library("lavaan", exclude = "cfa")

  # Inform conflicts
  # conflict_msg <- paste0(
  #   cli::rule(left = cli::style_bold("Conflicts")),
  #   "\n",
  #   bullets_mask("plFA", "lavaan", "cfa")
  # )
  # rlang::inform(conflict_msg, class = "packageStartupMessage")
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
