#' Internal mypackage functions
#'
#' Internal mypackage functions
#' @details
#' These functions are not intended to be called by the user.
#' @name mypackage-internal
#' @keywords internal
NULL

# ================= Manipulation of formulae for setThreshold =============== #

# Function to replace all formula terms of the form cy(x) or pbc(x) with
# cy(x, cyclic = TRUE) and pbc(cyclic = TRUE)
# The input "formula" is a formula
# The input which determine which of cy() or pbc() to deal with
# The output "formula" is a deparsed formula
cy_and_pbc <- function(formula, which = c("cy", "pbc")) {
  which <- match.arg(which)
  # Extract the terms object form the formula
  terms_obj <- terms(formula)
  # The attribute "term.labels" gives the labels (names) of the terms
  term_labels <- attr(terms(formula), "term.labels")
  # Turn the formula into a character string, so that we manipulate it
  formula <- deparse(formula)
  # Pattern to replace
  pattern <- paste0(which, "(")
  # Loop over all the terms and add cyclic = TRUE were required
  for (label in term_labels) {
    condition  <- any(grepl(pattern, label, fixed = TRUE))
    if (condition) {
      replace_with <- add_cyclic_equals_true(label)
      formula <- gsub(label, replace_with, formula, fixed = TRUE)
    }
  }
  return(formula)
}

# Function to convert cy(x) to cy(x, cyclic = TRUE), for example
add_cyclic_equals_true <- function(x) {
  # Split the formula into bits
  bits <- unlist(strsplit(x, "[[:<:]]", perl= TRUE))
  # Find the bits between "(" and ")"
  brackets <- which(bits == "(" | bits == ")")
  between_brackets <- (brackets[1] + 1):(brackets[2] - 1)
  # Create the new argument between ()
  # "x" becomes "x, cyclic = TRUE", for example
  variable_name <- paste0(bits[between_brackets], collapse = "")
  variable_name <- paste0(variable_name, ", cyclic = TRUE")
  # Create the text around "x, cyclic = TRUE"
  start <- paste0(bits[1:brackets[1]], collapse = "")
  end <- bits[brackets[2]]
  # Create the term new
  new_term <- paste0(start, variable_name, end)
  return(new_term)
}
