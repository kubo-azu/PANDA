# Column names evaluated inside dplyr data masks. Declaring them here keeps
# R CMD check's static analysis distinct from ordinary R object lookup.
utils::globalVariables(c(
  "Count", "Meth", "Methylation", "Num_Meth", "Num_Total", "Pattern",
  "Pct", "ReadID", "Total"
))
