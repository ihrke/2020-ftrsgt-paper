map(list.files(
  path = "lib",
  pattern = "*\\.cpp",
  full.names = TRUE
),
function(.x) {
  message(sprintf("\t\t\tRcpp::sourceCpp('%s')",.x))
  Rcpp::sourceCpp(.x)
})