#' @title Extract and compile sra model 
#'
#' @importFrom rstan stan_model
#' @export
sra <- function() {
  stan_model(file = system.file("extdata/stan", "sra.stan", package = "sra"), model_name = paste0("sra v", packageVersion("sra")))
}
