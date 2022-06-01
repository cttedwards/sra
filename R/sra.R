#' @title Extract and compile sra model 
#'
#' @description The lastest version of the seabird risk assessment model is extracted and compiled
#' 
#' @importFrom utils packageVersion
#' @importFrom rstan stan_model
#' @export
sra <- function() {
  stan_model(file = system.file("extdata/stan", "sra_v0.0.2.stan", package = "sra"), model_name = paste0("sra v", packageVersion("sra")))
}
