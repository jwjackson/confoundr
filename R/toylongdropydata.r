#' Artifical data set used to test the functionality of confoundr.
#'
#' The toy_long_dropoutY data set contains 28,410 records and 16
#' variables. These variables include time-varying exposures,
#' outcomes, and covariates, along with strata and censoring
#' indicators. Toy data are removed after s equals one.
#' Time-varying inverse-probability-of-exposure
#' weights and censoring weights are available as well.
#'
#' @docType data
#' @usage data(toy_long_dropoutY)
#'
#' @format A data frame with 2,847 rows and 13 variables:
#' \describe{
#'     \item{uid}{subject ID}
#'     \item{time}{time of observation}
#'     \item{a}{exposure measurement at time t}
#'     \item{l}{covariate measurement at time t}
#'     \item{m}{covariate measurement at time t}
#'     \item{n}{covariate measurement at time t}
#'     \item{o}{covariate measurement at time t}
#'     \item{p}{covariate measurement at time t}
#'     \item{s}{censoring indicator at time t}
#'     \item{h}{exposure history at time t}
#'     \item{hx}{grouped exposure history, by p_0, at time t}
#'     \item{wa}{inverse probability of exposure and censoring weight at time t}
#'     \item{wax}{cumulative inverse probability of exposure weight at time t}
#'     \item{wsx}{cumulative inverse probability of censoring weight at time t}
#'     \item{e5}{propensity score strata at time t}
#'
#'     }
#'
"toy_long_dropoutY"
