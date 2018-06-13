#' Artifical data set used to illustrate the functionality of confoundr.
#'
#' The example_sml data set contains 10,000 records and 38
#' variables. These variables include time-varying exposures,
#' outcomes, and covariates, along with strata and censoring
#' indicators. Time-varying inverse-probability-of-exposure
#' weights and censoring weights are available as well.
#'
#' @docType data
#' @usage data(example_sml)
#'
#' @format A data frame with 10,000 rows and 38 variables:
#' \describe{
#'     \item{X1}{row label, can be ignored}
#'     \item{id}{subject ID}
#'
#'     \item{a_0}{exposure measurement at time 0}
#'     \item{a_1}{exposure measurement at time 1}
#'     \item{a_2}{exposure measurement at time 2}
#'
#'     \item{l_0}{covariate measurement at time 0}
#'     \item{l_1}{covariate measurement at time 1}
#'     \item{l_2}{covariate measurement at time 2}
#'
#'     \item{m_0}{covariate measurement at time 0}
#'     \item{m_1}{covariate measurement at time 1}
#'     \item{m_2}{covariate measurement at time 2}
#'
#'     \item{n_0}{covariate measurement at time 0}
#'     \item{n_1}{covariate measurement at time 1}
#'     \item{n_2}{covariate measurement at time 2}
#'
#'     \item{o_0}{covariate measurement at time 0}
#'     \item{o_1}{covariate measurement at time 1}
#'     \item{o_2}{covariate measurement at time 2}
#'
#'     \item{p_0}{covariate measurement at time 0}
#'     \item{p_1}{covariate measurement at time 1}
#'     \item{p_2}{covariate measurement at time 2}
#'
#'     \item{wa_0}{inverse probability of exposure weight at time 0}
#'     \item{wa_1}{inverse probability of exposure weight at time 1}
#'     \item{wa_2}{inverse probability of exposure weight at time 2}
#'
#'     \item{wax_0}{cumulative inverse probability weight of exposure at time 0}
#'     \item{wax_1}{cumulative inverse probability weight of exposure at time 1}
#'     \item{wax_2}{cumulative inverse probability weight of exposure at time 2}
#'
#'     \item{wsx_0}{cumulative inverse probability of censoring weight at time 0}
#'     \item{wsx_1}{cumulative inverse probability of censoring weight at time 1}
#'     \item{wsx_2}{cumulative inverse probability of censoring weight at time 2}
#'
#'     \item{e5_0}{propensity score strata at time 0}
#'     \item{e5_1}{propensity score strata at time 1}
#'     \item{e5_2}{propensity score strata at time 2}
#'
#'     \item{h_0}{exposure history at time 0}
#'     \item{h_1}{exposure history at time 1}
#'     \item{h_2}{exposure history at time 2}
#'
#'     \item{s_0}{censoring indicator at time 0}
#'     \item{s_1}{censoring indicator at time 1}
#'     \item{s_2}{censoring indicator at time 2}
#'
#'     }
#'
"example_sml"
