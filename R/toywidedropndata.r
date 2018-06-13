#' Artifical data set used to test the functionality of confoundr.
#'
#' The toy_wide_censN data set contains 1,000 records and 52
#' variables. These variables include time-varying exposures,
#' outcomes, and covariates, along with strata and censoring
#' indicators. Time-varying inverse-probability-of-exposure
#' weights and censoring weights are available as well.
#'
#' @docType data
#' @usage data(toy_wide_dropoutN)
#'
#' @format A data frame with 1,000 rows and 52 variables:
#' \describe{
#'
#'     \item{uid}{subject ID}
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
#'     \item{s_0}{censoring indicator at time 0}
#'     \item{s_1}{censoring indicator at time 1}
#'     \item{s_2}{censoring indicator at time 2}
#'
#'     \item{haone_0}{exposure history at time 0}
#'     \item{haone_1}{exposure history at time 1}
#'     \item{haone_2}{exposure history at time 2}
#'
#'     \item{haoneg_0}{grouped by p_0, exposure history at time 0}
#'     \item{haoneg_1}{grouped by p_0, exposure history at time 1}
#'     \item{haoneg_2}{grouped by p_0, exposure history at time 2}
#'
#'     \item{hatwo_0}{a joint history given a,s at time 0}
#'	   \item{hatwo_1}{a joint history given a,s at time 1}
#'     \item{hatwo_2}{a joint history given a,s at time 2}
#'
#'     \item{hatwog_0}{grouped by p_0, a joint history given a,s at time 0}
#'	   \item{hatwog_1}{grouped by p_0, a joint history given a,s at time 1}
#'     \item{hatwog_2}{grouped by p_0, a joint history given a,s at time 2}
#'
#'     \item{hstwo_0}{s joint history given a,s at time 0}
#'	   \item{hstwo_1}{s joint history given a,s at time 1}
#'     \item{hstwo_2}{s joint history given a,s at time 2}
#
#'     \item{hstwog_0}{grouped by p_0, s joint history given a,s at time 0}
#'	   \item{hstwog_1}{grouped by p_0, s joint history given a,s at time 1}
#'     \item{hstwog_2}{grouped by p_0, s joint history given a,s at time 2}
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
#'
#'     }
#'
"toy_wide_dropoutN"
