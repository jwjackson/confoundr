#' @title Simulated data, loosely based on the Clinical Antipsychotic Trial of Intervention Effectiveness (CATIE) study.
#'
#' @description 1430 completely hypothetical persons with schizophrenia
#' randomized to one of five antipsychotics and followed for
#' up to 18 months. Note that the simulation did not build
#' in sequential randomization, as was done in the trial.
#'
#' @docType data
#' @usage data(catie_sim)
#'
#' @format A dataframe with 8,632 rows and 62 variables:
#' \describe{
#'      \item{CATIEID}{person id}
#'      \item{time}{month of study visit (0, 1, 3, 6, 9, 12, 15, 18)}
#'
#'      \item{td}{tardive diskinesia stratum}
#'      \item{zprcort}{ziprasidone cohort stratum}
#'      \item{race}{1:white, 2:black, 3:other}
#'      \item{age.grp}{1:18-24, 2:25-34, 3:45-44, 4:45-54, 5:55-67}
#'      \item{educ.bin}{high school graduate}
#'      \item{site.ro}{study site, research only}
#'      \item{site.sh}{study site, state mental health}
#'      \item{site.uc}{study site, university center}
#'      \item{site.va}{study site, veterans affairs}
#'      \item{treat.grp}{antipsychotic 1:ola, 2:que, 3:ris, 4:per, 5:zip}
#'
#'      \item{cs14}{drug use scale}
#'      \item{cs16}{clinical global impressions (CGI) severity scale}
#'      \item{calg1}{calgary depression scale}
#'      \item{weight}{in lbs}
#'      \item{epsmean}{Simpson-Agnes extrapyramidal symptoms}
#'      \item{qoltot}{quality of life total score}
#'      \item{pansstotal}{positive and negative syndrome scale (PANSS)}
#'      \item{phase.change.vis}{switch to new antipsychotic}
#'
#'      \item{white}{race dummy variable, white}
#'      \item{black}{race dummy variable, black}
#'      \item{other}{race dummy variable, other}
#'      \item{age.grp.1824}{age group dummy variable, 18-24 years}
#'      \item{age.grp.2534}{age group dummy variable, 25-34 years}
#'      \item{age.grp.3544}{age group dummy variable, 35-44 years}
#'      \item{age.grp.4554}{age group dummy variable, 45-54 years}
#'      \item{age.grp.5567}{age group dummy variable, 55-67 years}
#'      \item{Bpansstotal}{pansstotal at time 0}
#'      \item{Bcs14}{cs14 at time 0}
#'      \item{Bcs16}{cs16 at time 0}
#'      \item{Bcalg1}{calg1 at time 0}
#'      \item{Bqoltot}{qoltot at time 0}
#'
#'      \item{Chg.pansstotal}{change in pansstotal}
#'      \item{pct.gain}{percent weight gain}
#'      \item{phase.change.cum}{number of switches to antipsychotic}
#'      \item{phase.change.cum.rec}{time-varying version of ever switch to antipsychotic}
#'      \item{lead.pansstotal}{pansstotal at next visit}
#'
#'      \item{treat.grp.ola}{olanzapine arm dummy}
#'      \item{treat.grp.que}{quetiapine arm dummy}
#'      \item{treat.grp.ris}{risperidone arm dummy}
#'      \item{treat.grp.per}{perphenazine arm dummy}
#'      \item{treat.grp.zip}{ziprasidone arm dummy}
#'
#'      \item{studydisc}{last visit (1=yes, 0 otherwise)}
#'
#'      \item{num.x}{probability of treatment arm}
#'      \item{den.x}{probability of treatment arm given baseline covariates}
#'      \item{wx.b}{stabilized iptw for treatment arm}
#'      \item{num.po}{probability studydisc=1, given treat.grp & baseline covariates, common model}
#'      \item{den.po}{probability studydisc=1, given treat.grp & baseline and time-varying covariates, common model}
#'      \item{num.tr}{probability studydisc=1, given treat.grp & baseline covariates, treat.grp specific model}
#'      \item{den.tr}{probability studydisc=1, given treat.grp & baseline and time-varying covariates, treat.grp specific model}
#'
#'      \item{wpo}{stabilized ipcw, from common model, not truncated}
#'      \item{wtr}{stabilized ipcw, from treat.grp specific model, not truncated}
#'      \item{wpo}{stabilized ipcw, truncated 99th tile}
#'      \item{wtr}{stabilized ipcw, from treat.grp specific model, truncated 99th tile}
#'      \item{wpo}{stabilized ipcw, from common model, truncated 95th tile}
#'      \item{wtr}{stabilized ipcw, from treat.grp specific model, truncated 95th tile}
#'      \item{wpo}{stabilized ipcw, from common model, truncated 90th tile}
#'      \item{wtr}{stabilized ipcw, from treat.grp specific model, truncated 90th tile}
#'
#'      }
#'
#'
#' @keywords datasets
#'
#' @references Lieberman JA, Stroup TS, McEvoy JP, Swartz MS,
#' Rosenheck RA, Perkins DO, Keefe RS, Davis SM, Davis CE,
#' Lebowitz BD, Severe J, Hsiao JK; Clinical Antipsychotic Trials
#' of Intervention Effectiveness (CATIE) Investigators. Effectiveness
#' of antipsychotic drugs in patients with chronic schizophrenia.
#' N Engl J Med. 2005 Sep 22;353(12):1209-23. Epub 2005 Sep 19.
#' Erratum in: N Engl J Med. 2010 Sep 9;363(11):1092-3. PubMed PMID: 16172203.
#'
"catie_sim"
