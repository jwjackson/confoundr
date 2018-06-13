###################
## TEST DIAGNOSE ##
###################

context("diagnose function")

df.twdy <- toy_wide_dropoutY

test_that("diagnostic 1, scope recent, recency 0",{

  df.twdy.d <- diagnose(
    input=df.twdy,
    id="uid",
    diagnostic=1,
    scope="recent",
    recency=0,
    censoring="no",
    exposure="a",
    temporal.covariate=c("l","m","n"),
    static.covariate=c("o","p"),
    history="haone",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2)
  )

  #print(df.twdy.d[which(df.twdy.d$name.cov=="l"),])
  #nrow(df.twdy.d)

  E <- rep(1,23)
  H <- c(rep("H",5),rep("H0",3),rep("H00",3),rep("H01",3),rep("H1",3),rep("H10",3),rep("H11",3))
  name.cov <- c("l","m","n","o","p",rep(c("l","m","n"),6)) %>% factor(levels=c("l","m","n","o","p"))
  time.exposure <- c(rep(0,5),rep(1,3),rep(2,6),rep(1,3),rep(2,6))
  time.covariate <- c(rep(0,5),rep(1,3),rep(2,6),rep(1,3),rep(2,6))
  D<-c(0.02609263,0.21724389,0.17294715,0.34194962,0.24238544,0.21418539,0.15625715)
  SMD<-c(0.0559015,0.4696222,0.3729573,0.7463301,0.5127650,0.4382455,0.3660160)
  N<-c(1000, 511, 282, 193, 489, 169, 203)
  Nexp<-c(489, 209,  85, 110, 265,  80, 141)

  expect_equal(names(df.twdy.d),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
  expect_equal(nrow(df.twdy.d),23)
  expect_equivalent(df.twdy.d$name.cov,name.cov)
  expect_equivalent(df.twdy.d$time.exposure,time.exposure)
  expect_equivalent(df.twdy.d$time.covariate,time.covariate)
  expect_equal(df.twdy.d[which(df.twdy.d$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.d[which(df.twdy.d$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.d[which(df.twdy.d$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})
