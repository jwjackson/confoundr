###################################
## TEST BALANCE WITH APPLY SCOPE ##
###################################

context("balance function")

df.twdy <- toy_wide_dropoutY

test_that("diagnostic 1, recent, recency 1", {

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="no",
    scope="recent",
    recency=1,
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,6)
  name.cov <- rep("l",6) %>% factor(levels=c("l","m","n"))
  H <- c("H0","H00","H01","H1","H10","H11")
  time.exposure <- c(1,2,2,1,2,2)
  time.covariate <- c(0,1,1,0,1,1)
  D<-c(0.15209607,0.11621380,0.08740416,0.10434636,0.03426966,0.02276367)
  SMD<-c(0.33343944,0.27013382,0.17484151,0.22225472,0.06881350,0.04989958)
  N<-c(511,282,193,489,169,203)
  Nexp<-c(209, 85,110,265, 80,141)

  expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
  expect_equal(nrow(df.twdy.b),14)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"Nexp"],Nexp,tolerance=.0001,scale=1)

})

test_that("diagnostic 2, average over strata",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=2,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    strata="e5"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=2,
    approach="stratify",
    censoring="no",
    scope="average",
    average.over="strata",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    strata="e5"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,3)
  name.cov <- rep("l",3) %>% factor(levels=c("l","m","n"))
  time.exposure <- c(0,0,1)
  time.covariate <- c(1,2,2)
  D<-c(0.2729881,0.1713809,0.2339781)
  SMD<-c(0.5677838,0.3476923,0.5008615)
  Nexp<-c(489,372,396)
  N<-c(1000, 847, 847)

  expect_equal(names(df.twdy.b),c("time.exposure","time.covariate","name.cov","D","SMD","Nexp","N"))
  expect_equal(nrow(df.twdy.b),6)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"Nexp"],Nexp,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})

test_that("diagnostic 1, average over history",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="no",
    scope="average",
    average.over="history",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,6)
  name.cov <- rep("l",6) %>% factor(levels=c("l","m","n"))
  time.exposure <- c(0,1,1,2,2,2)
  time.covariate <- c(0,0,1,0,1,2)
  D<-c(0.02609263,0.12874646,0.22953811,0.07514628,0.07090188,0.21568466)
  SMD<-c(0.0559015,0.2790701,0.4907190,0.1601581,0.1554678,0.4693983)
  N<-c(1000,1000,1000, 847, 847, 847)

  expect_equal(names(df.twdy.b),c("time.exposure","time.covariate","name.cov","D","SMD","N"))
  expect_equal(nrow(df.twdy.b),15)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})

test_that("diagnostic 1, average over time",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="no",
    scope="average",
    average.over="time",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,3)
  name.cov <- rep("l",3) %>% factor(levels=c("l","m","n"))
  distance <- c(0,1,2)
  D<-c(0.15395702,0.10222001,0.07514628)
  SMD<-c(0.3316477,0.2223884,0.1601581)
  N<-c(2847,1847,847)

  expect_equal(names(df.twdy.b),c("distance","name.cov","D","SMD","N"))
  expect_equal(nrow(df.twdy.b),9)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"distance"],distance)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})

test_that("diagnostic 1, average within periods of distance",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="no",
    scope="average",
    average.over="distance",
    periods=list(0:1,2),
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  period.id <- c(1,2)
  period.start <- c(0,2)
  period.end <- c(1,2)
  name.cov <- c("l","l") %>% factor(levels=c("l","m","n"))
  D<- c(0.13359949,0.07514628)
  SMD<-c(0.2886562,0.1601581)
  N<-c(4694,847)

  expect_equal(names(df.twdy.b),c("period.id","period.start","period.end","name.cov","D","SMD","N"))
  expect_equal(nrow(df.twdy.b),6)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})

test_that("diagnostic 1, average over distance",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="no",
    scope="average",
    average.over="distance",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  period.id <- 1
  period.start <- 0
  period.end <- 2
  name.cov <- c("l") %>% factor(levels=c("l","m","n"))
  D <- 0.1246643
  SMD <- 0.2690139
  N <- 5541

  expect_equal(names(df.twdy.b),c("period.id","period.start","period.end","name.cov","D","SMD","N"))
  expect_equal(nrow(df.twdy.b),3)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" ),"N"],N,tolerance=.0001,scale=1)

})
