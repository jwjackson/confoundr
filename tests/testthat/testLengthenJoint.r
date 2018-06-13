###################################
## TEST LENGTHEN: JOINT EXPOSURE ##
###################################

context("Lengthen Function: Joint Exposure")

df.twdy <- toy_wide_dropoutY

test_that("Diagnostic 1, Exposure #1",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="hatwo"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1","1","1","1","1","1","1")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  hatwo <- c("H","H10","H10","H","H10","H10","H","H10")
  a <-  c(1,0,0,1,0,0,1,0)
  value.cov <- c(1,1,1,1,1,0,0,0)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hatwo","a","value.cov"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hatwo"],hatwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"a"],a)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)

})

test_that("Diagnostic 1, Exposure #2",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="s",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="hstwo"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1","1","1","1","1","1","1")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  hstwo <- c("H1","H100","H100","H1","H100","H100","H1","H100")
  s <-  c(0,1,1,0,1,1,0,1)
  value.cov <- c(1,1,1,1,1,0,0,0)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hstwo","s","value.cov"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hstwo"],hstwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"s"],s)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)

})

test_that("Diagnostic 2, Exposure #1, weight",{

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
    history="hatwo",
    weight.exposure="wax"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1")
  name.cov <- c("l","m")
  time.exposure <- c(0,0)
  time.covariate <- c(1,1)
  hatwo <- c("H","H")
  a <-  c(1,1)
  value.cov <-c(1,0)
  wax <- c(.946306,.946306)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hatwo","a","value.cov","wax"))
  expect_equal(nrow(df.twdy.t),5388)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hatwo"],hatwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"a"],a)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"wax"],wax,tolerance=.000001,scale=1)

})

test_that("Diagnostic 2, Exposure #1, stratification",{

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
    history="hatwo",
    strata="e5"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1")
  name.cov <- c("l","m")
  time.exposure <- c(0,0)
  time.covariate <- c(1,1)
  hatwo <- c("H","H")
  a <-  c(1,1)
  value.cov <- c(1,0)
  e5 <- c(3,3)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hatwo","a","value.cov","e5"))
  expect_equal(nrow(df.twdy.t),5388)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hatwo"],hatwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"a"],a)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)

})

test_that("Diagnostic 3, Exposure #1, weight",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=3,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="hatwo",
    weight.exposure="wax"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1","1","1","1","1","1","1")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  hatwo <- c("H","H10","H10","H","H10","H10","H","H10")
  a <-  c(1,0,0,1,0,0,1,0)
  value.cov <- c(1,1,1,1,1,0,0,0)
  wax <- c(0.946306,1.083579,1.083579,0.946306,1.083579,1.083579,0.946306,1.083579)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hatwo","a","value.cov","wax"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hatwo"],hatwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"a"],a)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"wax"],wax,tolerance=.000001,scale=1)

})

test_that("Diagnostic 3, Exposure #2, stratification",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=3,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="s",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="hstwo",
    strata="e5"
  )

  #print(check <- df.twdy.t[which(df.twdy.t$uid==1),])

  uid <- c("1","1","1","1","1","1","1","1")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  hstwo <- c("H1","H100","H100","H1","H100","H100","H1","H100")
  s <-  c(0,1,1,0,1,1,0,1)
  value.cov <- c(1,1,1,1,1,0,0,0)
  e5 <- c(3,4,4,3,4,4,3,4)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","hstwo","s","value.cov","e5"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"uid"],uid)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"hstwo"],hstwo)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"s"],s)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdy.t[which(df.twdy.t$uid==1),"e5"],e5)

})
