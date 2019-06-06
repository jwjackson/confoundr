#######################
## TEST OMIT.HISTORY ##
#######################

context("omit.history function")

df.twdy <- toy_wide_dropoutY

test_that("relative time omission",{

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

  df.twdy.o <- omit.history(
    input=df.twdy.t,
    omission="relative",
    distance=1,
    covariate.name="l"
  )

  #print(check <- df.twdy.o[which(df.twdy.o$uid==1),])

  uid <- c("1","1","1","1","1","1","1")
  name.cov <- c("l","l","m","m","m","p","p")
  time.exposure <- c(0,1,0,1,1,0,1)
  time.covariate <- c(0,1,0,0,1,0,0)
  haone <- c("H","H1","H","H1","H1","H","H1")
  a <-  c(1,0,1,0,0,1,0)
  value.cov <- c(1,1,1,1,0,0,0)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"uid"],uid)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"haone"],haone)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"a"],a)

})

test_that("fixed time omission",{

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

  df.twdy.o <- omit.history(
    input=df.twdy.t,
    omission="fixed",
    times=1,
    covariate.name="l"
  )

  #print(check <- df.twdy.o[which(df.twdy.o$uid==1),])

  uid <- c("1","1","1","1","1","1","1")
  name.cov <- c("l","l","m","m","m","p","p")
  time.exposure <- c(0,1,0,1,1,0,1)
  time.covariate <- c(0,0,0,0,1,0,0)
  haone <- c("H","H1","H","H1","H1","H","H1")
  a <-  c(1,0,1,0,0,1,0)
  value.cov <- c(1,1,1,1,0,0,0)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"uid"],uid)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"haone"],haone)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"a"],a)

})

test_that("same time omission",{

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

  df.twdy.o <- omit.history(
    input=df.twdy.t,
    omission="same.time",
    covariate.name="l"
  )

  #print(check <- df.twdy.o[which(df.twdy.o$uid==1),])

  uid <- c("1","1","1","1","1","1")
  name.cov <- c("l","m","m","m","p","p")
  time.exposure <- c(1,0,1,1,0,1)
  time.covariate <- c(0,0,0,1,0,0)
  haone <- c("H1","H","H1","H1","H","H1")
  a <-  c(0,1,0,0,1,0)
  value.cov <- c(1,1,1,0,0,0)

  expect_equal(names(df.twdy.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov"))
  expect_equal(nrow(df.twdy.t),13929)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"uid"],uid)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"haone"],haone)
  expect_equal(df.twdy.o[which(df.twdy.o$uid==1),"a"],a)

})

test_that("preserve name.cov format",{

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

  df.twdy.tfac <- df.twdy.t
  df.twdy.tfac$name.cov <- factor(df.twdy.t$name.cov,levels=unique(df.twdy.t$name.cov))

  df.twdy.o <- omit.history(
    input=df.twdy.tfac,
    omission="same.time",
    covariate.name="l"
  )

  expect_equal(class(df.twdy.o$name.cov),"factor")

})
