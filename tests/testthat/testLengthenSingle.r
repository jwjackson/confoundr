####################################
## TEST LENGTHEN: SINGLE EXPOSURE ##
####################################

context("Lengthen Function: One Exposure")

df.twdn <- toy_wide_dropoutN

test_that("Diagnostic 1, no censoring",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=1,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==1),])

  uid <- rep("1",15)
  name.cov <- c(rep("l",6),rep("m",6),rep("p",3))
  time.exposure <- c(0,1,1,2,2,2,0,1,1,2,2,2,0,1,2)
  time.covariate <- c(0,0,1,0,1,2,0,0,1,0,1,2,0,0,0)
  haone <- c("H","H1","H1","H10","H10","H10","H","H1","H1","H10","H10","H10","H","H1","H10")
  a <- c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0)
  value.cov <- c(1,1,1,1,1,1,1,1,0,1,0,0,0,0,0)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov"))
  expect_equal(nrow(df.twdn.t),15000)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"value.cov"],value.cov)

})

test_that("Diagnostic 1, censoring",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=1,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    censor="s"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==31),])

  uid <- c("31","31","31","31","31","31","31","31")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  haone <- c("H","H1","H1","H","H1","H1","H","H1")
  a <-  c(1,0,0,1,0,0,1,0)
  value.cov <- c(0,0,0,0,0,0,0,0)
  censor <- c(0,0,0,0,0,0,0,0)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","censor"))
  expect_equal(nrow(df.twdn.t),12478)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"censor"],censor)

})

test_that("Diagnostic 2, no censoring, weight",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=2,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    weight.exposure="wax"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==1),])

  uid <- rep("1",6)
  name.cov <- c("l","l","l","m","m","m")
  time.exposure <- c(0,0,1,0,0,1)
  time.covariate <- c(1,2,2,1,2,2)
  haone <- c("H","H","H1","H","H","H1")
  a <- c(1,1,0,1,1,0)
  value.cov <- c(1,1,1,0,0,0)
  wax <- rep(c(0.946306049370567,0.946306049370567,1.08357917112966),2)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","wax"))
  expect_equal(nrow(df.twdn.t),6000)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"wax"],wax,tolerance=.000001,scale=1)

})

test_that("Diagnostic 2, censoring, weight",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=2,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    censor="s",
    weight.exposure="wax",
    weight.censor="wsx"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==31),])

  uid <- c("31","31")
  name.cov <- c("l","m")
  time.exposure <- c(0,0)
  time.covariate <- c(1,1)
  haone <- c("H","H")
  a <-  c(1,1)
  value.cov <- c(0,0)
  censor <- c(0,0)
  wax <- c(1.0144498,1.0144498)
  wsx <- c(0.9751602,0.9751602)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","censor","wax","wsx"))
  expect_equal(nrow(df.twdn.t),4690)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"censor"],censor)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wax"],wax,tolerance=.000001,scale=1)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wsx"],wsx,tolerance=.000001,scale=1)

})

test_that("Diagnostic 2, no censoring, stratification",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=2,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    strata="e5"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==1),])

  uid <- rep("1",6)
  name.cov <- c("l","l","l","m","m","m")
  time.exposure <- c(0,0,1,0,0,1)
  time.covariate <- c(1,2,2,1,2,2)
  haone <- c("H","H","H1","H","H","H1")
  a <- c(1,1,0,1,1,0)
  value.cov <- c(1,1,1,0,0,0)
  e5 <- c(3,3,4,3,3,4)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","e5"))
  expect_equal(nrow(df.twdn.t),6000)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"e5"],e5)

})

test_that("Diagnostic 2, censoring, stratification",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=2,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    censor="s",
    weight.censor="wsx",
    strata="e5"
  )

  uid <- c("31","31")
  name.cov <- c("l","m")
  time.exposure <- c(0,0)
  time.covariate <- c(1,1)
  haone <- c("H","H")
  a <-  c(1,1)
  value.cov <- c(0,0)
  censor <- c(0,0)
  wsx <- c(0.9751602,0.9751602)
  e5 <- c(3,3)

  #print(check <- df.twdn.t[which(df.twdn.t$uid==31),])

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","censor","wsx","e5"))
  expect_equal(nrow(df.twdn.t),4690)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"censor"],censor)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wsx"],wsx,tolerance=.000001,scale=1)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"e5"],e5)

})

test_that("Diagnostic 3, no censoring, weight",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=3,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    weight.exposure="wax"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==1),])

  uid <- rep("1",15)
  name.cov <- c(rep("l",6),rep("m",6),rep("p",3))
  time.exposure <- c(0,1,1,2,2,2,0,1,1,2,2,2,0,1,2)
  time.covariate <- c(0,0,1,0,1,2,0,0,1,0,1,2,0,0,0)
  haone <- c("H","H1","H1","H10","H10","H10","H","H1","H1","H10","H10","H10","H","H1","H10")
  a <- c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0)
  value.cov <- c(1,1,1,1,1,1,1,1,0,1,0,0,0,0,0)
  wax <- c(0.946306,1.083579,1.083579,1.593498,1.593498,1.593498,0.946306,1.083579,1.083579,1.593498,1.593498,1.593498,0.946306,1.083579,1.593498)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","wax"))
  expect_equal(nrow(df.twdn.t),15000)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"wax"],wax,tolerance=.000001,scale=1)

})

test_that("Diagnostic 3, censoring, weight",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=3,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    weight.exposure="wax",
    censor="s",
    weight.censor="wsx"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==31),])

  uid <- c("31","31","31","31","31","31","31","31")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  haone <- c("H","H1","H1","H","H1","H1","H","H1")
  a <-  c(1,0,0,1,0,0,1,0)
  value.cov <- c(0,0,0,0,0,0,0,0)
  censor <- c(0,0,0,0,0,0,0,0)
  wax <- c(1.0144498,0.7868815,0.7868815,1.0144498,0.7868815,0.7868815,1.0144498,0.7868815)
  wsx <- c(1,0.9751602,0.9751602,1,0.9751602,0.9751602,1,0.9751602)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","censor","wax","wsx"))
  expect_equal(nrow(df.twdn.t),12478)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wax"],wax,tolerance=.000001,scale=1)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wsx"],wsx,tolerance=.000001,scale=1)

})

test_that("Diagnostic 3, no censoring, stratification",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=3,
    censoring="no",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    strata="e5"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==1),])

  uid <- rep("1",15)
  name.cov <- c(rep("l",6),rep("m",6),rep("p",3))
  time.exposure <- c(0,1,1,2,2,2,0,1,1,2,2,2,0,1,2)
  time.covariate <- c(0,0,1,0,1,2,0,0,1,0,1,2,0,0,0)
  haone <- c("H","H1","H1","H10","H10","H10","H","H1","H1","H10","H10","H10","H","H1","H10")
  a <- c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0)
  value.cov <- c(1,1,1,1,1,1,1,1,0,1,0,0,0,0,0)
  e5 <- c(3,4,4,4,4,4,3,4,4,4,4,4,3,4,4)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","e5"))
  expect_equal(nrow(df.twdn.t),15000)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==1),"e5"],e5)

})

test_that("Diagnostic 3, censoring, stratification",{

  df.twdn.t <- lengthen(
    input=df.twdn,
    diagnostic=3,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    history="haone",
    censor="s",
    weight.censor="wsx",
    strata="e5"
  )

  #print(check <- df.twdn.t[which(df.twdn.t$uid==31),])

  uid <- c("31","31","31","31","31","31","31","31")
  name.cov <- c("l","l","l","m","m","m","p","p")
  time.exposure <- c(0,1,1,0,1,1,0,1)
  time.covariate <- c(0,0,1,0,0,1,0,0)
  haone <- c("H","H1","H1","H","H1","H1","H","H1")
  a <-  c(1,0,0,1,0,0,1,0)
  value.cov <- c(0,0,0,0,0,0,0,0)
  censor <- c(0,0,0,0,0,0,0,0)
  e5 <-  c(3,2,2,3,2,2,3,2)
  wsx <- c(1,0.9751602,0.9751602,1,0.9751602,0.9751602,1,0.9751602)

  expect_equal(names(df.twdn.t),c("uid","name.cov","time.exposure","time.covariate","haone","a","value.cov","censor","wsx","e5"))
  expect_equal(nrow(df.twdn.t),12478)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"uid"],uid)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"name.cov"],name.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.exposure"],time.exposure)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"time.covariate"],time.covariate)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"haone"],haone)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"a"],a)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"value.cov"],value.cov)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"censor"],censor)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"wsx"],wsx,tolerance=.000001,scale=1)
  expect_equal(df.twdn.t[which(df.twdn.t$uid==31),"e5"],e5)

})
