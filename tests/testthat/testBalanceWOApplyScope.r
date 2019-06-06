######################################
## TEST BALANCE WITHOUT APPLY SCOPE ##
######################################

context("balance function")

df.twdy <- toy_wide_dropoutY

#diagnostic 1

test_that("diagnostic 1, no censoring, scope all",{

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
    scope="all",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,17)
  name.cov <- rep("l",17)  %>% factor(levels=c("l","m","p"))
  H <- c("H",rep("H0",2),rep("H00",3),rep("H01",3),rep("H1",2),rep("H10",3),rep("H11",3))
  time.exposure <- c(0,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,2)
  time.covariate <- c(0,0,1,0,1,2,0,1,2,0,1,0,1,2,0,1,2)
  D <- c(0.02609263,0.15209607,0.21724389,0.05392655,0.11621380,0.17294715,0.20766703,0.08740416,0.34194962,0.10434636,0.24238544,-0.01839888,0.03426966,0.21418539,0.05650881,0.02276367,0.15625715)
  SMD <- c(0.05590150,0.33343944,0.46962221,0.12520235,0.27013382,0.37295726,0.43471684,0.17484151,0.74633007,0.22225472,0.51276500,-0.04110115,0.06881350,0.43824548,0.11523424,0.04989958,0.36601604)
  N <- c(1000, 511, 511, 282, 282, 282, 193, 193, 193, 489, 489, 169, 169, 169, 203, 203, 203)
  Nexp <- c(489, 209, 209,  85,  85,  85, 110, 110, 110, 265, 265,  80,  80,  80, 141, 141, 141)

  expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
  expect_equal(nrow(df.twdy.b),41)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)



})

test_that("diagnostic 1, censoring, scope all",{

  df.twdy.t <- lengthen(
    input=df.twdy,
    diagnostic=1,
    censoring="yes",
    id="uid",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone",
    temporal.covariate=c("l","m"),
    static.covariate="p",
    censor="s"
  )

  df.twdy.b <- balance(
    input=df.twdy.t,
    diagnostic=1,
    approach="none",
    censoring="yes",
    scope="all",
    times.exposure=c(0,1,2),
    times.covariate=c(0,1,2),
    exposure="a",
    history="haone"
  )

  #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
  #nrow(df.twdy.b)

  E <- rep(1,17)
  name.cov <- rep("l",17)  %>% factor(levels=c("l","m","p"))
  H <- c("H",rep("H0",2),rep("H00",3),rep("H01",3),rep("H1",2),rep("H10",3),rep("H11",3))
  time.exposure <- c(0,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,2)
  time.covariate <- c(0,0,1,0,1,2,0,1,2,0,1,0,1,2,0,1,2)
  D <- c(0.02609263,0.13873884,0.22327564,0.05386251,0.11139746,0.16770827,0.19251701,0.05605442,0.31510204,0.12189932,0.27740694,-0.05938834,0.02489331,0.18492176,-0.02593126,-0.05947726,0.11648487)
  SMD <- c(0.05590150,0.30519938,0.48385999,0.12382796,0.25878327,0.36035888,0.40310289,0.11231286,0.68107384,0.25901317,0.58458062,-0.13265952,0.04991489,0.37730425,-0.05216420,-0.13248591,0.28258688)
  N<-c(1000, 475, 475, 270, 270, 270, 173, 173, 173, 372, 372, 150, 150, 150, 156, 156, 156)
  Nexp<-c(489, 193, 193,  83,  83,  83,  98,  98,  98, 203, 203,  74,  74,  74, 113, 113, 113)

  expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
  expect_equal(nrow(df.twdy.b),41)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
  expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
  expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)

})

#diagnostic 2

test_that("diagnostic 2, no censoring, weight, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=2,
      censoring="no",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      history="haone",
      exposure="a",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      weight.exposure="wax"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=2,
      approach="weight",
      censoring="no",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      weight.exposure="wax"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,4)
    name.cov <- rep("l",4)  %>% factor(levels=c("l","m","p"))
    H <- c("H","H","H0","H1")
    time.exposure <- c(0,0,1,1)
    time.covariate <- c(1,2,2,2)
    D<-c(0.2704053,0.1706721,0.1935464,0.2826376)
    SMD<-c(0.5624117,0.3462512,0.4054861,0.6096874)
    N<-c(1000.3453, 846.7623, 488.2676, 358.0988)
    Nexp<-c(489.1822, 371.5281, 207.5916, 195.2960)

    expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),8)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)

  })

  test_that("diagnostic 2, censoring, weight, stratify, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=2,
      censoring="yes",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      censor="s",
      history="haone",
      weight.exposure="wax",
      weight.censor="wsx"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=2,
      approach="weight",
      censoring="yes",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      weight.exposure="wax",
      weight.censor="wsx"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,4)
    name.cov <- rep("l",4)  %>% factor(levels=c("l","m","p"))
    H <- c("H","H","H0","H1")
    time.exposure <- c(0,0,1,1)
    time.covariate <- c(1,2,2,2)
    D<-c(0.2374073,0.1692529,0.2284099,0.3389303)
    SMD<-c(0.4907694,0.3434338,0.4778534,0.7431715)
    N<-c(847.1980,728.3102,378.8492,336.3361)
    Nexp<-c(418.5980,361.0482,163.0554,174.2889)

    expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),8)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)

})

  test_that("diagnostic 2, no censoring, stratify, scope all",{

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
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      strata="e5"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,7)
    name.cov <- rep("l",7) %>% factor(levels=c("l","m","n"))
    S <- c(1,2,3,3,3,4,5)
    time.exposure <- c(1,1,0,0,1,1,1)
    time.covariate <- c(2,2,1,2,2,2,2)
    D<-c(0.2405892,0.1935185,0.2729881,0.1713809,0.2397481,0.2871219,0.2139037)
    SMD<-c(0.5127960,0.4022628,0.5677838,0.3476923,0.4914239,0.6232376,0.4875769)
    N<-c(180, 174,1000, 847, 173, 154, 166)
    Nexp<-c(39,  54, 489, 372,  82,  89, 132)

    expect_equal(names(df.twdy.b),c("E","S","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),14)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"S"],S)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)
  })

  test_that("diagnostic 2, censoring, stratify, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=2,
      censoring="yes",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      strata="e5",
      weight.censor="wsx",
      censor="s"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=2,
      approach="stratify",
      censoring="yes",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      strata="e5",
      weight.censor="wsx"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,7)
    name.cov <- rep("l",7) %>% factor(levels=c("l","m","n"))
    S <- c(1,2,3,3,3,4,5)
    time.exposure <- c(1,1,0,0,1,1,1)
    time.covariate <- c(2,2,1,2,2,2,2)
    D<-c(0.2237924,0.2001997,0.2400822,0.1704699,0.3568130,0.3308707,0.2447526)
    SMD<-c(0.4717475,0.4186004,0.4963027,0.3459036,0.7391421,0.7126800,0.5709485)
    N<-c(153.6263,143.5396,846.8752,728.0074,153.6360,128.5053,148.7002)
    Nexp<-c( 31.66276, 44.50004,418.29001,360.94672, 70.78702, 74.55275,118.33543)

    expect_equal(names(df.twdy.b),c("E","S","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),14)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"S"],S)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)


  })

  #diagnostic 3

  test_that("diagnostic 3, no censoring, weight, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=3,
      censoring="no",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      weight.exposure="wax"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=3,
      approach="weight",
      censoring="no",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      weight.exposure="wax",
      sort.order="alphabetical",
      loop="no",
      ignore.missing.metric="no",
      metric="SMD",
      sd.ref="no"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,17)
    name.cov <- rep("l",17)  %>% factor(levels=c("l","m","p"))
    H <- c("H",rep("H0",2),rep("H00",3),rep("H01",3),rep("H1",2),rep("H10",3),rep("H11",3))
    time.exposure <- c(0,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,2)
    time.covariate <- c(0,0,1,0,1,2,0,1,2,0,1,0,1,2,0,1,2)
    D<-c(0.0132437982, 0.0281118841, 0.0473004498,-0.0007154549, 0.0219660364,-0.0602237719, 0.0950079294,-0.0144350060, 0.1212943137, 0.0351917311, 0.0522696694,-0.0457638834,-0.0577264801,-0.1349804239, 0.0409262647,-0.0889460751,-0.0077944626)
    SMD<-c(0.028373838 , 0.061484448 , 0.101989826 ,-0.001662109 , 0.051121193 ,-0.129995531 , 0.198281755 ,-0.028867581 , 0.265780729 , 0.074945005 , 0.110601626 ,-0.102210629 ,-0.115926434 ,-0.276297658 , 0.083479654 ,-0.194916554 ,-0.018217258 )
    N<-c(1000.3453, 521.9669, 521.9669, 289.7041, 289.7041, 289.7041, 191.3426, 191.3426, 191.3426, 476.0281, 476.0281, 184.6182, 184.6182, 184.6182, 193.8664, 193.8664, 193.8664)
    Nexp<- c(489.18220, 222.78259, 222.78259,  84.95693,  84.95693,  84.95693, 113.88562, 113.88562, 113.88562, 258.93368, 258.93368,  85.55982,  85.55982,  85.55982, 132.49058, 132.49058, 132.49058)

    expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),41)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)

  })


  test_that("diagnostic 3, censoring, weight, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=3,
      censoring="yes",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      weight.exposure="wax",
      weight.censor="wsx",
      censor="s"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=3,
      approach="weight",
      censoring="yes",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      weight.exposure="wax",
      weight.censor="wsx",
      sort.order="alphabetical",
      loop="no",
      ignore.missing.metric="no",
      metric="SMD",
      sd.ref="no"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l"),])
    #nrow(df.twdy.b)

    E <- rep(1,17)
    name.cov <- rep("l",17)  %>% factor(levels=c("l","m","p"))
    H <- c("H",rep("H0",2),rep("H00",3),rep("H01",3),rep("H1",2),rep("H10",3),rep("H11",3))
    time.exposure <- c(0,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,2)
    time.covariate <- c(0,0,1,0,1,2,0,1,2,0,1,0,1,2,0,1,2)
    D<-c(0.013243798,0.015715645,0.055516907,0.006727681,0.021543744,-0.066020583,0.101234266,-0.100786889,0.077260268,0.057415674,0.088087932,-0.115198220,-0.071602909,-0.151319468,-0.085470009,-0.216231479,-0.072440449)#c(0.013243798, 0.016157784, 0.051980229,-0.002844623, 0.015815967,-0.068254599, 0.104956374,-0.076899766, 0.115652754, 0.050064019, 0.092905944,-0.067090916,-0.069583841,-0.170534731,-0.007140200,-0.178699307,-0.075917283)
    SMD<-c(0.02837384,0.03448724,0.11996607,0.01547963,0.05013113,-0.14204427,0.21187469,-0.20193135,0.16708899,0.12208466,0.18550120,-0.25686111,-0.14360121,-0.30916580,-0.17182266,-0.48300684,-0.17433090)#c(0.013243798,0.015715645,0.055516907,0.006727681,0.021543744,-0.066020583,0.101234266,-0.100786889,0.077260268,0.057415674,0.088087932,-0.115198220 -0.071602909,-0.151319468,-0.085470009,-0.216231479,-0.072440449) #c(0.028373838, 0.035459449, 0.112330984,-0.006544283, 0.036793102,-0.146820420, 0.219116045,-0.154032969, 0.250922669, 0.106379786, 0.195775471,-0.149630727,-0.139548490,-0.348361620,-0.014358615,-0.398632453,-0.183401291)
    N<-c(1000.3453,440.6936,440.6936,224.0152,224.0152,224.0152,147.8302,147.8302,147.8302,401.0093,401.0093,183.3202,183.3202,183.3202,178.0509,178.0509,178.0509)#c(1000.3453, 488.2676, 488.2676, 280.6058, 280.6058, 280.6058, 169.4771, 169.4771, 169.4771, 358.0988, 358.0988, 172.7356, 172.7356, 172.7356, 145.8214, 145.8214, 145.8214)
    Nexp<-c(489.18220,187.57279,187.57279,66.34637,66.34637,66.34637,84.29028,84.29028,84.29028,215.73767,215.73767,85.60948,85.60948,85.60948,121.38114,121.38114,121.38114)#c(489.18220, 207.59162, 207.59162,  83.51905,  83.51905,  83.51905, 100.29219, 100.29219, 100.29219, 195.29603, 195.29603,  81.27641,  81.27641,  81.27641, 102.34637, 102.34637, 102.34637)

    expect_equal(names(df.twdy.b),c("E","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),41)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"H"],H)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l"),"Nexp"],Nexp,tolerance=.0001,scale=1)

  })

  test_that("diagnostic 3, no censoring, stratify, scope all",{

    df.twdy.t <- lengthen(
      input=df.twdy,
      diagnostic=3,
      censoring="no",
      id="uid",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      temporal.covariate=c("l","m"),
      static.covariate="p",
      strata="e5"
    )

    df.twdy.b <- balance(
      input=df.twdy.t,
      diagnostic=3,
      approach="stratify",
      censoring="no",
      scope="all",
      times.exposure=c(0,1,2),
      times.covariate=c(0,1,2),
      exposure="a",
      history="haone",
      strata="e5",
      loop="no"
    )

    #print(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),])
    #nrow(df.twdy.b)

    E <- rep(1,13)
    S <- rep(1,13)
    H <- c(rep("H0",2),rep("H00",3),rep("H1",2),rep("H10",3),rep("H11",3))
    name.cov <- rep("l",13)  %>% factor(levels=c("l","m","p"))
    time.exposure <- c(1,1,2,2,2,1,1,2,2,2,2,2,2)
    time.covariate <- c(0,1,0,1,2,0,1,0,1,2,0,1,2)
    D<-c(-0.08079590, 0.09918601, 0.02222222, 0.18888889, 0.21111111,-0.12615385, 0.09076923,-0.03584229,-0.13261649, 0.34767025, 1.00000000,-0.66666667, 0.00000000)
    SMD<-c(-0.22004338,  0.25714206 , 0.05467434 , 0.50009388, 0.58817876,-0.28716705, 0.24566384,-0.08072658,-0.27792191, 0.96505533,NA,NA, 0.00000000)
    N<-c(138 ,138,108,108,108, 63, 63, 40, 40, 40,  4,  4,  4)
    Nexp<-c(31, 31  ,18  ,18,18,13,13, 9, 9, 9, 1, 1, 1)


    expect_equal(names(df.twdy.b),c("E","S","H","name.cov","time.exposure","time.covariate","D","SMD","N","Nexp"))
    expect_equal(nrow(df.twdy.b),186)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"S"],S)
    expect_equivalent(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"H"],H)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"name.cov"],name.cov)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"time.exposure"],time.exposure)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"time.covariate"],time.covariate)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1 ),"D"],D,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"SMD"],SMD,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"N"],N,tolerance=.0001,scale=1)
    expect_equal(df.twdy.b[which(df.twdy.b$name.cov=="l" & df.twdy.b$S==1),"Nexp"],Nexp,tolerance=.0001,scale=1)

  })
