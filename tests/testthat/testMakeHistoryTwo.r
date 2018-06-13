##########################
## TEST MAKEHISTORY TWO ##
##########################

context("makehistory.two Function")

df.twdn   <- toy_wide_dropoutN
df.twdy   <- toy_wide_dropoutY

test_that("Formats and Names under No Dropout",{

  df.twdn.h <- makehistory.two(
    input=df.twdn,
    id="uid",
    times=c(0,1,2),
    exposure.a="a",
    exposure.b="s",
    name.history.a="ha",
    name.history.b="hs"
  )

  df.twdn.hg <- makehistory.two(
    input=df.twdn,
    id="uid",
    times=c(0,1,2),
    exposure.a="a",
    exposure.b="s",
    name.history.a="hag",
    name.history.b="hsg",
    group="p_0"
  )

  expect_equal(ncol(df.twdn.h),eval(ncol(df.twdn)+6))
  expect_equal(nrow(df.twdn.h),nrow(df.twdn))
  expect_equal(df.twdn.h$uid,df.twdn$uid)

  expect_equal(class(df.twdn.h$ha_0),"character")
  expect_equal(class(df.twdn.h$ha_1),"character")
  expect_equal(class(df.twdn.h$ha_2),"character")
  expect_equal(df.twdn.h$ha_0,df.twdn$hatwo_0)
  expect_equal(df.twdn.h$ha_1,df.twdn$hatwo_1)
  expect_equal(df.twdn.h$ha_2,df.twdn$hatwo_2)
  expect_equal(class(df.twdn.h$hs_0),"character")
  expect_equal(class(df.twdn.h$hs_1),"character")
  expect_equal(class(df.twdn.h$hs_2),"character")
  expect_equal(df.twdn.h$hs_0,df.twdn$hstwo_0)
  expect_equal(df.twdn.h$hs_1,df.twdn$hstwo_1)
  expect_equal(df.twdn.h$hs_2,df.twdn$hstwo_2)

  expect_equal(ncol(df.twdn.hg),eval(ncol(df.twdn)+6))
  expect_equal(nrow(df.twdn.hg),nrow(df.twdn))
  expect_equal(df.twdn.hg$uid,df.twdn$uid)

  expect_equal(class(df.twdn.hg$hag_0),"character")
  expect_equal(class(df.twdn.hg$hag_1),"character")
  expect_equal(class(df.twdn.hg$hag_2),"character")
  expect_equal(df.twdn.hg$hag_0,df.twdn$hatwog_0)
  expect_equal(df.twdn.hg$hag_1,df.twdn$hatwog_1)
  expect_equal(df.twdn.hg$hag_2,df.twdn$hatwog_2)
  expect_equal(class(df.twdn.hg$hsg_0),"character")
  expect_equal(class(df.twdn.hg$hsg_1),"character")
  expect_equal(class(df.twdn.hg$hsg_2),"character")
  expect_equal(df.twdn.hg$hsg_0,df.twdn$hstwog_0)
  expect_equal(df.twdn.hg$hsg_1,df.twdn$hstwog_1)
  expect_equal(df.twdn.hg$hsg_2,df.twdn$hstwog_2)

})

test_that("Formats and Names under Dropout",{

  df.twdy.h <- makehistory.two(
    input=df.twdy,
    id="uid",
    times=c(0,1,2),
    exposure.a="a",
    exposure.b="s",
    name.history.a="ha",
    name.history.b="hs"
  )

  df.twdy.hg <- makehistory.two(
    input=df.twdy,
    id="uid",
    times=c(0,1,2),
    exposure.a="a",
    exposure.b="s",
    name.history.a="hag",
    name.history.b="hsg",
    group="p_0"
  )

  expect_equal(ncol(df.twdy.h),eval(ncol(df.twdy)+6))
  expect_equal(nrow(df.twdy.h),nrow(df.twdy))
  expect_equal(df.twdy.h$uid,df.twdy$uid)

  expect_equal(class(df.twdy.h$ha_0),"character")
  expect_equal(class(df.twdy.h$ha_1),"character")
  expect_equal(class(df.twdy.h$ha_2),"character")
  expect_equal(df.twdy.h$ha_0,df.twdy$hatwo_0)
  expect_equal(df.twdy.h$ha_1,df.twdy$hatwo_1)
  expect_equal(df.twdy.h$ha_2,df.twdy$hatwo_2)
  expect_equal(class(df.twdy.h$hs_0),"character")
  expect_equal(class(df.twdy.h$hs_1),"character")
  expect_equal(class(df.twdy.h$hs_2),"character")
  expect_equal(df.twdy.h$hs_0,df.twdy$hstwo_0)
  expect_equal(df.twdy.h$hs_1,df.twdy$hstwo_1)
  expect_equal(df.twdy.h$hs_2,df.twdy$hstwo_2)

  expect_equal(ncol(df.twdy.hg),eval(ncol(df.twdy)+6))
  expect_equal(nrow(df.twdy.hg),nrow(df.twdy))
  expect_equal(df.twdy.hg$uid,df.twdy$uid)

  expect_equal(class(df.twdy.hg$hag_0),"character")
  expect_equal(class(df.twdy.hg$hag_1),"character")
  expect_equal(class(df.twdy.hg$hag_2),"character")
  expect_equal(df.twdy.hg$hag_0,df.twdy$hatwog_0)
  expect_equal(df.twdy.hg$hag_1,df.twdy$hatwog_1)
  expect_equal(df.twdy.hg$hag_2,df.twdy$hatwog_2)
  expect_equal(class(df.twdy.hg$hsg_0),"character")
  expect_equal(class(df.twdy.hg$hsg_1),"character")
  expect_equal(class(df.twdy.hg$hsg_2),"character")
  expect_equal(df.twdy.hg$hsg_0,df.twdy$hstwog_0)
  expect_equal(df.twdy.hg$hsg_1,df.twdy$hstwog_1)
  expect_equal(df.twdy.hg$hsg_2,df.twdy$hstwog_2)

})
