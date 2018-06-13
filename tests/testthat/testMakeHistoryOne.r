##########################
## TEST MAKEHISTORY ONE ##
##########################

context("makehistory.one Function")

df.twdn <- toy_wide_dropoutN
df.twdy <- toy_wide_dropoutY

test_that("Formats and Names under No Dropout",{

df.twdn.h <- makehistory.one(
  input=df.twdn,
  id="uid",
  times=c(0,1,2),
  exposure="a",
  name.history="ha"
)

df.twdn.hg <- makehistory.one(
  input=df.twdn,
  id="uid",
  times=c(0,1,2),
  exposure="a",
  name.history="hag",
  group="p_0"
)

expect_equal(nrow(df.twdn.h),nrow(df.twdn))
expect_equal(df.twdn.h$uid,df.twdn$uid)

expect_equal(class(df.twdn.h$ha_0),"character")
expect_equal(class(df.twdn.h$ha_1),"character")
expect_equal(class(df.twdn.h$ha_2),"character")
expect_equal(df.twdn.h$ha_0,df.twdn$haone_0)
expect_equal(df.twdn.h$ha_1,df.twdn$haone_1)
expect_equal(df.twdn.h$ha_2,df.twdn$haone_2)

expect_equal(class(df.twdn.hg$hag_0),"character")
expect_equal(class(df.twdn.hg$hag_1),"character")
expect_equal(class(df.twdn.hg$hag_2),"character")
expect_equal(df.twdn.hg$hag_0,df.twdn$haoneg_0)
expect_equal(df.twdn.hg$hag_1,df.twdn$haoneg_1)
expect_equal(df.twdn.hg$hag_2,df.twdn$haoneg_2)

})

test_that("Formats and Names under Dropout",{

  df.twdy.h <- makehistory.one(
    input=df.twdy,
    id="uid",
    times=c(0,1,2),
    exposure="a",
    name.history="ha"
  )

  df.twdy.hg <- makehistory.one(
    input=df.twdy,
    id="uid",
    times=c(0,1,2),
    exposure="a",
    name.history="hag",
    group="p_0"
  )

  expect_equal(ncol(df.twdy.h),eval(ncol(df.twdy)+3))
  expect_equal(nrow(df.twdy.h),nrow(df.twdy))
  expect_equal(df.twdy.h$uid,df.twdy$uid)

  expect_equal(class(df.twdy.h$ha_0),"character")
  expect_equal(class(df.twdy.h$ha_1),"character")
  expect_equal(class(df.twdy.h$ha_2),"character")
  expect_equal(df.twdy.h$ha_0,df.twdy$haone_0)
  expect_equal(df.twdy.h$ha_1,df.twdy$haone_1)
  expect_equal(df.twdy.h$ha_2,df.twdy$haone_2)

  expect_equal(class(df.twdy.hg$hag_0),"character")
  expect_equal(class(df.twdy.hg$hag_1),"character")
  expect_equal(class(df.twdy.hg$hag_2),"character")
  expect_equal(df.twdy.hg$hag_0,df.twdy$haoneg_0)
  expect_equal(df.twdy.hg$hag_1,df.twdy$haoneg_1)
  expect_equal(df.twdy.hg$hag_2,df.twdy$haoneg_2)

})
