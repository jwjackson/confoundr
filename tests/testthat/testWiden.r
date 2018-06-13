context("widen Function")

df.twdn <- toy_wide_dropoutN
df.twdy <- toy_wide_dropoutY
df.tldn <- toy_long_dropoutN
df.tldy <- toy_long_dropoutY

test_that("Formats and Names under No Dropout",{

  df.tldn.w <- widen(
    input=df.tldn,
    id="uid",
    time="time",
    exposure="s",
    history="h",
    covariate=c("a","l","m","n","o","p"),
    weight.exposure=c("wa","wax"),
    weight.censor="wsx",
    strata="e5",
    censor="s"
  )

  df.tldn.w$h_0 <- NULL
  df.tldn.w$h_1 <- NULL
  df.tldn.w$h_2 <- NULL

 expect_equal(class(df.tldn.w$a_0),"numeric")         #check numeric format preservation
 expect_equal(eval(ncol(df.twdn)-18),ncol(df.tldn.w)) #check expected num columns
 expect_equal(nrow(df.twdn),nrow(df.tldn.w))          #check expected num rows
 expect_true(all(names(df.tldn.w %in% df.twdn)))      #check expected names

 })

test_that("Formats and Names under Dropout",{

   df.tldy.w <- widen(
     input=df.tldy,
     id="uid",
     time="time",
     exposure="s",
     history="h",
     covariate=c("a","l","m","n","o","p"),
     weight.exposure=c("wa","wax"),
     weight.censor="wsx",
     strata="e5",
     censor="s"
   )

   df.tldy.w$h_0 <- NULL
   df.tldy.w$h_1 <- NULL
   df.tldy.w$h_2 <- NULL

   expect_equal(class(df.tldy.w$a_0),"numeric")         #check numeric format preservation
   expect_equal(eval(ncol(df.twdy)-18),ncol(df.tldy.w)) #check expected num columns
   expect_equal(nrow(df.twdy),nrow(df.tldy.w))          #check expected num rows
   expect_true(all(names(df.tldy.w %in% df.twdy)))      #check expected names

})
