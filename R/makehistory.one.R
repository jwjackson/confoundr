#' Create exposure history for a single time varying exposure
#'
#' @param input this is a matrix of
#' @param id vector of ids
#' @param times vector of times
#' @param group grouping variable (optional)
#' @param exposure
#' @param name.history
#'
#' @return A \code{\link{data.frame}} of the output stuff
#' @export
#'
#' @examples
#' x = 5
#' \dontrun{
#'
#' }
#' @import dplyr
#' @importFrom tidyr gather_ spread
makehistory.one <- function (input,id,times,group=NULL,exposure,name.history="h") {

  if (is.null(id)) {
    stop ("ERROR: 'id' is missing. Please specify a unique identifier")
  }

  if (is.null(input)) {
    stop ("ERROR: 'input' dataframe is missing")
  }

  if (is.null(exposure)) {
    stop ("ERROR: root name for exposure is missing")
  }

  if (is.null(times)) {
    stop ("ERROR: indices for exposure measurement times is missing. Please specify a numeric vector of times")
  }

  list.exposure <- paste(exposure,times,sep="_")

  if (any(!list.exposure %in% names(input))) {
    stop("ERROR: The exposure root name is misspelled, or some exposure measurements are missing from the input dataframe, or incorrect measurement times have been specified")
  }

  if (!all(!input$id %in% input$id[duplicated(input$id)])) {
    stop("ERROR: id does not uniquely identify each observation (i.e. each row). Please specify a unique identifier.")
  }

  cumpaste = function(x,.sep="") {
    Reduce(function(x1, x2) paste(x1,x2,sep=.sep),x,accumulate=TRUE)
  }

  if (is.null(group)) {

    input.temp <- input %>% ungroup() %>% select_(.dots=c(id,list.exposure)) %>%
      gather_(key_col="exp.name.time",value_col="exp.value",gather_cols=list.exposure) %>%
      separate_(col="exp.name.time",into=c("exp.name","exp.time"),sep="_") %>%
      mutate(exp.time=as.numeric(exp.time)) %>%
      arrange_(id,"exp.time") %>%
      group_by_(id) %>%
      mutate(his.name=name.history,
             his.time=exp.time,
             his.value=ifelse(his.time==first(his.time),
                              "H",
                              paste("H",lag(cumpaste(exp.value)),sep=""))
      ) %>%
      select_(.dots=c(id,"his.name","his.time","his.value")) %>%
      unite_(col="his.name.time",from=c("his.name","his.time"),sep="_") %>%
      spread(his.name.time,his.value)

    output <- left_join(input,input.temp,by=id)
    data.frame(output)

  } else if (!is.null(group)) {

    input.temp <- input %>% ungroup() %>% select_(.dots=c(id,group,list.exposure)) %>% rename_("GROUP"=group) %>%
      gather_(key_col="exp.name.time",value_col="exp.value",gather_cols=list.exposure) %>%
      separate_(col="exp.name.time",into=c("exp.name","exp.time"),sep="_") %>%
      mutate(exp.time=as.numeric(exp.time)) %>%
      arrange_(id,"exp.time") %>%
      group_by_(id) %>%
      mutate(his.name=name.history,
             his.time=exp.time,
             his.value=ifelse(his.time==first(his.time),
                              paste(GROUP,"H",sep=""),
                              paste(GROUP,"H",lag(cumpaste(exp.value)),sep=""))
      ) %>%
      select_(.dots=c(id,"his.name","his.time","his.value")) %>%
      unite_(col="his.name.time",from=c("his.name","his.time"),sep="_") %>%
      spread(his.name.time,his.value)

    output <- left_join(input,input.temp,by=id)
    data.frame(output)

  }

  return(output)

}
