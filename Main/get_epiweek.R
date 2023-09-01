# uses package lubridate
get_week <- function(dateTxt, useISOWeek=TRUE) {
  if (useISOWeek) {
    return( isoweek(ymd(dateTxt)) )
  } else {
    return( week(ymd(dateTxt)) )
  }
}


# uses package lubridate
get_epiweek <- function(dateTxt) {
  return( epiweek(ymd(dateTxt)) )
}
