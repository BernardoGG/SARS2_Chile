calcDecimalDate_fromTxt	<- function( dateTxt, sep="/", namedMonths=FALSE, dayFirst=FALSE) {
  els 	<- strsplit(dateTxt, sep)[[1]]
  if (dayFirst) {
    if (length(els) > 1) {
      els <- els[length(els):1]
    }
  }
  
  year 	<- as.integer(els[1])
  
  if (length(els)==1) {
    month <- 6  #7
    day	<- 15 #2
    decDate <- year + 0.5
  } else {
    
    if (length(els)==2) {
      if (nchar(els[2]) > 0) {
        if (namedMonths) {
          month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
        } else {
          month <- as.integer(els[2])
        }
        day	<- 15
        decDate <- calcDecimalDate(day,month,year)
      } else {
        month <- 6 #7
        day   <- 15 #2
        decDate <- year + 0.5
      }
    } else {
      if (namedMonths) {
        month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
      } else {
        month <- as.integer(els[2])
      }
      
      if (nchar(els[3]) > 0) {
        day 	<- as.integer(els[3])
      } else {
        day <- 15
      }
      decDate <- calcDecimalDate(day,month,year)
    }
  }
  
  
  return ( decDate )
}


calcDecimalDate_from_yymmdd	<- function( dateTxt, sep="/", ycutoff=15, defaultMonth=6, defaultDay=15 ) {
  els	<- strsplit(dateTxt, sep)[[1]]
  yy	<- as.integer(els[1])
  mm	<- as.integer(els[2])
  dd	<- as.integer(els[3])
  
  if (!is.finite(yy)) {
    return( -1 )
  } else {
    if (yy <= ycutoff) {
      yy <- yy+2000
    }
    if ((yy > ycutoff) & (yy < 99)) {
      yy <- yy+1900
    }
    
    if (!is.finite(mm)) {
      mm <- 0
    }
    if (!is.finite(dd)) {
      dd <- 0
    }
    return ( calcDecimalDate( dd, mm, yy, defaultMonth=defaultMonth, defaultDay=defaultDay ) )
  }
  
}