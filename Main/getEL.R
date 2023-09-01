getEl	<- function( line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
  els	<- strsplit(line, sep)[[1]]
  
  if (ind[1] != -1) {
    if (fromEnd) {
      ind <- length(els)-(ind-1)
    }
  }
  
  if (final) {
    return( els[length(els)] )
  } else {
    
    if (reconstruct) {
      if (ex[1] > 0) {
        if (fromEnd) {
          ex <- length(els)-(ex-1)
        }
        ind <- setdiff((1:length(els)),ex)
      }
      
      newLine <- els[ind[1]]
      if (length(ind) > 1) {
        for (i in 2:length(ind)) {
          newLine <- paste(newLine, els[ind[i]], sep=sep)
        }
      }
      return ( newLine )
    } else {
      if ( ind[1] == -1 ) {
        return( els )
      } else {
        return( els[ind] )
      }
    }
  }
}
#calulate decimal date 
calcDecimalDate	<- function(day, month, year, defaultMonth=6, defaultDay=15) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334)
  
  if (month==0) {
    if (defaultMonth >= 1) {
      month <- defaultMonth
    } else {
      month	<- ceiling(runif(1)*12)
    }
  }
  
  if (day==0) {
    if (defaultDay >= 1) {
      day	<- defaultDay
    } else {
      day	<- ceiling(runif(1)*30)
    }
  }
  
  dd	<- cd[month] + day - 1
  
  decDate <- year + (dd/365)
  
  return ( decDate )
}


invertDecimalDate <- function( decDate, formatAsTxt=FALSE, ddmmyy=FALSE ) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
  fractD<- cd/365
  
  year		<- floor(decDate)
  fractYear 	<- decDate-year
  month		<- which(fractD >= fractYear)[1]-1
  
  if (month > 0) {
    fractMonth  <- fractYear-fractD[month]
    day		<- round((fractMonth*365)+1)
  } else {
    month <- 1
    day   <- 1
  }
  
  res <- c(day,month,year)
  
  if (formatAsTxt) {
    if (month < 10) {
      mm  <- paste("0",month,sep="")
    } else {
      mm <- month
    }
    res <- paste(year,mm,day,sep="-")
  }
  
  if (ddmmyy) {
    months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    if (day < 10) {
      dd <- paste("0",day,sep="")
    } else {
      dd <- day
    }
    res <- paste(dd,months[month],year,sep="-")
  }
  return( res )
  
}
