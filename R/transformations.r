#
# transformations.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2012
# Contains: transform, donothing, msqrt,mlog, reciproce, mprobit, mlogit
#

donothing <- function(x, ...){ invisible(return(x)) }
msqrt     <- function(x, ...){ invisible(return(sqrt(x))) }
mlog      <- function(x, ...){ invisible(return(log(x, ...))) }
reciproce <- function(x, ...){ invisible(return(1/x)) }
mprobit   <- function(x, ...){ require(VGAM); invisible(return(probit(x, ...))) }
mlogit    <- function(x, ...){ require(VGAM); invisible(return(logit(x, ...))) }

transform <- function(matrix, transformations=c("nothing","log","sqrt","reciprocal","probit","logit"), ... , verbose=TRUE){
  options <- c("nothing","log","sqrt","reciprocal","probit","logit")
  chosen  <- pmatch(transformations, options)
  methods <- c(donothing,msqrt,mlog,reciproce,mprobit,mlogit)

  res <- vector("list",length(methods))
  idx <- 1
  for(x in chosen){
    cat("Applying",x,"a",options[x],"transformation to the data\n")
    res[[idx]] <- apply(matrix,1,function(t){
      (methods[x][[1]])(t, ...)
    })
  idx <- idx+1
  }
  transformations
}

