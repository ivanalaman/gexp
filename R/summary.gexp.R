summary.gexp <- function(object,
                         digits=3L,
                         ...){

  cat("Database \n")

  print(object$dfm,
        digits = digits,
        ...) 

}
