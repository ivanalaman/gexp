summary.gexp <- function(object,
                         digits = NULL,
                         ...){

  cat("Database \n")

  print(object$dfm,
        digits = digits,
        ...) 

}
