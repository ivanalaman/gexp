print.gexp <- function(x,
                       digits = 3L,
                       ...){
  cat("Design Matrix \n")

  print(x$X,
        digits = digits,
        ...)

  if(!is.null(x$Z)){ 
    cat("\n Plot Design Matrix \n")

    print(x$Z,
          digits = digits,
          ...)  
  }

  cat("\n Response Variable \n")

  print(x$Y,
        digits = digits,
        ...)  

  cat("\n Database \n")

  print(x$dfm,
        digits = digits,
        ...)  

}
