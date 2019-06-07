makeTreatments <- function(fl,
                           fe,
                           quali,
                           quanti,
                           posquanti)
{
  if (is.null(fl)) {#A função faz os fatores automaticamente!
    lfactors <- lapply(fe, 
                       function(x) as.matrix(x))

    names(lfactors) <- paste0("X", 
                              1:length(fe))

    lfaclevels <- as.list(tolower(names(lfactors)))

    lfacnumberlevels <- lapply(lfactors, 
                               function(x) 1:dim(x)[1])

    factors <- mapply(function(x, y) paste(x, y, sep = ""), 
                      lfaclevels, 
                      lfacnumberlevels, 
                      SIMPLIFY = FALSE)

    names(factors) <- names(lfactors)

  } else {
    if(!is.list(fl)){
      stop("This argument must be a list. See examples!")
    }
    if(length(fl) != length(fe)){
      stop("You must give names to all factors!")
    }

    if(quanti) {
      factors <- lapply(fl, 
                        as.ordered)
    } else if (!quanti & !quali) {
      # híbrido
      factors <- fl
      factors[posquanti] <- lapply(fl[posquanti], 
                                   as.ordered)
    } else {
      factors <- fl
    }
  }    
  return(factors)
}
