makeContrasts <- function(factors, 
                          quali, 
                          quanti, 
                          posquanti)
{
  if(quali){#S� qualitativos
    contrast <- lapply(factors,
                       function(x)diag(length(x))) 
  } else if(quanti){#S� quantitativos
    contrast <- lapply(factors,
                       function(x)contr.poly(length(x)))                            
  } else {#h�brido
    contrast <- lapply(factors,
                       function(x)diag(length(x)))
    contrast[posquanti] <- lapply(factors[posquanti],
                                  function(x)contr.poly(length(x)))   
  }
  return(contrast)
}
