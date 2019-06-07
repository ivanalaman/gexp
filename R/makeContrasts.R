makeContrasts <- function(factors, 
                          quali, 
                          quanti, 
                          posquanti)
{
  if(quali){#Só qualitativos
    contrast <- lapply(factors,
                       function(x)diag(length(x))) 
  } else if(quanti){#Só quantitativos
    contrast <- lapply(factors,
                       function(x)contr.poly(length(x)))                            
  } else {#híbrido
    contrast <- lapply(factors,
                       function(x)diag(length(x)))
    contrast[posquanti] <- lapply(factors[posquanti],
                                  function(x)contr.poly(length(x)))   
  }
  return(contrast)
}
