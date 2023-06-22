makeContrasts <- function(factors, 
                          quali, 
                          quanti, 
                          posquanti)
{
  if(quali){#only quantitatives
    contrast <- lapply(factors,
                       function(x)diag(length(x))) 
  } else if(quanti){#only quantitatives
    contrast <- lapply(factors,
                       function(x)contr.poly(length(x)))                            
  } else {#hybrido
    contrast <- lapply(factors,
                       function(x)diag(length(x)))
    contrast[posquanti] <- lapply(factors[posquanti],
                                  function(x)contr.poly(length(x)))   
  }
  return(contrast)
}
