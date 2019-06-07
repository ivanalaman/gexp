makeXBeta <- function(cformula, 
                      dados, 
                      mu, 
                      fe, 
                      blke,
                      rowe, 
                      cole, 
                      inte, 
                      contrasts)
{
  X <- model.matrix(eval(parse(text = cformula)),
                    dados,
                    contrasts.arg = contrasts)

  if(length(mu) == 1){#univariado
    betas <- as.matrix(c(mu,
                         blke,
                         rowe,
                         cole,
                         unlist(fe),
                         inte))
  } else {#multivariado
    betas <- rbind(mu,
                   blke,
                   rowe,
                   cole,
                   do.call('rbind', 
                           fe),
                   inte)
  } 
  XB <- list(X = X,
             B = betas,
             XB = X%*%betas)
  return(XB)
}
