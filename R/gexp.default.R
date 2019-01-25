gexp.default <- function(mu        = NULL, 
                         err       = NULL,   
                         errp      = NULL,   
                         r         = 2L,  
                         factorsl  = NULL,
                         blocksl   = NULL,
                         rowsl     = NULL,
                         colsl     = NULL,
                         ef        = NULL, 
                         eint      = NULL, 
                         eb        = NULL, 
                         erow      = NULL, 
                         ecol      = NULL, 
                         contrasts = NULL,
                         type      = c('CRD','RCBD','LSD','FE','SPE'),
                         nrand     = 1,    
                         round     = 2L,   
                         random    = FALSE,  
                         ...)
{

  cl <- match.call()     

  switch(match.arg(type),
         CRD = {
           if(any(!is.null(eint) | !is.null(eb) | !is.null(erow) | !is.null(ecol))) stop("Only ef must be specified as effect of treatment!")

           result <- gexp.crd(mu        = mu, 
                              err       = err,   
                              r         = r,   
                              factorsl  = factorsl,  
                              ef        = ef,
                              contrasts = contrasts, 
                              round     = round, 
                              random    = random)
           class(result) <- c('gexp.crd','gexp','list')
         },
         RCBD= {
           if(is.null(eb)) stop("You must specify the block effect with eb argument!")

           result <- gexp.rcbd(mu        = mu,
                               err       = err,
                               r         = r,
                               factorsl  = factorsl,
                               blocksl   = blocksl,
                               ef        = ef,
                               eb        = eb,
                               contrasts = contrasts,
                               round     = round,
                               random    = random)
           class(result) <- c('gexp.rcbd','gexp','list')
         },
         LSD = {
           if(length(ef)!=1) stop("Use only one factor this design!")

           result <- gexp.lsd(mu        = mu, 
                              err       = err,
                              factorsl  = factorsl,
                              rowsl     = rowsl,
                              colsl     = colsl,
                              ef        = ef,  
                              erow      = erow, 
                              ecol      = ecol, 
                              nrand     = nrand,
                              contrasts = contrasts,
                              round     = round,
                              random    = random)
           class(result) <- c('gexp.lsd','gexp','list')
         },
         FE  = {
           result <- gexp.fe(mu        = mu,
                             err       = err, 
                             r         = r,
                             factorsl  = factorsl,
                             blocksl   = blocksl,
                             rowsl     = rowsl,
                             colsl     = colsl,
                             ef        = ef, 
                             eint      = eint, 
                             eb        = eb, 
                             erow      = erow, 
                             ecol      = ecol, 
                             nrand     = nrand, 
                             contrasts = contrasts,
                             round     = round, 
                             random    = random)
           class(result) <- c('gexp.fe','gexp','list')   
         },
         SPE = {
           result <- gexp.spe(mu        = mu,
                              err       = err,
                              errp      = errp,
                              r         = r,   
                              factorsl  = factorsl,
                              blocksl   = blocksl,
                              rowsl     = rowsl,
                              colsl     = colsl,
                              ef        = ef, 
                              eint      = eint, 
                              eb        = eb, 
                              erow      = erow, 
                              ecol      = ecol, 
                              nrand     = nrand,
                              contrasts = contrasts,
                              round     = round,
                              random    = random) 
           class(result) <- c('gexp.spe','gexp','list') 
         }
  )
  result$call <- cl  
  return(result)
}
