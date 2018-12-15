gexp.default <- function(mu         = NULL, 
                         sigma      = 1L,   
                         sigmap     = 1L,   
                         r          = 2L,  
                         newfactors = NULL,
                         ef         = NULL, 
                         einter     = NULL, 
                         eb         = NULL, 
                         erow       = NULL, 
                         ecol       = NULL, 
                         contrasts  = NULL,
                         type       = c('CRD','RCBD','LSD','FE','SPE'),
                         nrand      = 1,    
                         rd         = 2L,   
                         randomized = FALSE,  
                         ...)
{

  cl <- match.call()     

  switch(match.arg(type),
         CRD = {

           if(any(!is.null(einter) | !is.null(eb) | !is.null(erow) | !is.null(ecol))) stop("Only ef must be specified as effect of treatment!")

           result <- gexp.crd(mu         = mu, 
                              sigma      = sigma,   
                              r          = r,   
                              newfactors = newfactors,  
                              ef         = ef,
                              contrasts  = contrasts, 
                              rd         = rd, 
                              randomized = randomized)
           class(result) <- c('gexp.crd','gexp','list')
         },
         RCBD= {

           if(is.null(eb)) stop("You must specify the block effect with eb argument!")

           result <- gexp.rcbd(mu         = mu,
                               sigma      = sigma,
                               r          = r,
                               newfactors = newfactors, 
                               ef         = ef,
                               eb         = eb,
                               contrasts  = contrasts,
                               rd         = rd,
                               randomized = randomized)
           class(result) <- c('gexp.rcbd','gexp','list')
         },
         LSD = {

           if(length(ef)!=1) stop("Use only one factor this design!")

           result <- gexp.lsd(mu         = mu, 
                              sigma      = sigma,
                              newfactors = newfactors, 
                              ef         = ef,  
                              erow       = erow, 
                              ecol       = ecol, 
                              nrand      = nrand,
                              contrasts  = contrasts,
                              rd         = rd,
                              randomized = randomized)
           class(result) <- c('gexp.lsd','gexp','list')
         },
         FE  = {

           result <- gexp.fe(mu         = mu,
                             sigma      = sigma, 
                             r          = r,
                             newfactors = newfactors, 
                             ef         = ef, 
                             einter     = einter, 
                             eb         = eb, 
                             erow       = erow, 
                             ecol       = ecol, 
                             nrand      = nrand, 
                             contrasts  = contrasts,
                             rd         = rd, 
                             randomized = randomized)
           class(result) <- c('gexp.fe','gexp','list')   
         },
         SPE = {

           result <- gexp.spe(mu         = mu,
                              sigma      = sigma,
                              sigmap     = sigmap,
                              r          = r,   
                              newfactors = newfactors, 
                              ef         = ef, 
                              einter     = einter, 
                              eb         = eb, 
                              erow       = erow, 
                              ecol       = ecol, 
                              nrand      = nrand,
                              contrasts  = contrasts,
                              rd         = rd,
                              randomized = randomized) 
           class(result) <- c('gexp.spe','gexp','list') 
         }
  )
  result$call <- cl  
  return(result)
}
