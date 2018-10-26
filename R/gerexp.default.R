gerexp.default <- function(mu     = NULL, 
                           sigma  = 1L,   
                           sigmap = 1L,   
                           r      = 2L,   
                           ef     = NULL, 
                           einter = NULL, 
                           eb     = NULL, 
                           erow   = NULL, 
                           ecol   = NULL, 
                           type   = c('CRD','RCBD','LSD','FE','SPE'),
                           nrand  = 1,    
                           rd     = 2L,   
                           randomized = FALSE,  
                           ...)
{

  switch(match.arg(type),
         CRD = {

           if(any(!is.null(einter) | !is.null(eb) | !is.null(erow) | !is.null(ecol))) stop("Only f must be specified as effect of treatment!")

           result <- gerexp.crd(mu = mu, 
                                sigma = sigma,   
                                r = r,   
                                ef = ef, 
                                rd = rd,   
                                randomized = randomized)
          class(result) <- c('gerexp.crd','gerexp','list')
         },
         RCBD= {

           if(is.null(eb)) stop("You must specify the block effect with eb argument!")

           result <- gerexp.rcbd(mu = mu,
                                 sigma = sigma,
                                 ef = ef,
                                 eb = eb,
                                 rd = rd,
                                 randomized = randomized)
           class(result) <- c('gerexp.rcbd','gerexp','list')
         },
         LSD = {

           if(length(ef)!=1) stop("Use only one factor this design!")

           result <- gerexp.lsd(mu = mu, 
                                sigma = sigma, 
                                ef = ef,  
                                erow = erow, 
                                ecol = ecol, 
                                nrand = nrand,    
                                rd = rd,   
                                randomized = randomized)
           class(result) <- c('gerexp.lsd','gerexp','list')
         },
         FE  = {

           result <- gerexp.fe(mu = mu,
                               sigma = sigma, 
                               r = r,   
                               ef = ef, 
                               einter = einter, 
                               eb = eb, 
                               erow = erow, 
                               ecol = ecol, 
                               nrand = nrand,    
                               rd = rd,   
                               randomized = randomized)
           class(result) <- c('gerexp.fe','gerexp','list')   
         },
         SPE = {

           result <- gerexp.spe(mu = mu,
                                sigma = sigma,
                                sigmap = sigmap,
                                r = r,   
                                ef = ef, 
                                einter = einter, 
                                eb = eb, 
                                erow = erow, 
                                ecol = ecol, 
                                nrand = nrand,    
                                rd = rd,   
                                randomized = randomized) 
           class(result) <- c('gerexp.spe','gerexp','list') 
         }
  )

  cl <- match.call()

  result$call <- cl

  return(result)

}
