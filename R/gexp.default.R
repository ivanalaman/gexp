gexp.default <- function(mu        = 26, 
                         err       = NULL,   
                         errp      = NULL,   
                         r         = 5L,  
                         fl        = NULL,
                         blkl      = NULL,
                         rowl      = NULL,
                         coll      = NULL,
                         fe        = NULL,
                         inte      = NULL,
                         blke      = NULL,
                         rowe      = NULL,
                         cole      = NULL,
                         contrasts = NULL,
                         type      = c('CRD','RCBD','LSD','FE','SPE'),
                         nrand     = 1L,    
                         round     = 2L,   
                         random    = FALSE,  
                         ...)
{

  #cl <- match.call()

  switch(match.arg(type),
         CRD = {
           if(any(!is.null(inte) | !is.null(blke) | !is.null(rowe) | !is.null(cole))) stop("Only fe must be specified as effect of treatment!")

           if(is.null(fe)) fe <- list(f1 = rep(1,3)) 

           result <- gexp.crd(mu        = mu, 
                              err       = err,   
                              r         = r,   
                              fl        = fl,
                              fe        = fe,
                              contrasts = contrasts, 
                              round     = round, 
                              random    = random,
                              ...)
           class(result) <- c('gexp.crd', 'gexp', 'list')
         },
         RCBD= {
           if(is.null(blke))  blke <- rep(1,3) 

           result <- gexp.rcbd(mu        = mu,
                               err       = err,
                               r         = r,
                               fl        = fl,
                               blkl      = blkl,
                               fe        = fe,
                               blke      = blke,
                               contrasts = contrasts,
                               round     = round,
                               random    = random,
                               ...)
           class(result) <- c('gexp.rcbd', 'gexp', 'list')
         },
         LSD = {

           if(is.null(fe)) fe <- list(f1 = rep(1,3))
           if(length(fe)!=1) stop("Use only one factor this design!")  
           if(is.null(rowe)) rowe <- unlist(fe)  
           if(is.null(cole)) cole <- rowe 

           result <- gexp.lsd(mu        = mu, 
                              err       = err,
                              fl        = fl,
                              rowl      = rowl,
                              coll      = coll,
                              fe        = fe,
                              rowe      = rowe,
                              cole      = cole,
                              nrand     = nrand,
                              contrasts = contrasts,
                              round     = round,
                              random    = random,
                              ...)
           class(result) <- c('gexp.lsd', 'gexp', 'list')
         },
         FE  = {

           if(is.null(fe)) fe <- list(f1 = rep(1,2),
                                      f2 = rep(1,3))
           if(is.null(inte)) inte <- rep(1,
                                         do.call('prod',
                                                 lapply(fe,length)))

           result <- gexp.fe(mu        = mu,
                             err       = err, 
                             r         = r,
                             fl        = fl,
                             blkl      = blkl,
                             rowl      = rowl,
                             coll      = coll,
                             fe        = fe,
                             inte      = inte,
                             blke      = blke,
                             rowe      = rowe,
                             cole      = cole,
                             nrand     = nrand, 
                             contrasts = contrasts,
                             round     = round, 
                             random    = random,
                             ...)
           class(result) <- c('gexp.fe', 'gexp', 'list')   
         },
         SPE = {

           if(is.null(fe)) fe <- list(f1 = rep(1,2),
                                      f2 = rep(1,3))
           if(is.null(inte)) inte <- rep(1,
                                         do.call('prod',
                                                 lapply(fe,length)))

           result <- gexp.spe(mu        = mu,
                              err       = err,
                              errp      = errp,
                              r         = r,   
                              fl        = fl,
                              blkl      = blkl,
                              rowl      = rowl,
                              coll      = coll,
                              fe        = fe,
                              inte      = inte,
                              blke      = blke,
                              rowe      = rowe,
                              cole      = cole,
                              nrand     = nrand,
                              contrasts = contrasts,
                              round     = round,
                              random    = random,
                              ...) 
           class(result) <- c('gexp.spe', 'gexp', 'list') 
         }
  )
  result$call <- match.call()
  return(result)
}
