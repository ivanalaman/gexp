gexp.default <- function(mu        = NULL, 
                         err       = NULL,   
                         errp      = NULL,   
                         r         = 2L,  
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
                         nrand     = 1,    
                         round     = 2L,   
                         random    = FALSE,  
                         ...)
{

  cl <- match.call()     

  switch(match.arg(type),
         CRD = {
           if(any(!is.null(inte) | !is.null(blke) | !is.null(rowe) | !is.null(cole))) stop("Only fe must be specified as effect of treatment!")

           result <- gexp.crd(mu        = mu, 
                              err       = err,   
                              r         = r,   
                              fl        = fl,
                              fe        = fe,
                              contrasts = contrasts, 
                              round     = round, 
                              random    = random)
           class(result) <- c('gexp.crd', 'gexp', 'list')
         },
         RCBD= {
           if(is.null(blke)) stop("You must specify the block effect with 'blke' argument!")

           result <- gexp.rcbd(mu        = mu,
                               err       = err,
                               r         = r,
                               fl        = fl,
                               blkl      = blkl,
                               fe        = fe,
                               blke      = blke,
                               contrasts = contrasts,
                               round     = round,
                               random    = random)
           class(result) <- c('gexp.rcbd', 'gexp', 'list')
         },
         LSD = {
           if(length(fe)!=1) stop("Use only one factor this design!")

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
                              random    = random)
           class(result) <- c('gexp.lsd', 'gexp', 'list')
         },
         FE  = {
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
                             random    = random)
           class(result) <- c('gexp.fe', 'gexp', 'list')   
         },
         SPE = {
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
                              random    = random) 
           class(result) <- c('gexp.spe', 'gexp', 'list') 
         }
  )
  result$call <- cl  
  return(result)
}
