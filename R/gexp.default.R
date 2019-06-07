gexp.default <- function(x         = NULL,
                         mu        = 26, 
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
                         type      = c('SIMPLE', 'FE', 'SPE'),
                         design    = c('CRD', 'RCBD', 'LSD'),
                         round     = 2L,
                         ...)
{
  toe <- match.arg(type)
  des <- match.arg(design)
  option <- paste(toe, 
                  des,
                  sep = '_')                 

  qualiquanti <- checkQualiQuanti(fl)

  obj <- list(mu = mu, 
              err = err, 
              errp = errp, 
              r = r, 
              fl = fl, 
              blkl = blkl,
              rowl = rowl, 
              coll = coll, 
              fe = fe, 
              inte = inte, 
              blke = blke,
              rowe = rowe, 
              cole = cole, 
              contrasts = contrasts, 
              round = round, 
              qualiquanti = qualiquanti) 
  
  class(obj) <- tolower(option)
  
  res <- gexp(obj,
              ...)

}

checkQualiQuanti <- function(fl){
  if(is.null(fl)){
    quali <- TRUE 
    quanti <- FALSE
    posquanti <- NULL 
  }else{
    quanti <- all(lapply(fl, 
                         function(x) is.numeric(x)) == TRUE)
    quali <- all(lapply(fl, 
                        function(x) is.numeric(x)) != TRUE)
    # Se não for nem quanti nem quali é pq é um híbrido! Neste caso,
    # vamos encontrar as posições.
    posquanti <- which(unlist(lapply(fl, 
                                     is.numeric)) == TRUE)  #em qual posição estão os quanti
  }
  res <- list(quali = quali,
              quanti = quanti,
              posquanti = posquanti)
  return(res)
}
