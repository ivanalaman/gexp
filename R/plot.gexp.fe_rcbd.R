plot.gexp.fe_rcbd <- function(x,
                              main       = NULL,
                              sub        = NULL,
                              colgrid    = 'red',
                              coltext    = 'blue',
                              ltygrid    = 'dotted',
                              lwdgrid    = par('lwd'),
                              xleftimg   = par()$usr[1],
                              ybottomimg = par()$usr[3],
                              xrightimg  = par()$usr[2],
                              ytopimg    = par()$usr[4],
                              dynamic    = FALSE,
                              random     = TRUE,
                              ...)
{
  mfactors <- attr(x$X,
                   'contrasts')

  labelfac <- names(mfactors)[-1]

  labelblock <- names(mfactors)[1]

  levelsfac <- lapply(mfactors[-1],
                      rownames)

  levelsblock <- rownames(mfactors[[1]])

  labelinter <- suppressWarnings(do.call('interaction',
                                         levelsfac))

  levelsinter <- levels(labelinter)

  repp <- dim(x$X)[1]/(length(levelsinter)*length(levelsblock)) 

  if(is.null(main)){
    main = 'Factorial Structure \n Random Completely Block Design'
  }

  if(is.null(sub)){
    sub <- paste('Factors: ',
                 paste(labelfac,
                       collapse = ', '),
                 '\n',
                 paste('Levels: ',
                       paste(unlist(levelsfac),
                             collapse = ', '),
                       sep = ''),
                 '\n',
                 paste('Replication: ',
                       repp,
                       sep = ''),
                 '\n',
                 paste('Block: ',
                       length(levelsblock),
                       sep = '')) 
  }
  
  ifelse(random == FALSE,
         {
           treat <- rep(levelsinter,
                        repp*length(levelsblock))
         },
         {
           lblock <- rep(list(rep(levelsinter,
                                  repp)),
                         length(levelsblock))
          
           rtreat <- lapply(lblock,
                            sample)
           
           treat <- unlist(rtreat)
         })
   
  rowsquare <- length(levelsblock)
 
  columsquare <- dim(x$X)[1]/rowsquare

  aux_posxcentro <- 1/columsquare
  
  aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
  
  posxcentro <- seq(aux_posxcentro, 
                    aux_posxcentro1, 
                    by = 2/columsquare)

  aux_posycentro <- 1/rowsquare
  
  aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
  
  posycentro <- seq(aux_posycentro, 
                    aux_posycentro1, 
                    by = 2/rowsquare)

  if(!dynamic){ 
    op <- par('xaxs', 'yaxs') # Original par('xaxs', 'yaxs')
    
    par(xaxs='i', 
        yaxs='i')
    
    plot(1,
         type = 'n',
         xlim = c(0, 2),
         ylim = c(0, 2),
         axes = FALSE,
         xlab = '',
         ylab = '',
         main = main,
         sub = sub,
         ...)
   
    box()
    
    grid(nx = columsquare,
         ny = rowsquare,
         col = colgrid,
         lty = ltygrid,
         lwd = lwdgrid)

    text(x = rep(posxcentro, 
                 length(posycentro)),
         y = rep(posycentro, 
                 rep(length(posxcentro), 
                     length(posycentro))),
         treat,
         col = coltext,
         srt = 40)

    arrows(-0.05,
           seq(0, 
               2, 
               by = 2/rowsquare),
           -0.05,
           seq(2/rowsquare, 
               2, 
               by = 2/rowsquare),
           angle = 90,
           xpd = TRUE,
           code = 3,
           length = 0.06) 

    text(-0.08,
         posycentro,
         levelsblock,
         col = colgrid,
         xpd = TRUE,
         srt = 90)
   
    par(op) # Restoring the original par('xaxs', 'yaxs')
  } else {
    auxin <- tcltk::tk_choose.files()
   
    auxin1 <- gsub('[\\s\\S]*?\\.', 
                   '', 
                   auxin, 
                   perl = TRUE)
   
    auxin2 <- toupper(auxin1)

    switch(auxin2,
           PNG = {
             myimage <- png::readPNG(auxin)
           },
           JPEG = {
             myimage <- jpeg::readJPEG(auxin)
           },
           JPG = {
             myimage <- jpeg::readJPEG(auxin)
           }) 

    plot(1,
         type = 'n',
         xlab = '',
         ylab = '',
         axes = FALSE,
         main = main,
         sub = sub,
         ...)

    rasterImage(myimage, 
                xleft = xleftimg, 
                ybottom = ybottomimg, 
                xright = xrightimg, 
                ytop = ytopimg) 

    text(x = locator(),
         y = NULL,
         paste(labelblock, 
               1:rowsquare), 
         col = coltext)  

    text(x = locator(),
         y = NULL,
         treat,
         col = coltext)
  }    
} 
