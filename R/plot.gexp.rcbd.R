plot.gexp.rcbd <- function(x,
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
                           ...)
{

  aux <- update(x, random=TRUE) 
  aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
  labelblock <- names(aux1)[3] 
  nblocks <- length(levels(aux1[[labelblock]]))
 
  aux11 <- aux1[order(aux1[[labelblock]]), ]

  factors <- aux11[, 1]
  
 if (is.null(eval(getCall(x)$blkl))) {
    blocks <- paste(labelblock, 1:nblocks, sep = " ")
} else {
    blocks <- levels(aux1[[labelblock]])
}

  if(length(attr(aux$X, 'contrasts')[-1]) != 1){
    stop('Graphic option only for one factor!')
  }

  if(is.null(main)){
    main = 'Random Completely Block Design'
  }

  repp <- length(unique(aux1$r)) 

  if(is.null(sub)){

    aux_factors <- names(attr(x$X, 'contrasts'))
    Lfactors <- aux_factors[aux_factors!=labelblock]
    levelss <- paste(levels(x$dfm[[Lfactors]]),
                     collapse=',')

    sub <- paste('Factors:',
                 Lfactors,
                 '\n',
                 paste('Levels:',
                       levelss,
                       sep=''),
                 '\n',
                 paste('Replication:',
                       repp,
                       sep=''),
                 '\n',
                 paste('Block:',
                       nblocks,
                       sep=''))

  }

  aux_colsquare <- length(levels(factors))
  columsquare <- prod(aux_colsquare*repp)
  rowsquare <- nblocks

  aux_posxcentro <- 1/rowsquare
  aux_posxcentro1 <- aux_posxcentro + ((rowsquare - 1)*2/rowsquare)
  posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/rowsquare)

  aux_posycentro <- 1/columsquare
  aux_posycentro1 <- aux_posycentro + ((columsquare - 1)*2/columsquare)
  posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/columsquare)

  if(!dynamic){ 
    op <- par('xaxs', 'yaxs') # Original par('xaxs', 'yaxs')
    par(xaxs='i', yaxs='i')
    plot(1,
         type = 'n',
         xlim = c(0, 2),
         ylim = c(0, 2),
         axes = FALSE,
         xlab = '',
         ylab = '',
         main = main,
         sub  = sub,
         ...)
    box()
    grid(nx=columsquare,
         ny=rowsquare,
         col=colgrid,
         lty = ltygrid,
         lwd = lwdgrid)

    text(x = rep(posycentro, length(posxcentro)),
         y = rep(posxcentro, rep(length(posycentro), length(posxcentro))), 
         factors,
         col = coltext)

    arrows(-0.05,
           seq(0, 2, by=2/rowsquare),
           -0.05,
           seq(2/rowsquare, 2, by=2/rowsquare),
           angle=90,
           xpd=TRUE,
           code=3,
           length=0.06) 

    text(-0.08,
         posxcentro,
         blocks,
         col=colgrid,
         xpd=TRUE,
         srt=90)
    par(op) # Restoring the original par('xaxs', 'yaxs')
  } else {
    auxin <- tcltk::tk_choose.files()
    auxin1 <- gsub('[\\s\\S]*?\\.', '', auxin, perl=TRUE)
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
         sub  = sub,
         ...)

    rasterImage(myimage, 
                xleft = xleftimg, 
                ybottom = ybottomimg, 
                xright = xrightimg, 
                ytop = ytopimg) 

    tcltk::tkmessageBox(message='Click with the left button on block and end with the right button!')
     
    text(x = locator(),
         y = NULL,
         paste(labelblock, 1:rowsquare), 
         col = coltext) 

    tcltk::tkmessageBox(message='Now, click with the left button on experimental unit and end with the right button!')

    text(x = locator(),
         y = NULL,
         factors,
         col = coltext)
  }
}
