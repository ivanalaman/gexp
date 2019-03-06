plot.gexp.crd <- function(x,
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
  aux1 <- aux$dfm[ , -dim(aux$dfm)[2]]
  aux2 <- aux1[ , -dim(aux1)[2]]

  if(length(attr(aux$X, 'contrasts')) != 1){
  stop('Graphic option only for one factor!')
  }

  if(is.null(main)){
    main = 'Completely Random Design'
  }

  if(is.null(sub)){
  
    factors <- names(attr(x$X, 'contrasts'))
    levelss <- paste(levels(x$dfm[[factors]]),
                     collapse=',')
    repp <- eval(getCall(x)$r)

    sub <- paste('Factors:',
                 factors,
                 '\n',
                 paste('Levels:',
                       levelss,
                       sep=''),
                 '\n',
                 paste('Replication:',
                       repp,
                       sep=''))
  }

  aux_rowsquare <- eval(getCall(x)$fe)
  aux_rowsquare1 <- lapply(aux_rowsquare, length)
  rowsquare <- do.call('prod', aux_rowsquare1)
  columsquare <- eval(getCall(x)$r)

  aux_posxcentro <- 1/columsquare
  aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
  posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/columsquare)

  aux_posycentro <- 1/rowsquare
  aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
  posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

  if(!dynamic){ 
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

    text(x = rep(posxcentro, rep(length(posycentro), length(posxcentro))),
         y = rep(posycentro, length(posxcentro)),
         aux2,
         col = coltext)
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

    tcltk::tkmessageBox(message='Click with the left button on experimental unit and end with the right button!')

    text(x = locator(),
         y = NULL,
         aux2,
         col = coltext)
  }
}
