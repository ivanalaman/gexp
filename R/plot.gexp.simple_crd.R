plot.gexp.simple_crd <- function(x,
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
  if(length(attr(x$X, 
                 'contrasts')) != 1){
    stop('Graphic option only for one factor!')
  }

  mfactors <- attr(x$X,
                   'contrasts')

  labelfac <- names(mfactors)

  levelss <- attr(mfactors[[labelfac]],
                  'dimnames')[[1]]

  repp <- dim(x$X)[1]/length(levelss)

  if(is.null(main)){
    main = 'Completely Random Design'
  }

  if(is.null(sub)){
    s_levels <- paste('Levels: ',
                      paste(levelss,
                            collapse = ', '))

    sub <- paste('Factors: ',
                 labelfac,
                 '\n',
                 s_levels,
                 '\n',
                 'Replication: ',
                 repp)
  }

  ifelse(random == FALSE,
         treat <- rep(levelss,
                      repp),
         treat <- sample(rep(levelss,
                             repp)))

  rowsquare <- length(levelss)

  columsquare <- repp

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

    par(xaxs = 'i', 
        yaxs = 'i')

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
                 rep(length(posycentro),
                     length(posxcentro))),
         y = rep(posycentro,
                 length(posxcentro)),
         treat,
         col = coltext)

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

    tcltk::tkmessageBox(message='Click with the left button on experimental unit and end with the right button!')

    text(x = locator(),
         y = NULL,
         treat,
         col = coltext)
  }
}
