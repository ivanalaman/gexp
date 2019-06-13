plot.gexp.simple_rcbd <- function(x,
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
  if(length(attr(x$X, 'contrasts')[-1]) != 1){
    stop('Graphic option only for one factor!')
  }

  mfactors <- attr(x$X,
                   'contrasts')

  labelfac <- names(mfactors)[-1]

  labelblock <- names(mfactors)[1]

  levelsfac <- attr(mfactors[[labelfac]],
                    'dimnames')[[1]]

  levelsblock <- attr(mfactors[[labelblock]],
                      'dimnames')[[1]]

  repp <- dim(x$X)[1]/(length(levelsblock)*length(levelsfac))


  if(is.null(main)){
    main = 'Random Completely Block Design'
  }

  if(is.null(sub)){
    sub <- paste('Factors: ',
                 labelfac,
                 '\n',
                 paste('Levels: ',
                       paste(levelsfac,
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

  ifelse(random==FALSE,
         {
           treat <- rep(levelsfac,
                        repp*length(levelsblock))
         },
         {
           lblock <- rep(list(rep(levelsfac,
                                  repp)),
                         length(levelsblock))

           rtreat <- lapply(lblock,
                            sample)

           treat <- unlist(rtreat)
         })

  aux_colsquare <- length(levelsfac)

  columsquare <- prod(aux_colsquare*repp)

  rowsquare <- length(levelsblock)

  aux_posxcentro <- 1/rowsquare

  aux_posxcentro1 <- aux_posxcentro + ((rowsquare - 1)*2/rowsquare)

  posxcentro <- seq(aux_posxcentro, 
                    aux_posxcentro1, 
                    by = 2/rowsquare)

  aux_posycentro <- 1/columsquare

  aux_posycentro1 <- aux_posycentro + ((columsquare - 1)*2/columsquare)

  posycentro <- seq(aux_posycentro,
                    aux_posycentro1, 
                    by = 2/columsquare)

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

    text(x = rep(posycentro, 
                 length(posxcentro)),
         y = rep(posxcentro,
                 rep(length(posycentro),
                     length(posxcentro))), 
         treat,
         col = coltext)

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
         posxcentro,
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

    tcltk::tkmessageBox(message = 'Click with the left button on block and end with the right button!')

    text(x = locator(),
         y = NULL,
         paste(labelblock, 
               1:rowsquare), 
         col = coltext) 

    tcltk::tkmessageBox(message = 'Now, click with the left button on experimental unit and end with the right button!')

    text(x = locator(),
         y = NULL,
         treat,
         col = coltext)
  }
}
