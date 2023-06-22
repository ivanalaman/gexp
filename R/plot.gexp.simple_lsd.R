plot.gexp.simple_lsd <- function(x,
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
                 'contrasts')[-c(1:2)]) != 1){
    stop('Graphic option only for one factor!')
  }

  mfactors <- attr(x$X,
                   'contrasts')

  labelfac <- names(mfactors)[-c(1:2)]

  labelrow <- names(mfactors)[1]

  labelcol <- names(mfactors)[2]

  levelsfac <- attr(mfactors[[labelfac]],
                    'dimnames')[[1]]

  levelsrow <- attr(mfactors[[labelrow]],
                    'dimnames')[[1]]

  levelscol <- attr(mfactors[[labelcol]],
                    'dimnames')[[1]]

  repp <- dim(x$X)[1]/(length(levelsrow))

  if(is.null(main)){
    main = 'Latin Square Design'
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
                 paste('Rows: ',
                       length(levelsrow),
                       sep = ''),
                 '\n',
                 paste('Columns: ',
                       length(levelscol),
                       sep = ''))
  }

  ifelse(random == FALSE,
         {
           rtreat <- latin(length(levelsfac),
                           levelss = levelsfac,
                           nrand = 0)
           treat <- c(rtreat)
         },
         {
           rtreat <- latin(length(levelsfac),
                           levelss = levelsfac)
           treat <- c(rtreat)
         })

  rowsquare <- columsquare <- length(levelsrow)

  aux_posxcentro <- 1/rowsquare

  aux_posxcentro1 <- aux_posxcentro + ((rowsquare - 1)*2/rowsquare)

  posxcentro <- posycentro <- seq(aux_posxcentro, 
                                  aux_posxcentro1, 
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

    arrows(seq(0, 
               2, 
               by = 2/rowsquare),
           2.05,
           seq(2/rowsquare,
               2, 
               by = 2/rowsquare),
           2.05,
           angle = 90,
           xpd = TRUE,
           code = 3,
           length = 0.06) 

    text(-0.08,
         posxcentro,
         levelsrow,
         col = colgrid,
         xpd = TRUE,
         srt = 90)

    text(posxcentro,
         2.08,
         levelscol,
         col = colgrid,
         xpd = TRUE)

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

    tcltk::tkmessageBox(message = 'Click with the left button on row block and end with the right button!') 

    text(x = locator(),
         y = NULL,
         paste(labelrow, 
               1:rowsquare), 
         col = coltext)

    tcltk::tkmessageBox(message = 'Click with the left button on column block and end with the right button!') 

    text(x = locator(),
         y = NULL,
         paste(labelcol, 
               1:columsquare), 
         col = coltext) 

    tcltk::tkmessageBox(message = 'Now, click with the left button on experimental unit and end with the right button!') 

    text(x = locator(),
         y = NULL,
         treat,
         col = coltext)
  }
}
