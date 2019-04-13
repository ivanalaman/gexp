plot.gexp.lsd <- function(x,
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
  labelrow <- names(aux1)[1] 
  labelcol <- names(aux1)[2]
  nrows <- ncols <- length(levels(aux1[[labelrow]]))

  aux11 <- aux1[order(aux1[[labelrow]],
                      aux1[[labelcol]]), ]

  factors <- aux11[, 3]

  if(length(attr(aux$X, 'contrasts')[-c(1:2)]) != 1){
    stop('Graphic option only for one factor!')
  }

  if(is.null(main)){
    main = 'Latin Square Design'
  }

  if(is.null(eval(getCall(x)$rowl))){
    rows <- paste(labelrow, 1:nrows)
  } else {
    rows <- levels(aux1[[labelrow]])
  }

  if(is.null(eval(getCall(x)$coll))){
    cols <- paste(labelcol, 1:ncols)
  } else {
    cols <- levels(aux1[[labelcol]])
  }

  if(is.null(sub)){

    Lfactors <- names(aux1)[3]
    levelss <- paste(levels(x$dfm[[Lfactors]]),
                     collapse=',')

    sub <- paste('Factors:',
                 Lfactors,
                 '\n',
                 paste('Levels:',
                       levelss,
                       sep=''),
                 '\n',
                 paste('Rows:',
                       nrows,
                       sep=''),
                 '\n',
                 paste('Columns:',
                       ncols,
                       sep=''))

  }

  rowsquare <- columsquare <- nrows

  aux_posxcentro <- 1/rowsquare
  aux_posxcentro1 <- aux_posxcentro + ((rowsquare - 1)*2/rowsquare)
  posxcentro <- posycentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/rowsquare)

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

    arrows(seq(0, 2, by=2/rowsquare),
           2.05,
           seq(2/rowsquare, 2, by=2/rowsquare),
           2.05,
           angle=90,
           xpd=TRUE,
           code=3,
           length=0.06) 

    text(-0.08,
         posxcentro,
         rows,
         col=colgrid,
         xpd=TRUE,
         srt=90)

    text(posxcentro,
         2.08,
         cols,
         col=colgrid,
         xpd=TRUE)
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

    tcltk::tkmessageBox(message='Click with the left button on row block and end with the right button!') 

    text(x = locator(),
         y = NULL,
         paste(labelrow, 1:rowsquare), 
         col = coltext)

    tcltk::tkmessageBox(message='Click with the left button on column block and end with the right button!') 

    text(x = locator(),
         y = NULL,
         paste(labelcol, 1:columsquare), 
         col = coltext) 

    tcltk::tkmessageBox(message='Now, click with the left button on experimental unit and end with the right button!') 

    text(x = locator(),
         y = NULL,
         factors,
         col = coltext)
  }
}
