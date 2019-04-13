plot.gexp.fe <- function(x,
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
  if(is.null(getCall(x)$blke) & is.null(getCall(x)$rowe) & is.null(getCall(x)$cole)){ #é um DIC

    aux <- update(x, random=TRUE) 
    aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
    aux2 <- aux1[, -dim(aux1)[2]]
    auxinter <- paste('interaction(',
                      paste(names(aux2),
                            collapse=','),
                      ')',
                      sep='')
    aux2$inter <- with(aux2,
                       eval(parse(text=auxinter)))

    if(is.null(main)){
      main = 'Factorial Structure \n Completely Random Design'
    }

    Lfactors <- aux2$inter

    if(is.null(sub)){

      factors <- names(aux2)[-dim(aux2)[2]]
      levelss <- paste(levels(aux2$inter),
                       collapse=',')
      repp <- length(unique(aux1$r))

      sub <- paste('Factors:',
                   paste(factors,
                         collapse=','),
                   '\n',
                   paste('Levels:',
                         levelss,
                         sep=''),
                   '\n',
                   paste('Replication:',
                         repp,
                         sep=''))
    }

    rowsquare <- length(levels(aux2$inter))
    columsquare <- length(unique(aux1$r))  

    aux_posxcentro <- 1/columsquare
    aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
    posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/columsquare)

    aux_posycentro <- 1/rowsquare
    aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
    posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

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

      text(x = rep(posxcentro, rep(length(posycentro), length(posxcentro))),
           y = rep(posycentro, length(posxcentro)),
           Lfactors,
           col = coltext)
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
      text(x = locator(),
           y = NULL,
           Lfactors,
           col = coltext)
    }         
  } else if(!is.null(getCall(x)$blke) & is.null(getCall(x)$rowe) & is.null(getCall(x)$cole)){#é um DBC

    aux <- update(x, random=TRUE) 
    auxx <- aux$dfm
    dimen <- dim(auxx)

    aux1 <- auxx[, -dimen[2]]

    labelblock <- names(aux1)[dimen[2]-1]
    nblocks <- length(levels(aux1[[labelblock]]))

    aux11 <- aux1[order(aux1[[labelblock]]), ] 

    aux2 <- aux11[, -c((dimen[2]-1):(dimen[2]-2))]
    auxinter <- paste('interaction(',
                      paste(names(aux2),
                            collapse=','),
                      ')',
                      sep='')
    aux2$inter <- with(aux2,
                       eval(parse(text=auxinter)))

    if (is.null(eval(getCall(x)$blkl))) {
      blocks <- paste(labelblock, 1:nblocks, sep = " ")
    } else {
      blocks <- levels(aux1[[labelblock]])
    }

    if(is.null(main)){
      main = 'Factorial Structure \n Random Completely Block Design'
    }

    Lfactors <- aux2$inter

    if(is.null(sub)){

      factors <- names(aux2)[-dim(aux2)[2]]
      levelss <- paste(levels(aux2$inter),
                       collapse=',')
      repp <- length(unique(aux1$r))

      sub <- paste('Factors:',
                   paste(factors,
                         collapse=','),
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

    rowsquare <- nblocks
    columsquare <- dim(aux11)[1]/nblocks

    aux_posxcentro <- 1/columsquare
    aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
    posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/columsquare)

    aux_posycentro <- 1/rowsquare
    aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
    posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

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

      text(x = rep(posxcentro, length(posycentro)),
           y = rep(posycentro, rep(length(posxcentro), length(posycentro))),
           Lfactors,
           col = coltext,
           srt = 40)

      arrows(-0.05,
             seq(0, 2, by=2/rowsquare),
             -0.05,
             seq(2/rowsquare, 2, by=2/rowsquare),
             angle=90,
             xpd=TRUE,
             code=3,
             length=0.06) 

      text(-0.08,
           posycentro,
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

      text(x = locator(),
           y = NULL,
           paste(labelblock, 1:rowsquare), 
           col = coltext)  

      text(x = locator(),
           y = NULL,
           Lfactors,
           col = coltext)
    }    
  } else {#is a LSD
    aux <- update(x, random=TRUE) 
    aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
    labelrow <- names(aux1)[1] 
    labelcol <- names(aux1)[2]
    nrows <- ncols <- length(levels(aux1[[labelrow]])) 

    aux2 <- aux1[order(aux1[[labelrow]],
                       aux1[[labelcol]]), ]

    auxinter <- paste('interaction(',
                      paste(names(aux2)[-c(1:2)],
                            collapse=','),
                      ')',
                      sep='')
    aux2$inter <- with(aux2,
                       eval(parse(text=auxinter)))

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

    if(is.null(main)){
      main = 'Factorial Structure: Latin Square Design'
    }

    Lfactors <- aux2$inter

    if(is.null(sub)){

      factors <- paste(names(aux1)[-c(1:2)],
                       collapse=',')
      levelss <- paste(levels(aux2$inter),
                       collapse=',')

      sub <- paste('Factors:',
                   factors,
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
    rowsquare <- columsquare <- length(eval(getCall(x)$rowe))

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
           Lfactors,
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
      text(x = locator(),
           y = NULL,
           paste(labelrow, 1:rowsquare), 
           col = coltext)
      text(x = locator(),
           y = NULL,
           paste(labelcol, 1:columsquare), 
           col = coltext) 
      text(x = locator(),
           y = NULL,
           Lfactors,
           col = coltext) 
    }
  }
}
