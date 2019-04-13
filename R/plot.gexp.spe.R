plot.gexp.spe <- function(x,
                          main       = NULL,
                          sub        = NULL,
                          coltext    = 'blue',
                          srttext    = 30,
                          colgrid    = 'red', 
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

    aux <- update(x, random=FALSE) 
    aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
    aux2 <- aux1[, -dim(aux1)[2]]
    plott <- names(aux2)[1] 
    subplott <- names(aux2)[-1] 
    aux22 <- aux2[order(aux2[[plott]]),]

    auxx <- names(aux22)[-1]
    auxinter <- paste('interaction(',
                      paste(auxx,
                            collapse=','),
                      ')',
                      sep='')
    aux22$inter <- with(aux22,
                        eval(parse(text=auxinter)))



    if(is.null(main)){
      main = 'Split plot Structure \n Completely Random Design'
    }

    Lplot <- levels(aux2[, 1])
    nplot <- length(Lplot)
    Lsub <- levels(aux22$inter)
    nsub <- length(Lsub)
    repp <- length(unique(aux1$r)) 

    auxLplot <- sample(rep(Lplot, repp))

    auxmatrplot <- rep(list(NA), length(auxLplot))
    names(auxmatrplot) <- auxLplot

    matrplot <- lapply(auxmatrplot, function(x) sample(Lsub))


    if(is.null(sub)){

      sub <- paste('Plot:',
                   plott,
                   '\n',
                   'Levels Plot:',
                   paste(Lplot,
                         collapse=','),
                   '\n',
                   'Subplot:',
                   paste(subplott,
                         collapse=','),
                   '\n',
                   'Levels Subplot:',
                   paste(Lsub,
                         collapse=','),
                   '\n',
                   'Replication:',
                   repp)
    }

    rowsquare <- nplot
    columsquare <- repp  

    aux_posxcentro <- 1/columsquare
    aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
    posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/columsquare)

    aux_posycentro <- 2/rowsquare*0.9
    aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
    posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

    auxsub_posxcentro <- 1/(nsub*repp) 
    auxsub_posxcentro1 <- auxsub_posxcentro + ((nsub*repp - 1)*2/(nsub*repp))
    subposxcentro <- seq(auxsub_posxcentro, auxsub_posxcentro1, by=2/(nsub*repp))

    auxsub_posycentro <- 1/rowsquare
    auxsub_posycentro1 <- auxsub_posycentro + ((rowsquare - 1)*2/rowsquare)
    subposycentro <- seq(auxsub_posycentro, auxsub_posycentro1, by=2/rowsquare)


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
           col=1,
           lty = 1,
           lwd = 1) 

      grid(nx  = nsub*repp,
           ny  = 1,
           col = c(rep(colgrid, nsub-1), 1),
           lty = c(rep(ltygrid, nsub-1), 'solid'),
           lwd = lwdgrid)  

      text(x = rep(posxcentro, rep(length(posycentro), length(posxcentro))),
           y = rep(posycentro, length(posxcentro)),
           auxLplot,
           col = coltext)

      text(x = rep(subposxcentro, length(subposycentro)),
           y = rep(subposycentro, rep(length(subposxcentro), length(subposycentro))),
           unlist(matrplot),
           srt = srttext,
           col = colgrid) 
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

      tcltk::tkmessageBox(message='Click with the left button on plot and end with the right button!') 

      text(x = locator(),
           y = NULL,
           auxLplot,
           col = coltext)

      tcltk::tkmessageBox(message='Now, click with the left button on sub plot and end with the right button!') 

      text(x = locator(),
           y = NULL,
           unlist(matrplot),
           col = coltext) 
    }         
  } else if(!is.null(getCall(x)$blke) & is.null(getCall(x)$rowe) & is.null(getCall(x)$cole)){#é um DBC

    aux <- update(x, random=FALSE) 
    aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
    aux2 <- aux1[, -c(dim(aux1)[2]-1)]#tirando a repetição

    labelblock <- names(aux2)[dim(aux2)[2]]
    plott <- names(aux2)[1] 
    subplott <- names(aux2)[-c(1, dim(aux2)[2])] 

    aux22 <- aux2[order(aux2[[labelblock]], aux2[[plott]]), ]

    auxinter <- paste('interaction(',
                      paste(subplott,
                            collapse=','),
                      ')',
                      sep='')
    aux22$inter <- with(aux22,
                        eval(parse(text=auxinter)))

    nblock <- length(levels(aux2[[labelblock]]))

    if (is.null(eval(getCall(x)$blkl))) {
      blocks <- paste(labelblock, 1:nblock, sep = " ")
    } else {
      blocks <- levels(aux1[[labelblock]])
    }

    if(is.null(main)){
      main = 'Split plot Structure \n Random Completely Block Design'
    }

    Lplot <- levels(aux2[, 1])
    nplot <- length(Lplot)
    Lsub <- levels(aux22$inter)
    nsub <- length(Lsub)
    repp <- length(unique(aux1$r)) 

    auxLplot <- sample(rep(Lplot, repp))

    auxmatrblock <- rep(list(NA), nblock)
    names(auxmatrblock) <- 1:nblock

    matrblock <- lapply(auxmatrblock,
                        function(...) rep(list(NA),
                                          length(auxLplot)))

    for(i in 1:nblock){
      names(matrblock[[i]]) <- sample(auxLplot)
    }

    for(i in 1:nblock){
      for(j in 1:length(auxLplot)){
        matrblock[[i]][[j]] <- sample(Lsub)
      }
    } 

    if(is.null(sub)){

      sub <- paste('Plot:',
                   plott,
                   '\n',
                   'Levels Plot:',
                   paste(Lplot,
                         collapse=','),
                   '\n',
                   'Subplot:',
                   paste(subplott,
                         collapse=','),
                   '\n',
                   'Levels Subplot:',
                   paste(Lsub,
                         collapse=','),
                   '\n',
                   'Replication:',
                   repp,
                   '\n',
                   'Block:',
                   nblock)
    }

    rowsquare <- nblock
    columsquare <- nplot*repp 

    aux_posxcentro <- 1/columsquare
    aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
    posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/columsquare)

    aux_posycentro <- 2/rowsquare*0.9
    aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
    posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

    auxsub_posxcentro <- 1/(nsub*repp*nplot) 
    auxsub_posxcentro1 <- auxsub_posxcentro + ((nsub*repp*nplot - 1)*2/(nsub*repp*nplot))
    subposxcentro <- seq(auxsub_posxcentro, auxsub_posxcentro1, by=2/(nsub*repp*nplot))

    auxsub_posycentro <- 1/rowsquare
    auxsub_posycentro1 <- auxsub_posycentro + ((rowsquare - 1)*2/rowsquare)
    subposycentro <- seq(auxsub_posycentro, auxsub_posycentro1, by=2/rowsquare)

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
           cex.sub=0.8,
           ...)

      box()

      grid(nx=columsquare,
           ny=rowsquare,
           col=1,
           lty = 1,
           lwd = 1) 

      grid(nx  = nsub*columsquare,
           ny  = 1,
           col = c(rep(colgrid, nsub-1), 1),
           lty = c(rep(ltygrid, nsub-1), 'solid'),
           lwd = lwdgrid)  

      text(x = rep(posxcentro, length(posycentro)),
           y = rep(posycentro, rep(length(posxcentro), length(posycentro))),
           unlist(lapply(matrblock, names)),
           col = coltext)

      text(x = rep(subposxcentro, length(subposycentro)),
           y = rep(subposycentro, rep(length(subposxcentro), length(subposycentro))),
           unlist(matrblock),
           srt = srttext,
           col = colgrid)  

      arrows(-0.05,
             seq(0, 2, by=2/rowsquare),
             -0.05,
             seq(2/rowsquare, 2, by=2/rowsquare),
             angle=90,
             xpd=TRUE,
             code=3,
             length=0.06) 

      text(-0.08,
           subposycentro,
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

      tcltk::tkmessageBox(message='Click with the left button on plot and end with the right button!')   

      text(x = locator(),
           y = NULL,
           unlist(lapply(matrblock, names)), 
           col = coltext)

      tcltk::tkmessageBox(message='Click with the left button on sub plot and end with the right button!')    

      text(x = locator(),
           y = NULL,
           unlist(matrblock), 
           col = coltext)
    }    
  } else {#is a LSD
    aux <- update(x, random=TRUE) 
    aux1 <- aux$dfm[, -dim(aux$dfm)[2]]
    labelrow <- names(aux1)[1] 
    labelcol <- names(aux1)[2]
    nrows <- ncols <- length(levels(aux1[[labelrow]]))  

    aux2 <- aux1[order(aux1[[labelrow]],
                       aux1[[labelcol]]),]

    auxinter <- paste('interaction(',
                      paste(names(aux2)[-c(1:3)],
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
      main = 'Split plot Structure: Latin Square Design'
    }

    Lplot <- levels(aux2[, 3])
    nplot <- length(Lplot)
    Lsub <- levels(aux2$inter)
    nsub <- length(Lsub)

    auxmatr <- rep(list(NA), nrows)
    names(auxmatr) <- 1:nrows

    matr <- lapply(auxmatr,
                   function(...) rep(list(NA),
                                     ncols))
    for(i in 1:ncols){
      names(matr[[i]]) <- 1:ncols
    }

    plott <- names(aux2)[3]
    auxlsd <- matrix(aux2[[plott]], nrow=nrows, byrow=T)
    lsd <- apply(auxlsd, 1, unique)

    for(i in 1:nrows){
      for(j in 1:ncols){
        matr[[i]][[j]] <- list(sample(Lsub))
        names(matr[[i]][[j]]) <- lsd[i, j] 
      }
    }

    if(is.null(sub)){

      subplott <- names(aux2)[-c(1:3, dim(aux2)[2])]

      sub <- paste('Plot:',
                   plott,
                   '\n',
                   'Levels Plot:',
                   paste(Lplot,
                         collapse=','),
                   '\n',
                   'Subplot:',
                   paste(subplott,
                         collpase=','),
                   '\n',
                   'Levels Subplot:',
                   paste(Lsub,
                         collapse=','), 
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
    posxcentro <- seq(aux_posxcentro, aux_posxcentro1, by=2/rowsquare)

    aux_posycentro <- 2/rowsquare*0.9
    aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
    posycentro <- seq(aux_posycentro, aux_posycentro1, by=2/rowsquare)

    auxsub_posxcentro <- 1/(nsub*nplot) 
    auxsub_posxcentro1 <- auxsub_posxcentro + ((nsub*nplot - 1)*2/(nsub*nplot))
    subposxcentro <- seq(auxsub_posxcentro, auxsub_posxcentro1, by=2/(nsub*nplot))

    auxsub_posycentro <- 1/rowsquare
    auxsub_posycentro1 <- auxsub_posycentro + ((rowsquare - 1)*2/rowsquare)
    subposycentro <- seq(auxsub_posycentro, auxsub_posycentro1, by=2/rowsquare)

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
           cex.sub = 0.8,
           ...)
      box()

      grid(nx=columsquare,
           ny=rowsquare,
           col=1,
           lty = 1,
           lwd = 1)

      grid(nx  = nsub*columsquare,
           ny  = 1,
           col = c(rep(colgrid, nsub-1), 1),
           lty = c(rep(ltygrid, nsub-1), 'solid'),
           lwd = lwdgrid)   

      text(x = rep(posxcentro, length(posycentro)),
           y = rep(posycentro, rep(length(posxcentro), length(posycentro))),
           unlist(lapply(matr, function(x)lapply(x, names))),
           col = coltext)

      text(x = rep(subposxcentro, length(subposycentro)),
           y = rep(subposycentro, rep(length(subposxcentro), length(subposycentro))),
           unlist(matr),
           srt = srttext,
           col = colgrid)  

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
    }else{
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

      tcltk::tkmessageBox(message='Click with the left button on plot and end with the right button!')     

      text(x = locator(),
           y = NULL,
           unlist(lapply(matr, function(x)lapply(x, names))), 
           col = coltext)

      tcltk::tkmessageBox(message='Click with the left button on sub plot and end with the right button!')     

      text(x = locator(),
           y = NULL,
           unlist(matr),
           srt = srttext,
           col = colgrid)   
    }
  }
}
