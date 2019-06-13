plot.gexp.spe_lsd <- function(x,
                              main       = NULL,
                              sub        = NULL,
                              colgrid    = 'red',  
                              coltext    = 'blue',
                              srttext    = 30,
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

  labelfac <- names(mfactors)[-c(1:2)]

  labelrow <- names(mfactors)[1]

  labelcol <- names(mfactors)[2]

  levelsfac <- lapply(mfactors[-c(1:2)],
                      rownames)

  levelsrow <- rownames(mfactors[[1]])

  levelscol <- rownames(mfactors[[2]])

  labelinter <- suppressWarnings(do.call('interaction',
                                         levelsfac[-1]))

  levelsinter <- levels(labelinter)

  if(is.null(main)){
    main = 'Split plot Structure: Latin Square Design'
  }

  if(is.null(sub)){
    sub <- paste('Plot: ',
                 labelfac[1],
                 '\n',
                 'Levels Plot: ',
                 paste(levelsfac[[1]],
                       collapse = ', '),
                 '\n',
                 'Subplot: ',
                 labelfac[-1],
                 '\n',
                 'Levels Subplot: ',
                 paste(unlist(levelsfac[-1]),
                       collapse = ', '),
                 '\n',
                 paste('Rows: ',
                       labelrow,
                       sep = ''),
                 '\n',
                 paste('Columns: ',
                       labelcol,
                       sep = ''))
  }

  ifelse(random == FALSE,
         {
           treat <- rep(levelsinter,
                        length(levelsrow)^2)

           Labelsplot <- latin(n = length(levelsrow),
                               levelss = levelsfac[[1]],
                               nrand = 0)
         },
         {
           lplot <- rep(list(levelsinter),
                        length(levelsrow)^2)

           rtreat <- lapply(lplot,
                            sample)

           treat <- unlist(rtreat)

           Labelsplot <- latin(n = length(levelsrow),
                               levelss = levelsfac[[1]]) 
         }) 

  rowsquare <- columsquare <- length(levelsrow)

  aux_posxcentro <- 1/rowsquare

  aux_posxcentro1 <- aux_posxcentro + ((rowsquare - 1)*2/rowsquare)

  posxcentro <- seq(aux_posxcentro, 
                    aux_posxcentro1, 
                    by = 2/rowsquare)

  aux_posycentro <- 2/rowsquare*0.9

  aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)

  posycentro <- seq(aux_posycentro, 
                    aux_posycentro1, 
                    by = 2/rowsquare)

  auxsub_posxcentro <- 1/(length(levelsinter)*length(levelsfac[[1]])) 

  auxsub_posxcentro1 <- auxsub_posxcentro + ((length(levelsinter)*length(levelsfac[[1]]) - 1)*2/(length(levelsinter)*length(levelsfac[[1]])))

  subposxcentro <- seq(auxsub_posxcentro, 
                       auxsub_posxcentro1, 
                       by = 2/(length(levelsinter)*length(levelsfac[[1]])))

  auxsub_posycentro <- 1/rowsquare

  auxsub_posycentro1 <- auxsub_posycentro + ((rowsquare - 1)*2/rowsquare)

  subposycentro <- seq(auxsub_posycentro, 
                       auxsub_posycentro1, 
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
         cex.sub = 0.8,
         ...)

    box()

    grid(nx = columsquare,
         ny = rowsquare,
         col = 1,
         lty = 1,
         lwd = 1)

    grid(nx = length(levelsinter)*columsquare,
         ny = 1,
         col = c(rep(colgrid, 
                     length(levelsinter) - 1),
                 1),
         lty = c(rep(ltygrid, 
                     length(levelsinter) - 1),
                 'solid'),
         lwd = lwdgrid)   

    text(x = rep(posxcentro, 
                 length(posycentro)),
         y = rep(posycentro, 
                 rep(length(posxcentro), 
                     length(posycentro))),
         Labelsplot,
         col = coltext)

    text(x = rep(subposxcentro, 
                 length(subposycentro)),
         y = rep(subposycentro,
                 rep(length(subposxcentro), 
                     length(subposycentro))),
         treat,
         srt = srttext,
         col = colgrid)  

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
           levelsrow, 
           col = coltext)

      tcltk::tkmessageBox(message='Click with the left button on column block and end with the right button!')     

      text(x = locator(),
           y = NULL,
           levelscol, 
           col = coltext) 

      tcltk::tkmessageBox(message='Click with the left button on plot and end with the right button!')     

      text(x = locator(),
           y = NULL,
           Labelsplot, 
           col = coltext)

      tcltk::tkmessageBox(message='Click with the left button on sub-plot and end with the right button!')     
      text(x = locator(),
           y = NULL,
           treat,
           srt = srttext,
           col = colgrid)   
    }
} 
