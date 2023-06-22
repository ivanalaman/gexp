plot.gexp.spe_rcbd <- function(x,
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

  labelblock <- names(mfactors)[1]

  labelfac <- names(mfactors)[-1]

  levelsblock <- rownames(mfactors[[1]])

  levelsfac <- lapply(mfactors[-1],
                      rownames)

  labelinter <- suppressWarnings(do.call('interaction',
                                         levelsfac[-1]))

  levelsinter <- levels(labelinter)

  repp <- dim(x$X)[1]/(length(levelsinter)*length(levelsfac[[1]])*length(levelsblock)) 

  if(is.null(main)){
    main = 'Split plot Structure \n Random Completely Block Design'
  }

  if(is.null(sub)){
    sub <- paste('Plot:',
                 labelfac[1],
                 '\n',
                 'Levels Plot: ',
                 paste(levelsfac[[1]],
                       collapse = ', '),
                 '\n',
                 'Subplot: ',
                 paste(labelfac[-1],
                       collapse = ', '),
                 '\n',
                 'Levels Subplot: ',
                 paste(paste(unlist(levelsfac[-1]),
                             collapse = ', '),
                       collapse = ', '),
                 '\n',
                 'Replication: ',
                 repp,
                 '\n',
                 'Block: ',
                 length(levelsblock))
  }

  ifelse(random == FALSE,
         {
           treat <- rep(levelsinter,
                        repp*length(levelsblock)*length(levelsfac[[1]]))

           Labelsplot <- rep(levelsfac[[1]],
                             repp)
         },
         {
           lplot <- rep(list(levelsinter),
                        length(levelsfac[[1]])*repp*length(levelsblock))

           rtreat <- lapply(lplot,
                            sample)

           treat <- unlist(rtreat)

           lblock <- rep(list(rep(levelsfac[[1]],
                                  repp)),
                         length(levelsblock))

           rabelsplot <- lapply(lblock,
                                sample)

           Labelsplot <- unlist(rabelsplot)
         }) 

  rowsquare <- length(levelsblock)

  columsquare <- length(levelsfac[[1]])*repp 

  aux_posxcentro <- 1/columsquare

  aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)

  posxcentro <- seq(aux_posxcentro, 
                    aux_posxcentro1, 
                    by = 2/columsquare)

  aux_posycentro <- 2/rowsquare*0.9

  aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)

  posycentro <- seq(aux_posycentro, 
                    aux_posycentro1, 
                    by = 2/rowsquare)                                 

  auxsub_posxcentro <- 1/(length(levelsinter)*repp*length(levelsfac[[1]])) 

  auxsub_posxcentro1 <- auxsub_posxcentro + ((length(levelsinter)*repp*length(levelsfac[[1]]) - 1)*2/(length(levelsinter)*repp*length(levelsfac[[1]])))

  subposxcentro <- seq(auxsub_posxcentro, 
                       auxsub_posxcentro1, 
                       by = 2/(length(levelsinter)*repp*length(levelsfac[[1]])))

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
                     length(levelsinter)-1), 
                 1),
         lty = c(rep(ltygrid, 
                     length(levelsinter)-1), 
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

    text(-0.08,
         subposycentro,
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

    tcltk::tkmessageBox(message = 'Click with the left button on plot and end with the right button!')   

    text(x = locator(),
         y = NULL,
         Labelsplot, 
         col = coltext)

    tcltk::tkmessageBox(message = 'Click with the left button on sub plot and end with the right button!')    

    text(x = locator(),
         y = NULL,
         treat, 
         col = coltext)
  }   
} 
