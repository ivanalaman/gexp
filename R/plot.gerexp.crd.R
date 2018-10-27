plot.gerexp.crd <- function(x,
                            newlevels = NULL,
                            colgrid = 'red',
                            coltext = 'blue',
                            ltygrid = 'dotted',
                            lwdgrid = par('lwd'),
                            xleftimg = par()$usr[1],
                            ybottomimg = par()$usr[3],
                            xrightimg = par()$usr[2],
                            ytopimg = par()$usr[4],
                            angleimg = 0,
                            dynamic = FALSE,
                            ...)
{

  aux <- update(x,randomized=TRUE) 
  aux1 <- aux$dfm[,-dim(aux$dfm)[2]]
  aux2 <- aux1[,-dim(aux1)[2]]

  if(!is.null(newlevels)){  
    levels(aux2) <- newlevels
  }

  factors <- aux2

  aux_rowsquare <- eval(getCall(x)$ef)
  aux_rowsquare1 <- lapply(aux_rowsquare,length)
  rowsquare <- do.call('prod',aux_rowsquare1)
  columsquare <- eval(getCall(x)$r)

  aux_posxcentro <- 1/columsquare
  aux_posxcentro1 <- aux_posxcentro + ((columsquare - 1)*2/columsquare)
  posxcentro <- seq(aux_posxcentro,aux_posxcentro1,by=2/columsquare)

  aux_posycentro <- 1/rowsquare
  aux_posycentro1 <- aux_posycentro + ((rowsquare - 1)*2/rowsquare)
  posycentro <- seq(aux_posycentro,aux_posycentro1,by=2/rowsquare)

  if(!dynamic){ 
    par(xaxs='i',yaxs='i')
    plot(1,
         type = 'n',
         xlim = c(0,2),
         ylim = c(0,2),
         axes = FALSE,
         xlab = '',
         ylab = '',
         ...)
    box()
    grid(nx=columsquare,
         ny=rowsquare,
         col=colgrid,
         lty = ltygrid,
         lwd = lwdgrid)

    text(x = rep(posxcentro,rep(length(posycentro),length(posxcentro))),
         y = rep(posycentro,length(posxcentro)),
         factors,
         col = coltext)
  } else {

    auxin <- tcltk::tk_choose.files()
    auxin1 <- gsub('[\\s\\S]*?\\.','',auxin,perl=TRUE)
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
         ...)

    rasterImage(myimage, 
                xleft = xleftimg, 
                ybottom = ybottomimg, 
                xright = xrightimg, 
                ytop = ytopimg) 
    text(x = locator(),
         y = NULL,
         factors,
         col = coltext)
  }
}
