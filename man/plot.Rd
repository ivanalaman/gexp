%===============================================================================
% Name   : plot.gexp
% Author : I.B.Allaman
% Date   : 07/06/2019 16:00
% Version: v2
%===============================================================================
\name{plot}
\alias{plot.gexp.simple_crd}
\alias{plot.gexp.simple_rcbd}
\alias{plot.gexp.simple_lsd}
\alias{plot.gexp.fe_crd}
\alias{plot.gexp.fe_rcbd}
\alias{plot.gexp.fe_lsd} 
\alias{plot.gexp.spe_crd}
\alias{plot.gexp.spe_rcbd} 
\alias{plot.gexp.spe_lsd}

\title{
  Plot methods for \code{gexp} objects
}

\description{
  These are methods for objects of class
  \code{gexp.simple_crd} - Completely Randomized Design (CRD),
  \code{gexp.simple_rcbd} - Randomized Complete Block Design (RCBD),
  \code{gexp.simple_lsd} - Latin Squares Design (LSD),
  \code{gexp.fe_crd} - Factorial Experiment (FE) to CRD,
  \code{gexp.fe_rcbd} - FE to RCBD and
  \code{gexp.fe_lsd} - FE to LSD and
  \code{gexp.spe_crd} - Split-plot Experiment (SPE) to CRD,
  \code{gexp.spe_rcbd} - SPE to RCBD,
  \code{gexp.spe_lsd} - SPE to LSD.
  The main objective of these methods is to produce an experimental croqui with randomized treatments according with the design and type of experiment. It can be very useful in experiment planning.
}

\usage{
\method{plot}{gexp.simple_crd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.simple_rcbd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.simple_lsd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.fe_crd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.fe_rcbd}(x,
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
     random     = TRUE, \dots) 

\method{plot}{gexp.fe_lsd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.spe_crd}(x,
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
     random     = TRUE, \dots) 

\method{plot}{gexp.spe_rcbd}(x,
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
     random     = TRUE, \dots)

\method{plot}{gexp.spe_lsd}(x,
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
     random     = TRUE, \dots)                          
}

\arguments{
  \item{x}{A class \code{gexp} object.}
  \item{main}{An overall title for the \code{plot}.}
  \item{sub}{An sub title for the \code{plot}.}
  \item{coltext}{A color to the text on the \code{plot}.}
  \item{srttext}{The string rotation in degrees. See \code{srt} argument of the par function.}
  \item{colgrid}{A color to the grid on the \code{plot}.}
  \item{ltygrid}{A \code{lty} to the grid on the \code{plot}.}
  \item{lwdgrid}{A \code{lwd} to the grid on the \code{plot}.}
  \item{xleftimg}{A vector (or scalar) of left x positions.}
  \item{ybottomimg}{A vector (or scalar) of bottom y positions.}
  \item{xrightimg}{A vector (or scalar) of right x positions.}
  \item{ytopimg}{A vector (or scalar) of top y positions.}
  \item{dynamic}{A logical argument to plot experimental design using image.}
  \item{random}{It is a logical argument when the purpose is to plan experiments so that randomisation of treatments occurs in the experimental units. \code{TRUE} is the default.} 
  \item{\dots}{Further arguments (required by generic).}
}

\seealso{
  \code{\link[graphics]{plot.default}}, \code{\link[graphics]{rasterImage}}
}

\examples{
#! CRD
crd <- gexp()
plot(crd)  # Default

plot(crd,  # Changing some arguments
     main = '',
     sub = '',
     coltext = 'black',
     colgrid = 'darkred',
     ltygrid = 'solid',
     lwdgrid = 3)

#! FE - CRD
crd_fe <- gexp(type = 'FE')
plot(crd_fe)

#! SPE - CRD
split_crd <- gexp(type = 'SPE')
plot(split_crd)

#! RCBD
rcbd <- gexp(design = 'RCBD')
plot(rcbd)

#! FE - RCBD
fe_rcbd <- gexp(r = 2,
                type = 'FE', 
                design = 'RCBD')
plot(fe_rcbd)

fe_rcbd1 <- gexp(r = 1,
                 blke = c(1, 2),
                 blkl = list(Blk = c('B1', 'B2')), 
                 type = 'FE',
                 design = 'RCBD')
plot(fe_rcbd1)

#! SPE - RCBD
split_rcbd <- gexp(r = 2,
                   type = 'SPE',
                   design = 'RCBD')
plot(split_rcbd)

#! LSD
lsd <- gexp(r = 1,
            design = 'LSD')
plot(lsd)

#! FE - LSD
fe_lsd <- gexp(r = 1,
               fe = list(f1 = rep(1, 2), 
                         f2 = rep(1, 2)),
               rowe = c(1, 3, 2, 1),
               cole = c(2, 2, 1, 1),
               rowl = list(Row = paste('r',
                                       1:4,
                                       sep = '')),
               coll = list(Col = paste('c',
                                       1:4,
                                       sep = '')),
               type = 'FE',
               design = 'LSD')
plot(fe_lsd)

#! SPE - LSD
split_lsd <- gexp(r = 1,
                  type = 'SPE',
                  design = 'LSD')
plot(split_lsd)

\dontrun{
#! Using images in plannig
# CRD
# Dynamic
# Open picture 'crd.jpg' whem requested
crd_i <- gexp(r = 3,
              fe = list(f1 = c(1, 1)),
              fl = list(Treat = LETTERS[1:2]))
plot(crd_i)  # Default

# dynamic plot require 'jpeg' package!
plot(crd_i,
     dynamic = TRUE,
     xleftimg = 0.6,
     ybottomimg = .6,
     xrightimg = 1.4,
     ytopimg = 1.4)

# RCBD
rcbd_i <- gexp(fe = list(f1 = c(5, 1, 1)),
               fl = list(TR = LETTERS[1:3]),
               blke = c(1, 2, 3),
               blkl = list(BLK = paste('B',
                                       1:3,
                                       sep = '')),
               design = 'RCBD')

plot(rcbd_i,
     dynamic=TRUE)

# LSD
lsd <- gexp(r = 1,
            fe = list(f1 = c(1, 1, 10)),
            rowe = c(1, 1, 1),
            cole = c(1, 1, 1),
            rowl = list(Row = paste('r',
                                    1:3,
                                    sep = '')),
            coll = list(Col = paste('c',
                                    1:3,
                                    sep = '')),
            design = 'LSD')

#TODO: dynamic plot require 'png' package!
plot(lsd,
     dynamic=TRUE)
}
}
