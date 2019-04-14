%===============================================================================
% Name   : plot.gexp
% Author : J.C.Faria
% Date   : 17/02/2019 20:41:53
% Version: v1
%===============================================================================
\name{plot}
\alias{plot.gexp.crd}
\alias{plot.gexp.rcbd}
\alias{plot.gexp.lsd}
\alias{plot.gexp.fe}
\alias{plot.gexp.spe}

\title{
  Plot methods for \code{gexp} objects
}

\description{
  These are methods for objects of class
  \code{gexp.crd} - Completely Randomized Design (CRD),
  \code{gexp.rcbd} - Randomized Complete Block Design (RCBD),
  \code{gexp.lsd} - Latin Squares Design (LSD),
  \code{gexp.fe} - Factorial Experiment (FE) and
  \code{gexp.spe} - Split-plot Experiment (SPE).
  The main objective of these methods is to produce an experimental croqui with randomized treatments according with the design and type of experiment. It can be very useful in experiment planning.
}

\usage{
\method{plot}{gexp.crd}(x,
     main= NULL,
     sub= NULL,
     colgrid='red',
     coltext='blue',
     ltygrid='dotted',
     lwdgrid=par('lwd'),
     xleftimg=par()$usr[1],
     ybottomimg=par()$usr[3],
     xrightimg=par()$usr[2],
     ytopimg=par()$usr[4],
     dynamic=FALSE, \dots)

\method{plot}{gexp.rcbd}(x,
     main=NULL,
     sub=NULL,
     colgrid='red',
     coltext='blue',
     ltygrid='dotted',
     lwdgrid=par('lwd'),
     xleftimg=par()$usr[1],
     ybottomimg=par()$usr[3],
     xrightimg=par()$usr[2],
     ytopimg=par()$usr[4],
     dynamic=FALSE, \dots)

\method{plot}{gexp.lsd}(x,
     main=NULL,
     sub=NULL,
     colgrid='red',
     coltext='blue',
     ltygrid='dotted',
     lwdgrid=par('lwd'),
     xleftimg=par()$usr[1],
     ybottomimg=par()$usr[3],
     xrightimg=par()$usr[2],
     ytopimg=par()$usr[4],
     dynamic=FALSE, \dots)

\method{plot}{gexp.fe}(x,
     main=NULL,
     sub=NULL,
     colgrid='red',
     coltext='blue',
     ltygrid='dotted',
     lwdgrid=par('lwd'),
     xleftimg=par()$usr[1],
     ybottomimg=par()$usr[3],
     xrightimg=par()$usr[2],
     ytopimg=par()$usr[4],
     dynamic=FALSE, \dots)

\method{plot}{gexp.spe}(x,
     main=NULL,
     sub=NULL,
     coltext='blue',
     srttext=30,
     colgrid='red',
     ltygrid='dotted',
     lwdgrid=par('lwd'),
     xleftimg=par()$usr[1],
     ybottomimg=par()$usr[3],
     xrightimg=par()$usr[2],
     ytopimg=par()$usr[4],
     dynamic=FALSE, \dots)
}

\arguments{
  \item{x}{A \code{gexp.crd}, \code{gexp.rcbd}, \code{gexp.lsd}, \code{gexp.fe} or \code{gexp.spe} object.}
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
  \item{\dots}{Further arguments (required by generic).}
}

\seealso{
   \code{\link[graphics]{plot.default}}, \code{\link[graphics]{rasterImage}}
}

\examples{
#! Completely Randomized Design (CRD)
#! 1 factor - CRD
crd <- gexp(mu=1,
            r=3,
            fe=list(f1=c(1, 1, 5, 1, 1)),
            fl=list(Treat=LETTERS[1:5]),
            round=2)
summary(crd)
plot(crd)  # Default

plot(crd,  # Changing some arguments
     main='',
     sub='',
     coltext='black',
     colgrid='darkred',
     ltygrid='solid',
     lwdgrid=3)

#! Factorial Experiment (FE) - CRD
# 2 factors (f1,f2, level^factor): 3^1 * 2^1 = 6 experimental units
crd_fe <- gexp(mu=0,
               r=2,
               fe=list(f1=c(1, 1, 5),
                       f2=c(1, 1)),
               fl=list(A=paste('a',
                               1:3,
                               sep=''),
                       B=paste('b',
                               1:2,
                               sep='')),
               round=0,
               inte = rep(1, 6),
               type = 'FE')
summary(crd_fe)
plot(crd_fe)

#! Split-plot Experiment (SPE) - CRD
split_crd <- gexp(mu=30,
                  fe=list(f1=c(1, 1),
                          f2=c(2, 3)),
                  fl=list(P=paste('p',
                                  1:2,
                                  sep=''),
                          SP=paste('sp',
                                   1:2,
                                   sep='')),
                  inte=c(1, 15, 1, 1),
                  round=1,
                  type='SPE')
summary(split_crd)
plot(split_crd)

#! Randomized Complete Block Design (RCBD)
# 1 factor, 3 blocks
rcbd <- gexp(mu=0,
             fe=list(f1=c(5, 1, 1)),
             fl=list(TR=LETTERS[1:3]),
             blke=c(1, 2, 3),
             blkl=list(BLK=paste('B',
                                 1:3,
                                 sep='')),
             round=1,
             type='RCBD')
summary(rcbd)
plot(rcbd)

#! Factorial Experiment (FE) - RCBD
fe_rcbd <- gexp(mu=30,
                fe=list(f1=c(1, 1, 1),
                        f2=c(2, 3)),
                blke=c(1, 3),
                inte=c(1, 15, 1, 1, 5, 1),
                round=1,
                type='FE')
summary(fe_rcbd)
plot(fe_rcbd)

fe_rcbd1 <- gexp(mu=30,
                 fe=list(f1=c(1, 1, 1),
                         f2=c(2, 3)),
                 blke=c(1, 3),
                 blkl=list(Blk=c('B1', 'B2')), 
                 inte=c(1, 15, 1, 1, 5, 1),
                 round=1,
                 type='FE')
summary(fe_rcbd1)
plot(fe_rcbd1)

#! Split-plot Experiment (SPE) - RCBD
split_rcbd <- gexp(mu=30,
                   r = 2,
                   fe=list(f1=c(1, 1),
                           f2=c(2, 3),
                           f3=c(1, 1, 1)),
                   fl=list(A=paste('a',
                                   1:2,
                                   sep=''),
                           B=paste('b',
                                   1:2,
                                   sep=''),
                           C=paste('c',
                                   1:3,
                                   sep='')),
                   blke=c(1, 2, 3),
                   blkl=list(BLK=paste('B',
                                       1:3,
                                       sep='')),
                   inte=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                          1, 2, 1, 1,
                          1, 1, 1, 1, 1, 1,
                          1, 1, 3, 3, 3, 3),
                   round=1,
                   type='SPE')
summary(split_rcbd)
plot(split_rcbd)

#! Latin Square Design (LSD)
lsd <- gexp(mu=30,
            fe=list(f1=c(1, 1, 10)),
            rowe=c(1, 1, 1),
            cole=c(1, 1, 1),
            rowl=list(Row=paste('r',
                                1:3,
                                sep='')),
            coll=list(Col=paste('c',
                                1:3,
                                sep='')),
            round=1,
            type='LSD')
summary(lsd)
plot(lsd)

#! Factorial Experiment (FE) - LSD
fe_lsd <- gexp(mu=30,
               fe=list(f1=c(1, 1),
                       f2=c(2, 3)),
               rowe=c(1, 3, 2, 1),
               cole=c(2, 2, 1, 1),
               rowl=list(Row=paste('r',
                                   1:4,
                                   sep='')),
               coll=list(Col=paste('c',
                                   1:4,
                                   sep='')),
               inte=c(1, 15, 1, 1),
               round=1,
               type='FE')
summary(fe_lsd)
plot(fe_lsd)

#! Split-plot Experiment (SPE) - LSD
split_lsd <- gexp(mu=30,
                  fe=list(f1=c(1, 1, 2),
                          f2=c(2, 3, 1)),
                  fl=list(P=paste('p',
                                  1:3,
                                  sep=''),
                          SP=paste('sp',
                                   1:3,
                                   sep='')),
                  inte=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                  rowe = c(1, 1, 1),
                  cole = c(1, 1, 1),
                  rowl=list(Row=paste('r',
                                      1:3,
                                      sep='')),
                  coll=list(Col=paste('c',
                                      1:3,
                                      sep='')),
                  round=1,
                  type='SPE')
summary(split_lsd)
plot(split_lsd)

\dontrun{
#! Using images in plannig
# CRD
crd <- gexp(mu=1,
            r=2,
            fe=list(f1=c(1, 1, 5, 1, 1)),
            fl=list(Treat=LETTERS[1:5]),
            round=2)
summary(crd)
plot(crd)  # Default

# Dynamic
# Open picture 'crd.jpg' whem requested
crd_i <- update(crd,
                r=3,
                fe=list(f1=c(1, 1)),
                fl=list(Treat=LETTERS[1:2]))
summary(crd_i)
plot(crd_i)  # Default

# dynamic plot require 'jpeg' package!
plot(crd_i,
     dynamic=TRUE,
     xleftimg=0.6,
     ybottomimg=.6,
     xrightimg=1.4,
     ytopimg=1.4)

# RCBD
rcbd <- gexp(mu=0,
             fe=list(f1=c(5, 1, 1)),
             fl=list(TR=LETTERS[1:3]),
             blke=c(1, 2, 3),
             blkl=list(BLK=paste('B',
                                 1:3,
                                 sep='')),
             round=1,
             type='RCBD')
summary(rcbd)
plot(rcbd)

rcbd_i <- update(rcbd,
                 r=1,
                 fe=list(f1=c(1, 1, 1, 1)),
                 blke=c(1, 1, 1))
plot(rcbd_i,
     dynamic=TRUE)

# LSD
lsd <- gexp(mu=30,
            fe=list(f1=c(1, 1, 10)),
            rowe=c(1, 1, 1),
            cole=c(1, 1, 1),
            rowl=list(Row=paste('r',
                                1:3,
                                sep='')),
            coll=list(Col=paste('c',
                                1:3,
                                sep='')),
            round=1,
            type='LSD')
summary(lsd)

#TODO: dynamic plot require 'png' package!
plot(lsd)
plot(lsd,
     dynamic=TRUE)
}
}

