%===============================================================================
% Name   : summary.gexp
% Author : I.B.Allaman
% Date   : 07/06/2019 16:00
% Version: v2
%===============================================================================
\name{summary}
\alias{summary.gexp}

\title{
  Summary \code{gexp} objects.
}

\description{
  Summary \code{gexp} objects.
}

\usage{
\method{summary}{gexp}(object,
        digits=3L, \dots)
}

\arguments{
  \item{object}{A class \code{gexp} object.}
  \item{digits}{Number of decimal digits in the results. The default is \code{3}.}
  \item{\dots}{Further arguments (required by generic).}
}

\examples{
#! CRD
#! 1 factor - CRD
crd <- gexp(mu = 1,
            r = 3,
            fe = list(f1 = c(1, 1, 5, 1, 1)),
            fl = list(Treat = LETTERS[1:5]),
            round = 2)
summary(crd)

#! FE - CRD
# 5 factors (f1..f5, level^factor): 3^1 * 2^1 * 3^1 * 2^1 * 4^1 * 2 = 360 experimental units
# 5 factors (f1..f5, level^factor): 3^2 * 2^2 * 5^1 = 180 * 2 = 360 experimental units
crd_fe <- gexp(mu = 0,
               r = 2,
               fe = list(f1 = c(1, 1, 5),
                         f2 = c(1, 1),
                         f3 = c(2, 2, 1),
                         f4 = c(1, 5),
                         f5 = c(1, 2, 3, 4, 5)),
               fl = list(A = paste('a',
                                   1:3,
                                   sep = ''),
                         B = paste('b',
                                   1:2,
                                   sep = ''),
                         C = paste('c',
                                   1:3,
                                   sep = ''),
                         D = paste('d',
                                   1:2,
                                   sep = ''),
                         E = paste('e',
                                   1:5,
                                   sep = '')),
               round = 0,
               type = 'FE')
summary(crd_fe)

#! SPE - CRD
split_crd <- gexp(mu = 30,
                  fe = list(f1 = c(1, 1),
                            f2 = c(2, 3)),
                  fl = list(P = paste('p',
                                      1:2,
                                      sep = ''),
                            SP = paste('sp',
                                       1:2,
                                       sep = '')),
                  inte = c(1, 15, 1, 1),
                  round = 1,
                  type = 'SPE')
summary(split_crd)

#! RCBD
# 1 factor, 3 blocks
rcbd <- gexp(mu = 0,
             fe = list(f1 = c(5, 1, 1)),
             fl = list(TR = LETTERS[1:3]),
             blke = c(1, 2, 3),
             blkl = list(BLK = paste('B',
                                     1:3,
                                     sep = '')),
             round = 1,
             design = 'RCBD')
summary(rcbd)

#! FE - RCBD
fe_rcbd <- gexp(mu = 30,
                fe = list(f1 = c(1, 1, 1),
                          f2 = c(2, 3)),
                blke = c(1, 3),
                round = 1,
                type = 'FE',
                design = 'RCBD')
summary(fe_rcbd)

#! SPE - RCBD
split_rcbd <- gexp(mu = 30,
                   fe = list(f1 = c(1, 1),
                             f2 = c(2, 3),
                             f3 = c(1, 1, 1)),
                   fl = list(A = paste('a',
                                       1:2,
                                       sep = ''),
                             B = paste('b',
                                       1:2,
                                       sep = ''),
                             C = paste('c',
                                       1:3,
                                       sep = '')),
                   blke = c(1, 2),
                   blkl = list(BLK = paste('B',
                                           1:2,
                                           sep = '')),
                   inte = c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                            1, 2, 1, 1,
                            1, 1, 1, 1, 1, 1,
                            1, 1, 3, 3, 3, 3),
                   round = 1,
                   type = 'SPE',
                   design = 'RCBD')
summary(split_rcbd)

#! LSD
lsd <- gexp(mu = 30,
            r = 1,
            fe = list(f1 = c(1, 1, 10)),
            rowe = c(1, 1, 1),
            cole = c(1, 1, 1),
            rowl = list(Row = paste('r',
                                    1:3,
                                    sep = '')),
            coll = list(Col = paste('c',
                                    1:3,
                                    sep = '')),
            round = 1,
            design = 'LSD')
summary(lsd)

#! FE - LSD
fe_lsd <- gexp(mu = 30,
               r = 1,
               fe = list(f1 = c(1, 1),
                         f2 = c(2, 3)),
               rowe = c(1, 3, 2, 1),
               cole = c(2, 2, 1, 1),
               rowl = list(Row = paste('r',
                                       1:4,
                                       sep = '')),
               coll = list(Col = paste('c',
                                       1:4,
                                       sep = '')),
               inte = c(1, 15, 1, 1),
               round = 1,
               type = 'FE',
               design = 'LSD')
summary(fe_lsd)

#! SPE - LSD
split_lsd <- gexp(mu = 30,
                  r = 1,
                  fe = list(f1 = c(1, 1, 2),
                            f2 = c(2, 3, 1)),
                  fl = list(P = paste('p',
                                      1:3,
                                      sep = ''),
                            SP = paste('sp',
                                       1:3,
                                       sep = '')),
                  inte = c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                  rowe = c(1, 1, 1),
                  cole = c(1, 1, 1),
                  rowl = list(Row = paste('r',
                                          1:3,
                                          sep = '')),
                  coll = list(Col = paste('c',
                                          1:3,
                                          sep = '')),
                  round = 1,
                  type = 'SPE',
                  design = 'LSD')
summary(split_lsd)
}
