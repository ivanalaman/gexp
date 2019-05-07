# library(R.utils)
# sourceDirectory('./R')
# library(debug)
# mtrace(gexp.fe)

#! Quantitative Factor(s) (QT)

#!. DIC - QT
set.seed(1)
dic <- gexp(mu=NULL,
            r=4,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30)), # b3
            fl=list(tra=ordered(seq(0,
                                    30,
                                    by=10))),
            contrasts=list(tra=contr.poly(4)))
summary(dic)
plot(dic)

#!. DBC - QT
set.seed(2)
dbc <- gexp(mu=NULL,
             r=2,
             blke=c(1, 2, 3),
             blkl=list(blo=c('B1', 'B2', 'B3')),
             fe=list(f1=c(10,   # Intercepto
                           2,   # b1
                           4,   # b2
                           6)), # b3
             fl=list(tra=ordered(seq(0,
                                     6,
                                     by=2))),
             contrasts=list(tra=contr.poly(4),
                            blo=diag(3)),
             type='RCBD')
summary(dbc)
plot(dbc)

#!. DQL - QT
set.seed(3)
ql <- gexp(mu=NULL,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30,   # b3
                          40)), # b4
            fl=list(tra=ordered(seq(0,
                                    40,
                                    by=10))),
            rowe=c(1, 2, 3, 4, 5),
            rowl=list(row=paste('r',
                                 1:5,
                                 sep='')),
            cole=c(5, 4, 3, 2, 1),
            coll=list(col=paste('c',
                                 1:5,
                                 sep='')),
            contrasts=list(tra=contr.poly(5),
                           row=diag(5),
                           col=diag(5)),
            type='LSD')
summary(ql)
plot(ql)

#!. FE/DQL/QT
set.seed(4)
fe_ql <- gexp(mu=NULL,
              fe=list(f1=c(2, 3),
                      f2=c(10,   # Intercepto
                            5,   # b1*
                            0,   # b2
                            0,   # b3
                            0)), # b4
#              inte=c(1, 15, rep(1, 6)),
              rowe=rep(1, 10),
              cole=rep(1, 10),
              fl=list(var=paste('v',
                                1:2,
                                sep=''),
                      tra=ordered(seq(0,
                                      40,
                                      by=10))),
              coll=list(col=paste('c',
                                  1:10,
                                  sep='')),
              rowl=list(row=paste('r',
                                  1:10,
                                  sep='')),
              contrasts=list(var=diag(2),
                             tra=contr.poly(5),
                             col=diag(10),
                             row=diag(10)),
              type='FE')

summary(fe_ql)
plot(fe_ql)

#!. PSD/QL/QT
set.seed(5)
spe_ql <- gexp(mu=NULL,
               fe=list(f1=c(2, 3, 1),
                       f2=c(100,   # Intercepto
                              1,   # b1
                              5,   # b2*
                              1)), # b3
               fl=list(p=paste('p',
                                1:3,
                                sep=''),
                       sp=ordered(seq(0,
                                      30,
                                      by=10))),
#               inte=c(1, 15, rep(1, 7)),
               rowe=c(1, 2, 3),
               cole=c(3, 2, 1),
               rowl=list(row=paste('r',
                                   1:3,
                                   sep='')),
               coll=list(col=paste('c',
                                   1:3,
                                   sep='')),
               round=1,
               contrasts=list(p=diag(3),
                              sp=contr.poly(4),
                              row=diag(3),
                              col=diag(3)),
               type='SPE')
summary(spe_ql)
plot(spe_ql)
