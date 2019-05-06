library(gexp)
#! Quantitative Factor(s) (QT)

#!. CRD - QT
level1 <- c(0, 10, 20, 30)
cont_crd <- contr.poly(level1)
crd <- gexp(mu=NULL,
            r=4,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30)), # b3
            fl=list(tra=ordered(level1)),
            contrasts=list(f1=cont_crd))
summary(crd)
plot(crd)

#!. RCBD - QT
level2 <- c(0, 2, 4, 6)
cont_rcbd <- contr.poly(level2)
rcbd <- gexp(mu=NULL,
             r=2,
             blke=c(1, 2, 3),
             blkl=list(Blk=c('B1', 'B2', 'B3')),
             fe=list(f1=c(10,   # Intercepto
                           2,   # b1
                           4,   # b2
                           6)), # b3
             fl=list(Dose=ordered(level2)),
             contrasts=list(f1=cont_rcbd,
                            Blk=diag(3)),
             type='RCBD')
summary(rcbd)
plot(rcbd)

#!. LSD - QT
level3 <- c(0, 10, 20, 30, 40)
cont_lsd <- contr.poly(level3)
lsd <- gexp(mu=NULL,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30,   # b3
                          40)), # b4
            fl=list(tra=ordered(level3)),
            rowe=c(1, 2, 3, 4, 5),
            rowl=list(Row=paste('r',
                                 1:5,
                                 sep='')),
            cole=c(5, 4, 3, 2, 1),
            coll=list(Col=paste('c',
                                 1:5,
                                 sep='')),
            contrasts=list(f1=cont_lsd,
                           Row=diag(5),
                           Col=diag(5)),
            type='LSD')
summary(lsd)
plot(lsd)

#!. FE - LSD - QT
level3 <- c(0, 10, 20, 30, 40)
cont_lsd <- contr.poly(level3)
#Efeito linear
fe_lsd <- gexp(mu=10,
               fe=list(f1=c(0, 3, 0, 0, 0),
                       f2=c(2, 3)), 
               inte=c(1, 15, rep(1, 8)),
               rowe=rep(1, 10),
               cole=rep(1, 10),
               fl=list(tra=ordered(level3),
                       B=paste('b',
                               1:2,
                               sep='')), 
               rowl=list(Row=paste('r',
                                   1:10,
                                   sep='')),
               coll=list(Col=paste('c',
                                   1:10,
                                   sep='')),
#               contrasts=list(f1=cont_lsd,
#                              f2=diag(2),
#                              Row=diag(10),
#                              Col=diag(10)),
               type='FE')
summary(fe_lsd)
plot(fe_lsd)

#!. SPE - LSD - QT
level1 <- c(0, 10, 20, 30)
cont_crd <- contr.poly(level1)
#Efeito quadrático
spe_lsd <- gexp(mu=10,
                fe=list(f1=c(0, 1, 3, 0),
                        f2=c(2, 3, 1)),
                fl=list(trat=ordered(level1),
                        SP=paste('sp',
                                 1:3,
                                 sep='')),
                inte=c(1, 15, rep(1, 10)),
                rowe=c(1, 1, 1, 1),
                cole=c(1, 1, 1, 1),
                rowl=list(Row=paste('r',
                                    1:4,
                                    sep='')),
                coll=list(Col=paste('c',
                                    1:4,
                                    sep='')),
                round=1,
#                contrasts=list(f1=cont_crd,
#                               f2=diag(3),
#                               Row=diag(4),
#                               Col=diag(4)),
                type='SPE')
summary(spe_lsd)
plot(spe_lsd)
