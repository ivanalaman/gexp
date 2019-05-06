#! Quantitative Factor(s) (QT)

#!. CRD - QT
level <- c(0, 10, 20, 30)
cont_crd <- contr.poly(level)
crd <- gexp(mu=NULL,
            r=4,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30)), # b3
            fl=list(tra=ordered(level)),
            contrasts=list(f1=cont_crd))
summary(crd)
plot(crd)

#!. RCBD - QT
level <- c(0, 2, 4, 6)
cont_rcbd <- contr.poly(level)
rcbd <- gexp(mu=NULL,
             r=2,
             blke=c(1, 2, 3),
             blkl=list(Blk=c('B1', 'B2', 'B3')),
             fe=list(f1=c(10,   # Intercepto
                           2,   # b1
                           4,   # b2
                           6)), # b3
             fl=list(Dose=ordered(level)),
             contrasts=list(f1=cont_rcbd,
                            Blk=diag(3)),
             type='RCBD')
summary(rcbd)
plot(rcbd)

#!. LSD - QT
level <- c(0, 10, 20, 30, 40)
cont_lsd <- contr.poly(level)
lsd <- gexp(mu=NULL,
            fe=list(f1=c(100,   # Intercepto
                          10,   # b1
                          20,   # b2
                          30,   # b3
                          40)), # b4
            fl=list(tra=ordered(level)),
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
fe_lsd <- gexp(mu=10,
               inte=c(1, 15, 1, 1, 1, 1),
               fe=list(f1=c(1, 1, 1),
                       f2=c(2, 3)),
               fl=list(A=paste('a',
                               1:3,
                               sep=''),
                       B=paste('b',
                               1:2,
                               sep='')),
               rowe=c(1, 3, 2, 1, 1, 1),
               coll=list(Col=paste('c',
                                   1:6,
                                   sep='')),
               cole=c(2, 2, 1, 1, 1, 1),
               rowl=list(Row=paste('r',
                                   1:6,
                                   sep='')),
#               contrasts=list(f1=diag(2),
#                              f2=diag(2),
#                              Row=diag(4),
#                              Col=diag(4)),
               type='FE')
summary(fe_lsd)
plot(fe_lsd)

#!. SPE - LSD - QT
spe_lsd <- gexp(mu=10,
                fe=list(f1=c(1, 1, 2),
                        f2=c(2, 3, 1)),
                fl=list(P=paste('p',
                                1:3,
                                sep=''),
                        SP=paste('sp',
                                 1:3,
                                 sep='')),
                inte=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                rowe=c(1, 1, 1),
                cole=c(1, 1, 1),
                rowl=list(Row=paste('r',
                                    1:3,
                                    sep='')),
                coll=list(Col=paste('c',
                                    1:3,
                                    sep='')),
                round=1,
                type='SPE')
summary(spe_lsd)
plot(spe_lsd)