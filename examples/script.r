#########       CARREGANDO AS FUNÇÕES       ###########
R.utils::sourceDirectory('../R')
library(debug)
mtrace.off()
mtrace('plot.gerexp.crd')
mtrace(gexp.spe)
######### TESTANDO AS FUNÇÕES NÃO GRÁFICAS ############

#! CRD
#! 1 factor
crd1 <- gexp(mu = 10,
             r = 3,
             ef = list(f1 = c(1, -1)),
             rd = 2)
print(crd1)
summary(crd1)

mod <- lm(Y1 ~ X1,
          data=crd1$dfm)
summary(mod)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! 5 factors
crd2 <- gexp(mu = 0,
             r = 2,
             ef = list(f1 = c(1, 1, 5),
                       f2 = c(1, 1),
                       f3 = c(2, 2, 1),
                       f4 = c(1, 5),
                       f5 = c(1, 2, 3, 4, 5)),
             rd = 0)
print(crd2)

#! CRD - Multivariated
#+ In this case it is necessary to inform the error!. The default is univariate case.

errornorm <- mvtnorm::rmvnorm(n = 9,
                              sigma=matrix(c(1, 0, 0, 1),
                                  ncol = 2))
crd3 <- gexp(mu = c(0, 2),
             error = errornorm,
             r = 3,
             ef = list(f1 = matrix(c(1, 1, 5, 1, 1, 1),
                                   ncol = 2)),
             rd = 4)
print(crd3)

#! CRD - With other contrasts
#+ Orthogonal polynomios
niveis <- c(0,10,20,30) 
contrcrd4 <- contr.poly(length(niveis))

#+ Linear
crd4 <- gexp(mu = NULL,
             r = 4,
             ef = list(f1 = c(2,10,0,0)),
             labelfactors = list(Dose = ordered(niveis)),
             contrasts = list(f1 = contrcrd4))
crd44 <- crd4$dfm
crd44$dose <- as.numeric(as.character(crd44$Dose))
summary(crd44)

plot(Y1 ~ dose,crd44)

medias  <- with(crd44,tapply(Y1,dose,mean))
plot(medias ~ names(medias))

regcrd <- lm(Y1 ~ Dose,
             crd44)
summary(regcrd)
fitted(regcrd)

x <- matrix(c(niveis,
                niveis^2,
                niveis^3),
              ncol=3)  
X <- cbind(rep(1,4),x)
theta <- coef(regcrd) # parâmetros estimados por meio da ortogonalidade a título de inferência!
Z <- cbind(rep(1,4),contrcrd4)
betas <- solve(X)%*%Z%*%theta#obtendo os betas para equação de predição
betas

#+ Quadrático
crd5 <- gexp(mu = NULL,
             r = 4,
             ef = list(f1 = c(2,1,5,0)),
             labelfactors = list(Dose=ordered(c(0,10,20,30))),
             contrasts = list(f1 = contrcrd4))
crd55 <- crd5$dfm
crd55$dose <- as.numeric(as.character(crd55$Dose))

plot(Y1 ~ dose,crd55)

medias  <- with(crd55,tapply(Y1,Dose,mean))
plot(medias ~ names(medias))

regcrd2 <- lm(Y1 ~ Dose,
              crd55)
summary(regcrd2)
fitted(regcrd2)

#+ Cúbico
crd6 <- gexp(mu = NULL,
             r = 4,
             ef = list(f1 = c(2,1,1,10)),
             labelfactors = list(Dose=ordered(c(0,10,20,30))),
             contrasts = list(f1 = contrcrd4))

crd66 <- crd6$dfm
crd66$dose <- as.numeric(as.character(crd66$Dose))

plot(Y1 ~ dose,crd66)

medias  <- with(crd66,tapply(Y1,Dose,mean))
plot(medias ~ names(medias))

regcrd3 <- lm(Y1 ~ Dose,
              crd66)
summary(regcrd3)
fitted(regcrd3)

#+ Not orthogonal polynomios
#+ Linear
contrcrd7 <- matrix(c(niveis,
                      niveis^2,
                      niveis^3),
                    ncol=3) 
crd7 <- gexp(mu = NULL,
             r = 4,
             ef = list(f1 = c(2,10,0,0)),
             labelfactors = list(Dose = ordered(c(0,10,20,30))),
             contrasts = list(f1 = contrcrd7))

crd77 <- crd7$dfm
crd77$dose <- as.numeric(as.character(crd77$Dose))

plot(Y1 ~ dose,crd77)

medias  <- with(crd77,tapply(Y1,Dose,mean))
plot(medias ~ names(medias))

regcrd4 <- lm(Y1 ~ dose + I(dose^2) + I(dose^3),
              crd77)
summary(regcrd4)
fitted(regcrd4)

#! CRD - With other errors
#+ Binomial error
errorbinom <- as.matrix(rbinom(n=15, size=5, prob=0.1))

crd8 <- gexp(mu = 20,
             error = errorbinom,
             r = 5,
             ef = list(f1 = c(1,4,1)))
summary(crd8)

mod8 <- aov(Y1 ~ X1, data = crd8$dfm)
par(mfrow=c(2,2))
plot(mod8)
shapiro.test(mod8$res)

#! RCBD
rcbd1 <- gexp(mu = 0,
              ef = list(f1 = c(5, 1, 1)),
              eb = c(2, 1),
              rd = 1,
              type = 'RCBD')
print(rcbd1)

mod <- lm(Y1 ~ Block + X1,
          data=rcbd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! RCBD without replication
rcbd2 <- gexp(mu = 0,
              r = 1,
              ef = list(f1 = c(5, 1, 1)),
              eb = c(2, 1),
              rd = 1,
              type = 'RCBD')
print(rcbd2)

mod <- lm(Y1 ~ Block + X1,
          data=rcbd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! RCBD - Multivariated
rcbd_errornorm <- mvtnorm::rmvnorm(n = 18,
                                   sigma = matrix(c(1, 0, 0, 1),
                                                  ncol = 2)) 
rcbd3 <- gexp(mu = c(0, 2),
              error = rcbd_errornorm,
              ef = list(f1 = matrix(c(1, 1, 5, 1, 1, 1),
                                    ncol = 2)),
              eb = matrix(c(2, 1, 1, 2, 1, 1),
                          ncol = 2),
              rd = 1,
              type = 'RCBD')
print(rcbd3)

#! RCBD - With other contrasts
#+ Orthogonal polynomios 
contrrcbd4 <- contr.poly(4) 

rcbd4 <- gexp(ef = list(f1 = c(1,3,0,0)),
              eb = c(1,2,3),
              r = 1,
              labelfactors = list(Dose = ordered(c(0,2,4,6))),
              labelblocks = list(bloco = c('B1','B2','B3')),
              contrasts = list(f1 = contrrcbd4,
                               bloco = diag(3)),
              type = 'RCBD')
rcbd44 <- rcbd4$dfm
rcbd44$dose <- as.numeric(as.character(rcbd44$Dose))

plot(Y1 ~ dose,rcbd44)

medias  <- with(rcbd44,tapply(Y1,Dose,mean))
plot(medias ~ names(medias))

regrcbd4 <- lm(Y1 ~ bloco + Dose,
               rcbd44)
summary(regrcbd4)
fitted(regrcbd4)

regrcbd4m <- lm(medias ~ poly(c(0,2,4,6),1))
summary(regrcbd4m)
fitted(regrcbd4m)

#! LSD 
lsd1 <- gexp(mu = 30,
             ef = list(f1 = c(1, 1, 10)),
             erow = c(1, 1, 1),
             ecol = c(1, 1, 1),
             rd = 1,
             type = 'LSD')
print(lsd1)

mod <- lm(Y1 ~ Row + Column + X1,
          data=lsd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! FE - CRD
fecrd1 <- gexp(mu = 30,
               ef = list(f1 = c(1, 1, 3),
                         f2 = c(1, 1)),
               einter = c(3, 1, 1, 1, 1, 5),
               rd = 1,
               type = 'FE')
print(fecrd1)

mod <- lm(Y1 ~ X1*X2,
          data=fecrd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)
summary(tuk)
plot(tuk)

#! FE - CRD: Tree factors
fecrd2 <- gexp(mu = 30,
               ef = list(f1 = c(1, 1, 3),
                         f2 = c(1, 1),
                         f3 = c(2, 1, 1, 1)),
               einter = rep(1,50),
               rd = 1,
               type = 'FE')
print(fecrd2)

#! FE - RCBD
fercbd1 <- gexp(mu = 30,
                ef = list(f1 = c(1, 1, 1), 
                          f2 = c(2, 3)),
                eb = c(1, 3),
                einter = c(1, 15, 1, 1, 5, 1),
                rd = 1,
                type = 'FE')
print(fercbd1)

mod <- lm(Y1 ~ Block + X1*X2,
          data=fercbd1$dfm)
anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk)

#! FE - LSD
felsd1 <- gexp(mu = 30,
               ef = list(f1 = c(1, 1), 
                         f2 = c(2, 3)),
               erow = c(1, 3, 2, 1),
               ecol = c(2, 2, 1, 1),
               einter = c(1, 15, 1, 1),
               rd = 1,
               type = 'FE')
print(felsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2,
          data=felsd1$dfm)
anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT PLOT - CRD: two levels
splitcrd1 <- gexp(mu = 30,
                  ef = list(f1 = c(1, 1), 
                            f2 = c(2, 3)),
                  einter = c(1, 15, 1, 1),
                  rd = 1,
                  type = 'SPE')
print(splitcrd1)

mod <- lm(Y1 ~ X1*X2 + X1:r,
          data=splitcrd1$dfm)  #X1:r erro(a) parcela
anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ X1*X2 + Error(X1:r),
            data=splitcrd1$dfm))
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT PLOT - CRD: tree levels
splitcrd2 <- gexp(mu = 30,
                  ef = list(f1 = c(1, 1), 
                            f2 = c(2, 3),
                            f3 = c(1,1,1)),
                  einter = c(1, 15, 1, 1,rep(1,24)),
                  rd = 1,
                  type = 'SPE')
print(splitcrd2)

#SPLIT - RCBD
splitrcbd1 <- gexp(mu = 30,
                   ef = list(f1 = c(1, 1), 
                             f2 = c(2, 3),
                             f3 = c(1, 1, 1)),
                   eb = c(1,2,3,3),
                   einter = c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                              1, 2, 1, 1,
                              1, 1, 1, 1, 1, 1,
                              1, 1, 3, 3, 3, 3),
                   rd = 1,
                   type = 'SPE')
print(splitrcbd1)

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block,
          data=splitrcbd1$dfm)  #X1:Block erro(a) parcela
anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block),
                             data=splitrcbd1$dfm))
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT - LSD
splitlsd1 <- gexp(mu = 30,
                  ef = list(f1 = c(1, 1, 2), 
                            f2 = c(2, 3, 1)),
                  einter = c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                  erow = c(1, 1, 1),
                  ecol = c(1, 1, 1),
                  rd = 1,
                  type = 'SPE')
print(splitlsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column,
          data=splitlsd1$dfm)  #X1:Row:Column erro(a) parcela
anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column),
                             data=splitlsd1$dfm))
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk) 

#! SPLIT - LSD: tree factors 
splitlsd2 <- gexp(mu = 30,
                  ef = list(f1 = c(1, 1, 2), 
                            f2 = c(2, 3, 1),
                            f3 = c(1,1)),
                  einter = rep(1,39),
                  erow = c(1, 1, 1),
                  ecol = c(1, 1, 1),
                  rd = 1,
                  type = 'SPE') 
summary(splitlsd2)

########## COM FATORES FORNECIDOS PELOS USUÁRIOS #######
crdf1 <- update(crd1,
                labelfactors = list('Species' = c('Oreochromis','Piaractus')))
print(crdf1)

crdf2 <- update(crd1,
                r = 2,
                ef = list(f1 = c(1,1,1)),
                labelfactors = list('Breed' = c('Dexter','Guernsey','Jersey')))
print(crdf2)

rcbdf1 <- update(rcbd1,
                 labelfactors = list('Variety' = c('A','B','C')))
print(rcbdf1)

lsdf1 <- update(lsd1,
                labelfactors = list('Horse' = c('English','QM','ML')),
                labellsdrows = list('Linha'=1:3),
                labellsdcols = list('Coluna'=1:3)) 
print(lsdf1)

fecrdf1 <- update(fecrd1,
                  labelfactors = list('Ration'=c('0','10','20'),
                                    'Chicken'=c('Cob','Ross')))
print(fecrdf1)

fercbdf1 <- update(fercbd1,
                   newfactors = list('Ration'=c('0','10','20'),
                                     'Chicken'=c('Cob','Ross')))
print(fercbdf1)

felsdf1 <- update(felsd1,
                  labelfactors = list('Ration' = c('0','10'),
                                    'Chicken' = c('Cob','Ross'))) 
print(felsdf1)

splitcrdf1 <- update(splitcrd1,
                     labelfactors = list('Forage'=c('A','B'),
                                       'Period' = c('I1','I2')))
print(splitcrdf1)

splitlsdf1 <- update(splitlsd1,
                     labelfactors = list('Ration'=c('0','10','20'),
                                       'Chicken'=c('Cob','Ross','Caneludo')),
                     labellsdrows = list('Linha'=1:3),
                     labellsdcols = list('Coluna'=1:3))
print(splitlsdf1)

##########   TESTANDO AS FUNÇÕES GRÁFICAS  #############
##! CRD - Only one factor
# Static
plot(crd1)
plot(crdf1)
plot(crdf2)

#- Dynamic
#+ Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
crdp1 <- update(crd1, 
                r = 3,
                ef = list(f1 = c(1,1)))
plot(crdp1,
     dynamic = TRUE)

#+ With real levels
plot(crdf1,
     dynamic = TRUE)

##! RCBD - Only one factor 
# Static
plot(rcbd1)

#- Dynamic
#+ Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
rcbd1p <- update(rcbd1, 
                 r = 1,
                 ef = list(f1 = c(1,1,1,1)),
                 eb = c(1,1,1))
plot(rcbd1p,
     dynamic = TRUE)

##! LSD - Only one factor 
# Static
plot(lsd1)
plot(lsdf1)

#- Dynamic
#+ Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
plot(lsdf1,
     dynamic = TRUE)

##! FE:CRD  
# Static
plot(fecrd1)

##! FE:RCBD 
# Static
plot(fercbd1)

##! FE:LSD 
# Static
plot(felsd1)

##! SPE:CRD
plot(splitcrd1)
plot(splitcrdf1)
plot(splitcrd2)

##! SPE: RCBD
plot(splitrcbd1)

##! SPE: LSD
plot(splitlsd1)


