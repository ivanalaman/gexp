#########       CARREGANDO AS FUNÇÕES       ###########
sourceDirectory('../R')

######### TESTANDO AS FUNÇÕES NÃO GRÁFICAS ############

#! CRD
#! 1 factor
crd1 <- gerexp(mu=0,
               sigma=1,
               r=3,
               ef=list(f1=c(1, 1, 5, 1, 1)),
               rd=2)
print(crd1)
summary(crd1)

mod <- lm(Y1 ~ X1,
          data=crd1$dfm)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! 5 factors
crd2 <- gerexp(mu=0,
               sigma=1,
               r=2,
               ef=list(f1=c(1, 1, 5),
                       f2=c(1, 1),
                       f3=c(2, 2, 1),
                       f4=c(1, 5),
                       f5=c(1, 2, 3, 4, 5)),
               rd=0)
print(crd2)

#! CRD - Multivariated
crd3 <- gerexp(mu=c(0, 2),
               sigma=matrix(c(1, 0, 0, 1),
                            ncol=2),
               r=3,
               ef=list(f1=matrix(c(1, 1, 5, 1, 1, 1),
                                 ncol=2),
                       f2=matrix(c(1, 3, 2, 2),
                                 ncol=2)),
               rd=0)
print(crd3)

#! RCBD
rcbd1 <- gerexp(mu=0,
                sigma=1,
                ef=list(f1=c(5, 1, 1)),
                eb=c(2, 1, 1),
                rd=1,
                type='RCBD')
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

#! RCBD - Multivariated
rcbd2 <- gerexp(mu=c(0, 2),
                sigma=matrix(c(1, 0, 0, 1),
                             ncol=2),
                ef=list(f1= matrix(c(1, 1, 5, 1, 1, 1),
                                   ncol=2)),
                eb=matrix(c(2, 1, 1, 2, 1, 1),
                          ncol=2),
                rd=1,
                type='RCBD')
print(rcbd2)

#! LSD
lsd1 <- gerexp(mu=30,
               sigma=1,
               ef=list(f1=c(1, 1, 10)),
               erow=c(1, 1, 1),
               ecol=c(1, 1, 1),
               rd=1,
               type='LSD')
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
fecrd1 <- gerexp(mu=30,
                 sigma=1,
                 ef=list(f1=c(1, 1, 3),
                         f2=c(1, 1)),
                 einter=c(3, 1, 1, 1, 1, 5),
                 rd=1,
                 type='FE')
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

#! FE - RCBD
fercbd1 <- gerexp(mu=30,
                  sigma=1,
                  ef=list(f1=c(1, 1, 1), 
                          f2=c(2, 3)),
                  eb=c(1, 3),
                  einter=c(1, 15, 1, 1, 5, 1),
                  rd=1,
                  type='FE')
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
felsd1 <- gerexp(mu=30,
                 sigma=1,
                 ef=list(f1=c(1, 1), 
                         f2=c(2, 3)),
                 erow=c(1, 3, 2, 1),
                 ecol=c(2, 2, 1, 1),
                 einter=c(1, 15, 1, 1),
                 rd=1,
                 type='FE')
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

#! SPLIT PLOT - CRD
splitcrd1 <- gerexp(mu=30,
                    sigma=1,
                    sigmap=1,
                    ef=list(f1=c(1, 1), 
                            f2=c(2, 3)),
                    einter=c(1, 15, 1, 1),
                    rd=1,
                    type='SPE')
print(splitcrd1)

mod <- lm(Y1 ~ X1*X2 + X1:r,
          data=splitcrd1$dfm)  #X1:r erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ X1*X2 + Error(X1:r),
            data=splitcrd1$dfm)
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#SPLIT - RCBD
splitrcbd1 <- gerexp(mu=30,
                     sigma=1,
                     sigmap=1,
                     ef=list(f1=c(1, 1), 
                             f2=c(2, 3),
                             f3=c(1, 1, 1)),
                     eb=c(1,2,3,3),
                     einter=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                              1, 2, 1, 1,
                              1, 1, 1, 1, 1, 1,
                              1, 1, 3, 3, 3, 3),
                     rd=1,
                     type='SPE')
print(splitrcbd1)

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block,
          data=splitrcbd1$dfm)  #X1:Block erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block),
            data=splitrcbd1$dfm)
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT - LSD
splitlsd1 <- gerexp(mu=30,
                    sigma=1,
                    sigmap=1,
                    ef=list(f1=c(1, 1, 2), 
                            f2=c(2, 3, 1)),
                    einter=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                    erow = c(1, 1, 1),
                    ecol = c(1, 1, 1),
                    rd=1,
                    type='SPE')
print(splitlsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column,
          data=splitlsd1$dfm)  #X1:Row:Column erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column),
            data=splitlsd1$dfm)
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk) 

##########   TESTANDO AS FUNÇÕES GRÁFICAS  #############
##! CRD - Only one factor
# Static
plot(crd1)

#- Dynamic
#+ Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
crd1p <- update(crd1, 
                r = 3,
                ef = list(f1 = c(1,1)))
plot(crd1p,
     dynamic = TRUE)

#+ With real levels
plot(crd1p,
     newlevels = c('Tilápia','Pacu'),
     dynamic = TRUE)
 
