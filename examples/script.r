#########       CARREGANDO AS FUNÇÕES       ###########
R.utils::sourceDirectory('../R')
library(debug)
mtrace.off()
mtrace('plot.gerexp.crd')
mtrace(gexp.spe)
######### TESTANDO AS FUNÇÕES NÃO GRÁFICAS ############

#! CRD
#! 1 factor
crd_1 <- gexp(mu=10,
              r=3,
              ef=list(f1=c(1, -1)),
              round=2)

summary(crd_1)
str(crd_1$dfm)

mod <- lm(Y1 ~ X1,
          data=crd_1$dfm)

summary(mod)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')

summary(tk)

#! 5 factors
crd2 <- gexp(mu=0,
             r=2,
             ef=list(f1=c(1, 1, 5),
                     f2=c(1, 1),
                     f3=c(2, 2, 1),
                     f4=c(1, 5),
                     f5=c(1, 2, 3, 4, 5)),
             round=0)

summary(crd2)
str(crd2$dfm)
head(crd2$dfm)
tail(crd2$dfm)

#! CRD - Multivariated
# In this case it is necessary to inform the error!
# The default is univariate case.

crd3 <- gexp(mu=c(5, 2),
#             err=mvtnorm::rmvnorm(n=9,
#                                  sigma=matrix(c(1, 0,
#                                                 0, 1),
#                                               ncol=2)),
             r=3,
             ef=list(f1=matrix(c(1, 1,
                                 5, 1,
                                 1, 1),
                               ncol=2,
                               byrow=TRUE)),
             round=4)

summary(crd3)

#! CRD - with other contrasts
# Orthogonal polynomials
level <- c(0, 10, 20, 30)
cont_crd4 <- contr.poly(length(level))

# Linear
crd4 <- gexp(mu=NULL,
             r=4,
             ef=list(f1=c(2, 10, 0, 0)),
             factorsl=list(Dose=ordered(level)),
             contrasts=list(f1=cont_crd4))

crd44 <- crd4$dfm
crd44$dose <- as.numeric(as.character(crd44$Dose))
summary(crd44)

plot(Y1 ~ dose,
     crd44)

md <- with(crd44,
           tapply(Y1,
                  dose,
                  mean))

plot(md ~ names(md))

reg_crd <- lm(Y1 ~ Dose,
             crd44)

summary(reg_crd)
fitted(reg_crd)

x <- matrix(c(level,
              level^2,
              level^3),
              ncol=3)  

X <- cbind(rep(1, 4),
           x)

theta <- coef(reg_crd)

Z <- cbind(rep(1, 4), cont_crd4)

(betas <- solve(X)%*%Z%*%theta)

# Quadratic
crd5 <- gexp(mu=NULL,
             r=4,
             ef=list(f1=c(2, 1, 5, 0)),
             factorsl=list(Dose=ordered(c(0, 10, 20, 30))),
             contrasts=list(f1=cont_crd4))

crd55 <- crd5$dfm
crd55$dose <- as.numeric(as.character(crd55$Dose))

plot(Y1 ~ dose, crd55)

medias  <- with(crd55,
                tapply(Y1,
                       Dose,
                       mean))

plot(medias ~ names(medias))

reg_crd2 <- lm(Y1 ~ Dose,
              crd55)

summary(reg_crd2)
fitted(reg_crd2)

# Cubic
crd6 <- gexp(mu=NULL,
             r=4,
             ef=list(f1=c(2, 1, 1, 10)),
             factorsl=list(Dose=ordered(c(0, 10, 20, 30))),
             contrasts=list(f1=cont_crd4))

crd66 <- crd6$dfm
crd66$dose <- as.numeric(as.character(crd66$Dose))

plot(Y1 ~ dose,
     crd66)

medias  <- with(crd66,
                tapply(Y1,
                       Dose,
                       mean))

plot(medias ~ names(medias))

reg_crd3 <- lm(Y1 ~ Dose,
               crd66)

summary(reg_crd3)
fitted(reg_crd3)

#! Not orthogonal polynomios
#! Linear
cont_crd7 <- matrix(c(level,
                      level^2,
                      level^3),
                    ncol=3) 

crd7 <- gexp(mu=NULL,
             r=4,
             ef=list(f1=c(2, 10, 0, 0)),
             factorsl=list(Dose=ordered(c(0, 10, 20, 30))),
             contrasts=list(f1=cont_crd7))

crd77 <- crd7$dfm
crd77$dose <- as.numeric(as.character(crd77$Dose))

plot(Y1 ~ dose, crd77)

medias  <- with(crd77,
                tapply(Y1,
                       Dose,
                       mean))

plot(medias ~ names(medias))

reg_crd4 <- lm(Y1 ~ dose + I(dose^2) + I(dose^3),
               crd77)

summary(reg_crd4)
fitted(reg_crd4)

#! CRD - with other errors
# Binomial error
e_binom <- as.matrix(rbinom(n=15,
                            size=5,
                            prob=0.1))

crd8 <- gexp(mu=20,
             err=e_binom,
             r=5,
             ef=list(f1=c(1, 4, 1)))

summary(crd8)

mod8 <- aov(Y1 ~ X1,
            data=crd8$dfm)

par(mfrow=c(2, 2))
plot(mod8)
shapiro.test(mod8$res)

#! RCBD
rcbd1 <- gexp(mu=0,
              ef=list(f1=c(5, 1, 1)),
              eb=c(2, 1),
              round=1,
              type='RCBD')

str(rcbd1)

mod <- lm(Y1 ~ Block + X1,
          data=rcbd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)

#! RCBD without replication
rcbd2 <- gexp(mu=0,
              r=1,
              ef=list(f1=c(5, 1, 1)),
              eb=c(2, 1),
              round=1,
              type='RCBD')

str(rcbd2)

mod <- lm(Y1 ~ Block + X1,
          data=rcbd1$dfm)

anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)

#! RCBD - Multivariated
rcbd3 <- gexp(mu=c(0, 2),
#              err=mvtnorm::rmvnorm(n=18,
#                                   sigma=matrix(c(1, 0, 
#                                                  0, 1),
#                                                ncol=2)),
              ef=list(f1=matrix(c(1, 1, 
                                  5, 1, 
                                  1, 1),
                                ncol=2,
                                byrow=TRUE)),
              eb=matrix(c(2, 1, 
                          1, 2, 
                          1, 1),
                        ncol=2,
                        byrow=TRUE),
              round=1,
              type='RCBD')

str(rcbd3)
summary(rcbd3)

#! RCBD - With other contrasts
# Orthogonal polynomios
cont_rcbd4 <- contr.poly(4)  s

rcbd4 <- gexp(ef=list(f1=c(1, 3, 0, 0)),
              eb=c(1, 2, 3),
              r=1,
              factorsl=list(Dose=ordered(c(0, 2, 4, 6))),
              blocksl=list(bloco=c('B1', 'B2', 'B3')),
              contrasts=list(f1=cont_rcbd4,
                             bloco=diag(3)),
              type='RCBD')

rcbd44 <- rcbd4$dfm
rcbd44$dose <- as.numeric(as.character(rcbd44$Dose))

plot(Y1 ~ dose, rcbd44)

medias  <- with(rcbd44,
                tapply(Y1,
                       Dose,
                       mean))

plot(medias ~ names(medias))

regrcbd4 <- lm(Y1 ~ bloco + Dose,
               rcbd44)

summary(regrcbd4)
fitted(regrcbd4)

reg_rcbd4m <- lm(medias ~ poly(x=c(0, 2, 4, 6),
                               deg=1))

summary(reg_rcbd4m)
fitted(reg_rcbd4m)

#! LSD 
lsd1 <- gexp(mu=30,
             ef=list(f1=c(1, 1, 10)),
             erow=c(1, 1, 1),
             ecol=c(1, 1, 1),
             round=1,
             type='LSD')

str(lsd1)

mod <- lm(Y1 ~ Row + Column + X1,
          data=lsd1$dfm)

anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')

summary(tk)

#! FE - CRD
fe_crd1 <- gexp(mu=30,
                ef=list(f1=c(1, 1, 3),
                        f2=c(1, 1)),
                eint=c(3, 1, 1, 1, 1, 5),
                round=1,
                type='FE')

str(fe_crd1)

mod <- lm(Y1 ~ X1*X2,
          data=fe_crd1$dfm)

anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)
summary(tk)

#! FE - CRD: Tree factors
fe_crd2 <- gexp(mu=30,
                ef=list(f1=c(1, 1, 3),
                        f2=c(1, 1),
                        f3=c(2, 1, 1, 1)),
                eint=rep(1,50),
                round=1,
                type='FE')

str(fe_crd2)

#! FE - RCBD
fe_rcbd1 <- gexp(mu=30,
                 ef=list(f1=c(1, 1, 1),
                         f2=c(2, 3)),
                 eb=c(1, 3),
                 eint=c(1, 15, 1, 1, 5, 1),
                 round=1,
                 type='FE')

str(fe_rcbd1)

mod <- lm(Y1 ~ Block + X1*X2,
          data=fe_rcbd1$dfm)

anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=1)

summary(tk)

#! FE - LSD
fe_lsd1 <- gexp(mu=30,
                ef=list(f1=c(1, 1),
                        f2=c(2, 3)),
                erow=c(1, 3, 2, 1),
                ecol=c(2, 2, 1, 1),
                eint=c(1, 15, 1, 1),
                round=1,
                type='FE')

str(fe_lsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2,
          data=fe_lsd1$dfm)

anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)

summary(tk)

#! SPLIT PLOT - CRD: two level
split_crd1 <- gexp(mu=30,
                   ef=list(f1=c(1, 1),
                           f2=c(2, 3)),
                   eint=c(1, 15, 1, 1),
                   round=1,
                   type='SPE')

str(split_crd1)

mod <- lm(Y1 ~ X1*X2 + X1:r,
          data=split_crd1$dfm)  #X1:r erro(a) parcela

anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ X1*X2 + Error(X1:r),
                             data=split_crd1$dfm))

summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tk)

#! SPLIT PLOT - CRD: tree level
split_crd2 <- gexp(mu=30,
                   ef=list(f1=c(1, 1),
                           f2=c(2, 3),
                           f3=c(1, 1, 1)),
                   eint=c(1, 15, 1, 1, rep(1, 24)),
                   round=1,
                   type='SPE')

str(split_crd2)

#SPLIT - RCBD
split_rcbd1 <- gexp(mu=30,
                    ef=list(f1=c(1, 1),
                            f2=c(2, 3),
                            f3=c(1, 1, 1)),
                    eb=c(1, 2, 3, 3),
                    eint=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                           1, 2, 1, 1,
                           1, 1, 1, 1, 1, 1,
                           1, 1, 3, 3, 3, 3),
                    round=1,
                    type='SPE')

str(split_rcbd1)

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block,
          data=split_rcbd1$dfm)  #X1:Block erro(a) parcela

anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block),
                             data=split_rcbd1$dfm))

summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)

summary(tk)

#! SPLIT - LSD
split_lsd1 <- gexp(mu=30,
                   ef=list(f1=c(1, 1, 2),
                           f2=c(2, 3, 1)),
                   eint=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                   erow=c(1, 1, 1),
                   ecol=c(1, 1, 1),
                   round=1,
                   type='SPE')

str(split_lsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column,
          data=split_lsd1$dfm)  #X1:Row:Column erro(a) parcela

anova(mod)

mod1 <- suppressWarnings(aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column),
                             data=split_lsd1$dfm))

summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=1)

summary(tk)

#! SPLIT - LSD: tree factors 
split_lsd2 <- gexp(mu=30,
                   ef=list(f1=c(1, 1, 2),
                           f2=c(2, 3, 1),
                           f3=c(1, 1)),
                   eint=rep(1, 39),
                   erow=c(1, 1, 1),
                   ecol=c(1, 1, 1),
                   round=1,
                   type='SPE')

summary(split_lsd2)

########## COM FATORES FORNECIDOS PELOS USUÁRIOS #######
crd_f1 <- update(crd_1,
                 factorsl=list('Species'=c('Oreochromis', 'Piaractus')))

str(crd_f1)

crd_f2 <- update(crd_1,
                 r=2,
                 ef=list(f1=c(1, 1, 1)),
                 factorsl=list('Breed'=c('Dexter', 'Guernsey', 'Jersey')))

str(crd_f2)

rcbd_f1 <- update(rcbd1,
                  factorsl=list('Variety'=c('A', 'B', 'C')))

str(rcbd_f1)

lsd_f1 <- update(lsd1,
                 factorsl=list('Horse'=c('English', 'QM', 'ML')),
                 rowsl=list('Linha'=1:3),
                 colsl=list('Coluna'=1:3))

str(lsd_f1)

fe_crd_f1 <- update(fe_crd1,
                    factorsl=list('Ration'=c('0', '10', '20'),
                                  'Chicken'=c('Cob', 'Ross')))

str(fe_crd_f1)

fe_rcbd_f1 <- update(fe_rcbd1,
                     newfactors=list('Ration'=c('0', '10', '20'),
                                     'Chicken'=c('Cob', 'Ross')))

str(fe_rcbd_f1)

fe_lsd_f1 <- update(fe_lsd1,
                    factorsl=list('Ration'=c('0', '10'),
                                  'Chicken'=c('Cob', 'Ross')))

str(fe_lsd_f1)

split_crd_f1 <- update(split_crd1,
                       factorsl=list('Forage'=c('A', 'B'),
                                     'Period'=c('I1', 'I2')))

str(split_crd_f1)

split_lsd_f1 <- update(split_lsd1,
                       factorsl=list('Ration'=c('0', '10', '20'),
                                     'Chicken'=c('Cob', 'Ross', 'Caneludo')),
                       rowsl=list('Linha'=1:3),
                       colsl=list('Coluna'=1:3))

str(split_lsd_f1)

##########   TESTANDO AS FUNÇÕES GRÁFICAS  #############
##! CRD - Only one factor
# Static
plot(crd_1)
plot(crd_f1)
plot(crd_f2)

# Dynamic
# Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
crd_p1 <- update(crd_1,
                 r=3,
                 ef=list(f1=c(1, 1)))

plot(crd_p1,
     dynamic=TRUE)

# With real level
plot(crd_f1,
     dynamic=TRUE)

##! RCBD - Only one factor 
# Static
plot(rcbd1)

# Dynamic
# Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
rcbd_1p <- update(rcbd1,
                  r=1,
                  ef=list(f1=c(1, 1, 1, 1)),
                  eb=c(1, 1, 1))
plot(rcbd_1p,
     dynamic=TRUE)

##! LSD - Only one factor 
# Static
plot(lsd1)
plot(lsd_f1)

# Dynamic
# Obs: Clique apenas uma vez com o botão esquerdo do mouse sobre a unidade experimental. Após terminar, clique com o botão esquerdo do mouse para encerrar!!
plot(lsd_f1,
     dynamic=TRUE)

##! FE:CRD  
# Static
plot(fe_crd1)

##! FE:RCBD 
# Static
plot(fe_rcbd1)

##! FE:LSD 
# Static
plot(fe_lsd1)

##! SPE:CRD
plot(split_crd1)
plot(split_crd_f1)
plot(split_crd2)

##! SPE: RCBD
plot(split_rcbd1)

##! SPE: LSD
plot(split_lsd1)
