## TESTANDO A FUNÇÃO ger.design
library(debug)
mtrace(gerexp)
mtrace.off()
## DIC
# com 1 fator
ddic <- gerexp(mu=0,
               sig=1,
               r=3,
               f=list(f1=c(1, 1, 5)),
               roundd=0)

mod <- lm(Y1 ~ X1,
          ddic)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

# com 5 fatores
ddic <- gerexp(mu=0,
               sig=1,
               r=2,
               f=list(f1=c(1, 1, 5),
                      f2=c(1,1),
                      f3=c(2,2,1),
                      f4=c(1,5),
                      f5=c(1,2,3,4,5)),
               roundd=0)
ddic 

## DIC - multivariado
ddic <- gerexp(mu=c(0,2),
               sig=matrix(c(1,0,0,1),ncol=2),
               r=3,
               f=list(f1 = matrix(c(1, 1, 5, 1, 1, 1),ncol=2),
                      f2 = matrix(c(1, 3, 2, 2),ncol=2)),
               roundd=0)
ddic

# DBC
ddbc <- gerexp(mu=0,
               sig=1,
               f=list(f1=c(5, 1, 1)),
               b=c(2, 1, 1),
               roundd=1,
               type='DBC')

mod <- lm(Y1 ~ Block + X1, ddbc)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

# DQL
ddql <- gerexp(mu=30,
               sig=1,
               f=list(f1=c(1, 1, 10)),
               row=c(1, 1, 1),
               col=c(1, 1, 1),
               roundd=1,
               type='DQL')

mod <- lm(Y1 ~ Row + Column + X1,
          ddql)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)
 
# FAT - DIC
dfatdic <- gerexp(mu=30,
                  sig=1,
                  f = list(f1=c(1, 1, 3),
                           f2=c(1, 1)),
                  inter=c(3, 1, 1, 1, 1, 5),
                  roundd=1,
                  type='FAT')

mod <- lm(Y1 ~ X1*X2, dfatdic)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)
summary(tuk)
plot(tuk)
 
# FAT - DBC
dfatdbc <- gerexp(mu=30,
                  sig=1,
                  f=list(f1=c(1, 1, 1), 
                         f2=c(2, 3)),
                  b=c(1, 3),
                  inter=c(1, 15, 1, 1, 5, 1),
                  roundd=1,
                  type='FAT')

mod <- lm(Y1 ~ Block + X1*X2, dfatdbc)
anova(mod)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk)

#FAT - DQL
dfatdql <- gerexp(mu=30,
                  sig=1,
                  f=list(f1=c(1, 1), 
                         f2=c(2, 3)),
                  row=c(1, 3, 2, 1),
                  col=c(2, 2, 1, 1),
                  inter=c(1, 15, 1, 1),
                  roundd=1,
                  type='FAT')

mod <- lm(Y1 ~ Row + Column + X1*X2, dfatdql)
anova(mod)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#SPLIT - DIC
dsplitdic <- gerexp(mu=30,
                    sig=1,
                    sigplot=1,
                    f=list(f1=c(1, 1), 
                           f2=c(2, 3)),
                    inter=c(1, 15, 1, 1),
                    roundd=1,
                    type='SPLIT')

mod <- lm(Y1 ~ X1*X2 + X1:r, dsplitdic)#X1:r erro(a) parcela
anova(mod)
mod1 <- aov(Y1 ~ X1*X2 + Error(X1:r), dsplitdic)
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)
 
#SPLIT - DBC
dsplitdbc <- gerexp(mu=30,
                    sig=1,
                    sigplot=1,
                    f=list(f1=c(1, 1), 
                           f2=c(2, 3),
                           f3=c(1,1,1)),
                    b=c(1,2,3,3),
                    inter=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                            1,2,1,1,
                            1,1,1,1,1,1,
                            1,1,3,3,3,3),
                    roundd=1,
                    type='SPLIT')

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block, dsplitdbc)#X1:Block erro(a) parcela
anova(mod)
mod1 <- aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block), dsplitdbc)
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#SPLIT - DQL
dsplitdql <- gerexp(mu=30,
                    sig=1,
                    sigplot=1,
                    f=list(f1=c(1, 1, 2), 
                           f2=c(2, 3, 1)),
                    inter=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                    row = c(1,1,1),
                    col = c(1,1,1),
                    roundd=1,
                    type='SPLIT')

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column, dsplitdql)#X1:Row:Column erro(a) parcela
anova(mod)
mod1 <- aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column), dsplitdql)
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk)


 
 
