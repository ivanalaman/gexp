#! Testing the main fuction: ger.design

#library(debug)
#mtrace(gerexp)
#mtrace.off()

#! DIC
#! 1 factor
dic1 <- gerexp(mu=0,
               s2=1,
               r=3,
               f=list(f1=c(1, 1, 5)),
               rd=0)
print(dic1)

mod <- lm(Y1 ~ X1,
          data=dic1)
anova(mod) 

par(mfrow=c(2,2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! 5 factors
dic2 <- gerexp(mu=0,
               s2=1,
               r=2,
               f=list(f1=c(1, 1, 5),
                      f2=c(1, 1),
                      f3=c(2, 2, 1),
                      f4=c(1, 5),
                      f5=c(1, 2, 3, 4, 5)),
               rd=0)
print(dic2)

#! DIC - Multivariated
dic3 <- gerexp(mu=c(0, 2),
               s2=matrix(c(1, 0, 0, 1),
                         ncol=2),
               r=3,
               f=list(f1=matrix(c(1, 1, 5, 1, 1, 1),
                                ncol=2),
                      f2=matrix(c(1, 3, 2, 2),
                                ncol=2)),
               rd=0)
print(dic3)

#! DBC
dbc1 <- gerexp(mu=0,
               s2=1,
               f=list(f1=c(5, 1, 1)),
               b=c(2, 1, 1),
               rd=1,
               type='DBC')
print(dbc1)

mod <- lm(Y1 ~ Block + X1,
          data=dbc1)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)

#! DBC - Multivariated
dbc2 <- gerexp(mu=c(0, 2),
               s2=matrix(c(1, 0, 0, 1),
                         ncol=2),
               f=list(f1= matrix(c(1, 1, 5, 1, 1, 1),
                                 ncol=2)),
               b=matrix(c(2, 1, 1, 2, 1, 1),
                        ncol=2),
               rd=1,
               type='DBC')
print(dbc2)

#! DQL
dql1 <- gerexp(mu=30,
               s2=1,
               f=list(f1=c(1, 1, 10)),
               erow=c(1, 1, 1),
               ecol=c(1, 1, 1),
               rd=1,
               type='DQL')
print(dql1)

mod <- lm(Y1 ~ Row + Column + X1,
          data=dql1)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)
plot(tk)
 
#! FAT - DIC
fatdic1 <- gerexp(mu=30,
                  s2=1,
                  f=list(f1=c(1, 1, 3),
                         f2=c(1, 1)),
                  inter=c(3, 1, 1, 1, 1, 5),
                  rd=1,
                  type='FAT')
print(fatdic1)

mod <- lm(Y1 ~ X1*X2,
          data=fatdic1)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)
summary(tuk)
plot(tuk)
 
#! FAT - DBC
fatdbc1 <- gerexp(mu=30,
                  s2=1,
                  f=list(f1=c(1, 1, 1), 
                         f2=c(2, 3)),
                  b=c(1, 3),
                  inter=c(1, 15, 1, 1, 5, 1),
                  rd=1,
                  type='FAT')
print(fatdbc1)

mod <- lm(Y1 ~ Block + X1*X2,
          data=fatdbc1)
anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk)

#! FAT - DQL
fatdql1 <- gerexp(mu=30,
                  s2=1,
                  f=list(f1=c(1, 1), 
                         f2=c(2, 3)),
                  erow=c(1, 3, 2, 1),
                  ecol=c(2, 2, 1, 1),
                  inter=c(1, 15, 1, 1),
                  rd=1,
                  type='FAT')
print(fatdql1)

mod <- lm(Y1 ~ Row + Column + X1*X2,
          data=fatdql1)
anova(mod)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT PLOT - DIC
splitdic1 <- gerexp(mu=30,
                    s2=1,
                    s2sp=1,
                    f=list(f1=c(1, 1), 
                           f2=c(2, 3)),
                    inter=c(1, 15, 1, 1),
                    rd=1,
                    type='SPLIT')
print(splitdic1)

mod <- lm(Y1 ~ X1*X2 + X1:r,
          data=splitdic1)  #X1:r erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ X1*X2 + Error(X1:r),
            data=splitdic1)
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)
 
#SPLIT - DBC
splitdbc1 <- gerexp(mu=30,
                    s2=1,
                    s2sp=1,
                    f=list(f1=c(1, 1), 
                           f2=c(2, 3),
                           f3=c(1, 1, 1)),
                    b=c(1,2,3,3),
                    inter=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                            1, 2, 1, 1,
                            1, 1, 1, 1, 1, 1,
                            1, 1, 3, 3, 3, 3),
                    rd=1,
                    type='SPLIT')
print(splitdbc1)

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block,
          data=splitdbc1)  #X1:Block erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block),
            data=splitdbc1)
summary(mod1)

par(mfrow=c(2,2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=2)

summary(tuk)
plot(tuk)

#! SPLIT - DQL
splitdql1 <- gerexp(mu=30,
                    s2=1,
                    s2sp=1,
                    f=list(f1=c(1, 1, 2), 
                           f2=c(2, 3, 1)),
                    inter=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                    erow = c(1, 1, 1),
                    ecol = c(1, 1, 1),
                    rd=1,
                    type='SPLIT')
print(splitdql1)

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column,
          data=splitdql1)  #X1:Row:Column erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column),
            data=splitdql1)
summary(mod1)

par(mfrow=c(2, 2))
plot(mod)

tuk <- TukeyC(mod,
              which='X1:X2',
              fl1=1)

summary(tuk)
plot(tuk)
