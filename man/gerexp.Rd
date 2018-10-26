\name{gerexp}
\alias{gerexp}
\alias{gerexp.default}
\alias{gerexp.crd}
\alias{gerexp.rcbd}
\alias{gerexp.lsd}
\alias{gerexp.fe}
\alias{gerexp.spe}

\docType{package}

\title{Generator of Experiments}

\description{O pacote gerexp foi elaborado com o intuito principal de gerar variáveis aleatórias de acordo com a metodologia experimental e o tipo de efeito de cada tratamento. É possível também planejar experimentos com o auxílio das funções gráficas.
}

\usage{
gerexp(x, ...)

\method{gerexp}{default}(
       mu = NULL,
       sigma  = 1L,  
       sigmap = 1L,  
       r      = 2L,  
       ef     = NULL,
       einter = NULL,
       eb     = NULL,
       erow   = NULL,
       ecol   = NULL,
       type   = c('CRD','RCBD','LSD','FE','SPE'),
       nrand  = 1,   
       rd     = 2L,
       randomized = FALSE, 
       , ...)

\method{gerexp}{crd}(
       mu,
       sigma,
       r,
       ef,
       rd,
       randomized)

\method{gerexp}{rcbd}(
       mu,
       sigma,
       ef,
       eb,
       rd,
       randomized)
                                                
\method{gerexp}{lsd}(
       mu, 
       sigma, 
       ef,  
       erow, 
       ecol, 
       nrand,    
       rd,   
       randomized)

\method{gerexp}{fe}(
       mu,
       sigma, 
       r,   
       ef, 
       einter, 
       eb, 
       erow, 
       ecol, 
       nrand,    
       rd,   
       randomized)

\method{gerexp}{spe}(
       mu,
       sigma,
       sigmap,
       r,   
       ef, 
       einter, 
       eb, 
       erow, 
       ecol, 
       nrand,    
       rd,   
       randomized)
}

\arguments{
  \item{mu}{É um escalar numérico ou um vetor se tiver mais de uma variável resposta que representa a média geral.}
  \item{sigma}{É um escalar numérico ou uma matrix de variância-covariância para dados multivariados do erro experimental cujo a média é zero.}
  \item{sigmap}{É um escalar numérico ou uma matrix de variância-covariância para dados multivariados do erro da parcela principal no caso de type ser igual a SPE. A média também é zero.}
  \item{r}{É um escalar que se refere ao nº de repetições.}
  \item{ef}{É uma lista em que cada dimensão se refere ao efeito de um fator que pode ser um vetor numérico ou uma matriz (caso multivariado).}
  \item{einter}{É um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos da interação.}
  \item{eb}{É um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos dos blocos.}
  \item{erow}{É um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos das linhas no caso de type ser igual a LSD.}
  \item{ecol}{É um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos das colunas no caso de type ser igual a LSD.}
  \item{type}{É um vetor de strings que contém o tipo de delineamento a ser utilizado. CRD é o default.}
  \item{nrand}{É um escalar numérico utilizado especificamente no LSD para embaralhar os tratamentos nas linhas e colunas.}
  \item{rd}{É um escalar numérico para arredondamentos.}
  \item{randomized}{É um argumento lógico quando o objetivo é planejar experimentos para que ocorra a casualização dos tratamentos nas unidades experimentais. FALSE é o default.}
  \item{\dots}{Further arguments (required by generic).}
}
\value{
  The method \code{gerexp} returns a list of class \code{gerexp} with the slots:

  \item{X}{É a matriz de incidência do delineamento.}
  \item{Z}{É a matriz de incidência do erro da parcela principal no caso de type ser igual a SPLIT.}
  \item{Y}{É um vetor com os valores da variável aleatória ou uma matriz no caso de mais de uma variável aleatória.}
  \item{dfm}{É um data.frame com todas as informações do experimento como tratamentos, repetições e a variável aleatória.}
}

\references{

  Ferreira, Daniel Furtado. 2008. \emph{Estatística Multivariada}. 1 ed. Lavras: Ed. UFLA. 

  Aquino, Luiz Henrique. \emph{Técnica Experimental com Animais I}. Apostila da disciplina ``Técnica Experimental com Animais'' da Universidade Federal de Lavras, 1992.

  Rencher, Alvin C. and Schaalje, Bruce G. 2007. \emph{Linear Models in Statistics, second edition}. Hoboken: John Wiley \& Sons.

  Naes, T.; Aastveit, A.H.; Sahni, N.S. 2007. "Analysis of split-plot designs: An Overview and Comparison of Methods". \emph{Qual. Reliab. Engng. Int.} 23, 801-820. 

}

\seealso{
  \code{\link{plot.gerexp.crd}}
}

\examples{
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
}
