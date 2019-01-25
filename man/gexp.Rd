\name{gexp}
\alias{gexp}
\alias{gexp.default}
\alias{gexp.crd}
\alias{gexp.rcbd}
\alias{gexp.lsd}
\alias{gexp.fe}
\alias{gexp.spe}

\docType{package}

\title{Generator of Experiments}

\description{A funçao genérica S2 gexp foi criada para possibilitar planejar, criar e modelar experimentos estruturados, ou seja, simulando os delineamentos experimentais, no ambiente R. Os delineamentos são: Completely Randomized Design (CRD), Randomized Complete Block Design (RCBD) e Latin Squares Design (LSD). Os tipos de experimentos são: Factorial Experiment (FE) e Split-plot Experiment (SPE). Os experimentos podem ser gerados com uma ou mais variável de resposta, no último caso, pode se impor uma estrutura de covariância estre elas. É também possível planejar experimentos com o auxílio das funções gráficas para uso em planejamento a partir de figuras ou fotos da área experimental.

 Ela pacote provê recursos computacionais úteis no planejamento e modelagem de experimentos estruturados no ambiente R.
}

\usage{
gexp(x, ...)

\method{gexp}{default}(
       mu     = NULL,
       err    = 1L,  
       errp   = 1L,  
       r      = 2L,  
       ef     = NULL,
       eint   = NULL,
       eb     = NULL,
       erow   = NULL,
       ecol   = NULL,
       type   = c('CRD','RCBD','LSD','FE','SPE'),
       nrand  = 1,   
       round  = 2L,
       random = FALSE, 
       ...)

\method{gexp}{crd}(
       mu,
       err,
       r,
       ef,
       round,
       random)

\method{gexp}{rcbd}(
       mu,
       err,
       ef,
       eb,
       round,
       random)
                                                
\method{gexp}{lsd}(
       mu, 
       err, 
       ef,  
       erow, 
       ecol, 
       nrand,    
       round,
       random)

\method{gexp}{fe}(
       mu,
       err, 
       r,   
       ef, 
       eint, 
       eb, 
       erow, 
       ecol, 
       nrand,    
       round,
       random)

\method{gexp}{spe}(
       mu,
       err,
       errp,
       r,   
       ef, 
       eint, 
       eb, 
       erow, 
       ecol, 
       nrand,    
       round,
       random)
}

\arguments{
  \item{mu}{É um escalar numérico, ou um vetor para dados multivariados (MD, mais de uma variável resposta), que representa a média de cada fator.}
  \item{err}{É um vetor, ou matrix para MD, que representa o erro experimental. O valor default é NULL, ou seja, para cada variável de resposta é adicionado 
  um erro normal com média 0 e variância 1: rmvnorm(sigma = diag(length(mu)).}
  \item{errp}{É um vetor, ou uma matrix para MD, do erro associado as parcelas no caso de type ser igual a SPE (Split-Plot Experiments). 
  O valor default é NULL, ou seja, para cada variável de resposta é adicionado um erro normal com média 0 e variância 1: rmvnorm(sigma = diag(length(mu)).}
  \item{r}{É um escalar do número de repetições do experimento.}
  \item{ef}{É um vetor numérico, ou uma matrix (MD). É uma lista do efeito de um fator.}
  \item{eint}{É um vetor numérico, ou uma matrix para MD, dos efeitos da interação.}
  \item{eb}{É um vetor numérico, ou uma matriz para MD, dos efeitos dos blocos.}
  \item{erow}{É um vetor numérico, ou uma matriz para MD, dos efeitos das linhas no caso de type ser igual a LSD (Latin Square Design).}
  \item{ecol}{É um vetor numérico, ou uma matriz para MD, dos efeitos das colunas no caso de type ser igual a LSD (Latin Square Design).}
  \item{type}{É um vetor de strings que contém o tipo de delineamento e tipo de experimento a ser utilizado: 
  Completely Randomized Design (CRD), Randomized Complete Block Design (RCBD), Latin Squares Design (LSD), 
  Factorial Experiment (FE) and Split-plot Experiment (SPE). CRD é o default.}
  \item{nrand}{É um escalar numérico utilizado especificamente no Latin Squares Design (LSD) para embaralhar os tratamentos nas linhas e colunas.}
  \item{round}{É um escalar numérico para arredondamentos da variável de resposta.}
  \item{random}{É um argumento lógico quando o objetivo é planejar experimentos para que ocorra a casualização dos tratamentos nas unidades experimentais.
  FALSE é o default.}
  \item{\dots}{Further arguments (required by generic).}
}
\value{
  The method \code{gexp} returns a list of class \code{gexp} with the slots:
  \item{X}{É a matriz de incidência do delineamento.}
  \item{Z}{É a matriz de incidência do erro da parcela principal no caso de type ser igual a SPLIT.}
  \item{Y}{É um vetor, ou uma matrix para MD,  com os valores da(s) variável(eis) aleatória(s).}
  \item{dfm}{É um data.frame com todas as informações do experimento: tratamentos, repetições e a variável aleatória.}
}

\references{
  Ferreira, Daniel Furtado. 2008. \emph{Estatística Multivariada}. 1 ed. Lavras: Ed. UFLA. 

  Aquino, Luiz Henrique. \emph{Técnica Experimental com Animais I}. Apostila da disciplina ``Técnica Experimental com Animais'' da Universidade Federal de Lavras, 1992.

  Rencher, Alvin C. and Schaalje, Bruce G. 2007. \emph{Linear Models in Statistics, second edition}. Hoboken: John Wiley \& Sons.

  Naes, T.; Aastveit, A.H.; Sahni, N.S. 2007. "Analysis of split-plot designs: An Overview and Comparison of Methods". \emph{Qual. Reliab. Engng. Int.} 23, 801-820. 
}

\seealso{
  \code{\link{plot.gexp.crd}}
}

\examples{
#! Completely Randomized Design (CRD)
#! 1 factor - nonsense
#! Yi = mu + ef + e
crd00 <- gexp(mu=0,
              r=2,
              ef=list(f1=c(1, 0, 0)),
              err=matrix(0, 
                         nrow=2*3),
              round=0)
summary(crd00)

#! 1 factor - nonsense also
#! Yi = mu + ef + e
crd01 <- gexp(mu=1,
              r=3,
              ef=list(f1=c(0, 2, 4, 6, 8)),
              err=matrix(0, 
                         nrow=3*5),
              round=2)
summary(crd01)

#! 1 factor - more useful
crd1 <- gexp(mu=1,
             r=3,
             ef=list(f1=c(1, 1, 5, 1, 1)),
             round=2)

crd1
summary(crd1)

mod <- lm(Y1 ~ X1,
          data=crd1$dfm)
anova(mod) 

par(mfrow=c(2, 2))
plot(mod)

library(TukeyC)
tk <- TukeyC(mod,
             which='X1')
summary(tk)

#! Factorial Experiment (FE) under Completely Randomized Design (CRD)
#! 5 factors (f1..f5, level^factor): 3^1 * 2^1 * 3^1 * 2^1 * 4^1 * 2 = 360 experimental units
#! 5 factors (f1..f5, level^factor): 3^2 * 2^2 * 5^1 = 180 * 2 = 360 experimental units
crd2 <- gexp(mu=0,
             r=2,
             ef=list(f1=c(1, 1, 5),
                     f2=c(1, 1),
                     f3=c(2, 2, 1),
                     f4=c(1, 5),
                     f5=c(1, 2, 3, 4, 5)),
             round=0)
summary(crd2)

#! Factorial Experiment (FE) under Completely Randomized Design (CRD) - Multivariated
#! Nonsense (to easily undertand the effects)
crd3 <- gexp(mu=c(0, 0),
             r=1,
             err=mvtnorm::rmvnorm(n=3^1 * 2^1 * 1,
                                  sigma=matrix(c(0, 0,
                                                 0, 0),
                                               ncol=2)),
                               #Y1 Y2
             ef=list(f1=matrix(c(0, 3,  #X1  X1
                                 1, 4,  #X2  X2
                                 2, 5), #X3  X3
                               ncol=2,
                               byrow=TRUE),

                               #Y1 Y2
                     f2=matrix(c(0, 2,  #X1  X1
                                 1, 3), #X2  X2
                               ncol=2,
                               byrow=TRUE)),
             round=1)
summary(crd3)

#! More useful
set.seed(30)
crd30 <- gexp(mu=c(0, 2),
              r=3,
              ef=list(f1=matrix(c(1, 1,
                                  5, 1,
                                  1, 1),
                                ncol=2,
                                byrow=TRUE),
                      f2=matrix(c(1, 3,
                                  2, 2),
                                ncol=2,
                                byrow=TRUE)),
              round=1)
summary(crd30)

set.seed(30)
crd31 <- gexp(mu=c(0, 2),
              err=mvtnorm::rmvnorm(n=3^1 * 2^1 * 3,
                                   sigma=matrix(c(1, 0, 
                                                  0, 1),
                                                ncol=2)),
              r=3,
              ef=list(f1=matrix(c(1, 1, 
                                  5, 1, 
                                  1, 1),
                                ncol=2,
                                byrow=TRUE),
                      f2=matrix(c(1, 3, 
                                  2, 2),
                                ncol=2,
                                byrow=TRUE)),
              round=1)
summary(crd31)

crd30$dfm[, 4:5] # Use of the dafault error
crd31$dfm[, 4:5] # Use of the user error (same as the default)
crd30$dfm[, 4:5] == crd31$dfm[, 4:5]

#! Randomized Complete Block Design (RCBD)
rcbd1 <- gexp(mu=0,
              ef=list(f1=c(5, 1, 1)),
              eb=c(1, 2, 3),
              round=1,
              type='RCBD')
summary(rcbd1)

mod <- lm(Y1 ~ Block + X1,
          data=rcbd1$dfm)
anova(mod) 
cv(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)

#! Randomized Complete Block Design (RCBD) - Multivariated
rcbd2 <- gexp(mu=c(0, 2),
              ef=list(f1= matrix(c(1, 1, 
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
summary(rcbd2) 

#! Latin Squares Design (LSD)
lsd1 <- gexp(mu=30,
             ef=list(f1=c(1, 1, 10)),
             erow=c(1, 1, 1),
             ecol=c(1, 1, 1),
             round=1,
             type='LSD')
summary(lsd1)

mod <- lm(Y1 ~ Row + Column + X1,
          data=lsd1$dfm)
anova(mod) 
cv(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1')
summary(tk)

#! Factorial Experiment (FE)
#! in
#! Completely Randomized Design (CRD)
fe_crd1 <- gexp(mu=30,
                ef=list(f1=c(1, 1, 3),
                        f2=c(1, 1)),
                eint=c(3, 1, 1, 1, 1, 5),
                round=1,
                type='FE')
summary(fe_crd1)

mod <- lm(Y1 ~ X1*X2,
          data=fe_crd1$dfm)
anova(mod) 
cv(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)
summary(tk)

#! Factorial Experiment (FE)
#! in
#! Randomized Complete Block Design (RCBD)
fe_rcbd1 <- gexp(mu=30,
                 ef=list(f1=c(1, 1, 1), 
                         f2=c(2, 3)),
                 eb=c(1, 3),
                 eint=c(1, 15, 1, 1, 5, 1),
                 round=1,
                 type='FE')
summary(fe_rcbd1)

mod <- lm(Y1 ~ Block + X1*X2,
          data=fe_rcbd1$dfm)
anova(mod)
cv(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=1)

summary(tk)

#! Factorial Experiment (FE)
#! in
#! Latin Squares Design (LSD)
fe_lsd1 <- gexp(mu=30,
                ef=list(f1=c(1, 1), 
                        f2=c(2, 3)),
                erow=c(1, 3, 2, 1),
                ecol=c(2, 2, 1, 1),
                eint=c(1, 15, 1, 1),
                round=1,
                type='FE')
summary(fe_lsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2,
          data=fe_lsd1$dfm)
anova(mod)
cv(mod)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)

summary(tk)

#! Split-plot Experiment (SPE)
#! in
#! Completely Randomized Design (CRD)
split_crd1 <- gexp(mu=30,
                   ef=list(f1=c(1, 1), 
                           f2=c(2, 3)),
                   eint=c(1, 15, 1, 1),
                   round=1,
                   type='SPE')
summary(split_crd1)

mod <- lm(Y1 ~ X1*X2 + X1:r,
          data=split_crd1$dfm)  #X1:r erro(a) parcela
anova(mod)
cv(mod)

mod1 <- aov(Y1 ~ X1*X2 + Error(X1:r),
            data=split_crd1$dfm)
summary(mod1)
cv(mod1)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)

summary(tk)
plot(tk)

#! Split-plot Experiment (SPE)
#! in
#! Randomized Complete Block Design (RCBD)
split_rcbd1 <- gexp(mu=30,
                    ef=list(f1=c(1, 1), 
                            f2=c(2, 3),
                            f3=c(1, 1, 1)),
                    eb=c(1,2,3,3),
                    eint=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                           1, 2, 1, 1,
                           1, 1, 1, 1, 1, 1,
                           1, 1, 3, 3, 3, 3),
                    round=1,
                    type='SPE')
summary(split_rcbd1)

mod <- lm(Y1 ~ Block + X1*X2*X3 + X1:Block,
          data=split_rcbd1$dfm)  #X1:Block erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Block + X1*X2*X3 + Error(X1:Block),
            data=split_rcbd1$dfm)
summary(mod1)
cv(mod1)

par(mfrow=c(2,2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=2)

summary(tk)

#! Split-plot Experiment (SPE)
#! in
#! Latin Squares Design (LSD)
split_lsd1 <- gexp(mu=30,
                   ef=list(f1=c(1, 1, 2), 
                           f2=c(2, 3, 1)),
                   eint=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                   erow = c(1, 1, 1),
                   ecol = c(1, 1, 1),
                   round=1,
                   type='SPE')
summary(split_lsd1)

mod <- lm(Y1 ~ Row + Column + X1*X2 + X1:Row:Column,
          data=split_lsd1$dfm)  #X1:Row:Column erro(a) parcela
anova(mod)

mod1 <- aov(Y1 ~ Row + Column + X1*X2 + Error(X1:Row:Column),
            data=split_lsd1$dfm)
summary(mod1)
cv(mod1)

par(mfrow=c(2, 2))
plot(mod)

tk <- TukeyC(mod,
             which='X1:X2',
             fl1=1)

summary(tk)
}
