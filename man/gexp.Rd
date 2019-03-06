%===============================================================================
% Name   : plot.gexp
% Author : J.C.Faria
% Date   : 17/02/2019 20:41:53
% Version: v1
%===============================================================================
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

\description{ O pacote provê recursos computacionais úteis no planejamento e
 modelagem de experimentos estruturados no ambiente R.

 A funçao genérica S3 \code{gexp} foi criada para possibilitar
 planejar, criar e modelar experimentos estruturados, ou seja, sob delineamento
 experimental. Na modelagem permite simular resultados de experimentos com
 possibilidade do usuário informar os efeitos e o(s) erro(s) aleatório(s).
 Os delineamentos são: Completely Randomized Design (CRD),
 Randomized Complete Block Design (RCBD) e
 Latin Squares Design (LSD). Os tipos de experimentos são: Factorial Experiment
 (FE) e Split-plot Experiment (SPE).

 Os experimentos podem ser gerados com uma
 ou mais variável de resposta, no último caso, pode se impor uma estrutura de
 covariância estre elas. É também possível planejar experimentos com o auxílio 
 das funções gráficas para uso em planejamento a partir de figuras ou fotos da 
 área experimental.
}

\usage{
gexp(x, ...)

\method{gexp}{default}(
       mu        = NULL,
       err       = NULL,
       errp      = NULL,
       r         = 2L,
       fl        = NULL,
       blkl      = NULL,
       rowl      = NULL,
       coll      = NULL,
       fe        = NULL,
       inte      = NULL,
       blke      = NULL,
       rowe      = NULL,
       cole      = NULL,
       contrasts = NULL,
       type      = c('CRD','RCBD','LSD','FE','SPE'),
       nrand     = 1,
       round     = 2L,
       random    = FALSE,
       ...)

\method{gexp}{crd}(
       mu,
       err,
       r,
       fl,
       fe,
       contrasts,
       round,
       random)

\method{gexp}{rcbd}(
       mu,
       err,
       r,
       fl,
       blkl,
       fe,
       blke,
       contrasts,
       round,
       random)

\method{gexp}{lsd}(
       mu, 
       err,
       fl,
       rowl,
       coll,
       fe,
       rowe,
       cole,
       nrand,
       contrasts,
       round,
       random)

\method{gexp}{fe}(
       mu,
       err,
       r,
       fl,
       blkl,
       rowl,
       coll,
       fe,
       inte,
       blke,
       rowe,
       cole,
       nrand,
       contrasts,
       round,
       random)

\method{gexp}{spe}(
       mu,
       err,
       errp,
       r,
       fl,
       blkl,
       rowl,
       coll,
       fe,
       inte,
       blke,
       rowe,
       cole,
       contrasts,
       nrand,
       round,
       random)
}

\arguments{
  \item{mu}{É um escalar numérico, ou um vetor para dados multivariados (MD, 
   mais de uma variável resposta), que representa a média de cada fator.}
  \item{err}{É um vetor, ou matrix para MD, que representa o erro experimental. 
   O valor default é NULL, ou seja, para cada variável de resposta é adicionado
   um erro normal com média 0 e variância 1: rmvnorm(sigma = diag(length(mu)).}
  \item{errp}{É um vetor, ou uma matrix para MD, do erro associado as parcelas 
   no caso de type ser igual a SPE (Split-Plot Experiments). O valor default é 
   NULL, ou seja, para cada variável de resposta é adicionado um erro normal 
   com média 0 e variância 1: rmvnorm(sigma = diag(length(mu)).}
  \item{r}{É um escalar do número de repetições do experimento.}
  \item{fl}{Lista de um vetor de caracteres, ou uma matrix (MD). É uma lista 
   dos nomes dos fatores.}
  \item{blkl}{Lista de um vetor de caracteres, ou uma matriz para MD, dos nomes 
   dos blocos}
  \item{rowl}{Lista de um vetor de caracteres, ou uma matriz para MD, dos nomes 
   das linhas no caso de type ser igual a LSD (Latin Square Design).}
  \item{coll}{Lista de um vetor de caracteres, ou uma matriz para MD, dos nomes 
   das colunas no caso de type ser igual a LSD (Latin Square Design).}
  \item{fe}{É um vetor numérico, ou uma matrix (MD). É uma lista do efeito de 
   um fator.}
  \item{inte}{É um vetor numérico, ou uma matrix para MD, dos efeitos da 
   interação.}
  \item{blke}{É um vetor numérico, ou uma matriz para MD, dos efeitos dos 
   blocos.}
  \item{rowe}{É um vetor numérico, ou uma matriz para MD, dos efeitos das 
   linhas no caso de type ser igual a LSD (Latin Square Design).}
  \item{cole}{É um vetor numérico, ou uma matriz para MD, dos efeitos das 
   colunas no caso de type ser igual a LSD (Latin Square Design).}
  \item {contrasts}{...}
  \item{type}{É um vetor de strings que contém o tipo de delineamento e tipo de 
   experimento a ser utilizado:
   Completely Randomized Design (CRD), Randomized Complete Block Design (RCBD),
   Latin Squares Design (LSD), Factorial Experiment (FE) and Split-plot 
   Experiment (SPE). CRD é o default.}
  \item{nrand}{É um escalar numérico utilizado especificamente no Latin Squares 
   Design (LSD) para embaralhar os tratamentos nas linhas e colunas.}
  \item{round}{É um escalar numérico para arredondamentos da variável de 
   resposta.}
  \item{random}{É um argumento lógico quando o objetivo é planejar experimentos 
   para que ocorra a casualização dos tratamentos nas unidades experimentais.
   FALSE é o default.}
  \item{\dots}{Further arguments (required by generic).}
}
\value{
  The method \code{gexp} returns a list of class \code{gexp} with the slots:
  \item{X}{É a matriz de incidência do delineamento.}
  \item{Z}{É a matriz de incidência do erro da parcela principal no caso de type 
   ser igual a SPLIT.}
  \item{Y}{É um vetor, ou uma matrix para MD, com os valores da(s) variável(eis) 
   aleatória(s).}
  \item{dfm}{É um data.frame com todas as informações do experimento: 
   tratamentos, repetições e a variável aleatória de resposta.}
}

\references{
  Ferreira, Daniel Furtado. 2008. \emph{Estatística Multivariada}. 1 ed. 
  Lavras: Ed. UFLA. 

  Aquino, Luiz Henrique. \emph{Técnica Experimental com Animais I}. Apostila da 
  disciplina ``Técnica Experimental com Animais'' da Universidade Federal de 
  Lavras, 1992.

  Rencher, Alvin C. and Schaalje, Bruce G. 2007. \emph{Linear Models in 
  Statistics, second edition}. Hoboken: John Wiley \& Sons.

  Naes, T.; Aastveit, A.H.; Sahni, N.S. 2007. "Analysis of split-plot designs: 
  An Overview and Comparison of Methods". \emph{Qual. Reliab. Engng. Int.} 23, 
  801-820. 
}

\seealso{
  \code{\link{plot.gexp.crd}}
}

\examples{
#!___________________________
#! Qualitative Factor(s) (QL)
#!___________________________

#! Completely Randomized Design (CRD)
#! 1 factor - CRD - QLF
# Nonsense(experimental error = 0)
# Yi = mu + fe + e
r <- 2  # (repet. number)
fln <- 3  # (factor levels number)

crd00 <- gexp(mu=0,
              r=r,
              fe=list(f1=c(1, 2, 3)),
              err=matrix(0, 
                         nrow=r*fln),
              round=0)
crd00
print(crd00)
str(crd00)
summary(crd00)

#! 1 factor - CRD - QL
# Nonsense(error is 0) 
# Yi = mu + fe + e
r <- 3  # (repet. number)
fln <- 5  # (factor levels number)

crd01 <- gexp(mu=1,
              r=r,
              fe=list(f1=c(0, 2, 4, 6, 8)),
              err=matrix(0,
                         nrow=r*fln),
              round=2)
summary(crd01)

#! 1 factor - CRD - QL
# Default error: rmvnorm(sigma = diag(ncol(as.matrix([[fe]]))))
crd_1f <- gexp(mu=1,
               r=3,
               fe=list(f1=c(1, 1, 5, 1, 1)),
               fl=list(Treat=LETTERS[1:5]),
               round=2)

summary(crd_1f)

#! Binomial error - CRD - QL
e_binom <- as.matrix(rbinom(n=15,
                            size=5,
                            prob=0.1))

crd_bin <- gexp(mu=20,
                err=e_binom,
                r=5,
                fe=list(f1=c(1, 4, 1)))

summary(crd_bin)

mod <- aov(Y1 ~ X1,
           data=crd_bin$dfm)

shapiro.test(mod$res)

#! Factorial Experiment (FE) - CRD - QL
# 5 factors (f1..f5, level^factor): 3^1 * 2^1 * 3^1 * 2^1 * 4^1 * 2 = 360 experimental units
# 5 factors (f1..f5, level^factor): 3^2 * 2^2 * 5^1 = 180 * 2 = 360 experimental units
crd_fe <- gexp(mu=0,
               r=2,
               fe=list(f1=c(1, 1, 5),
                       f2=c(1, 1),
                       f3=c(2, 2, 1),
                       f4=c(1, 5),
                       f5=c(1, 2, 3, 4, 5)),
               fl=list(A=paste('a',
                               1:3,
                               sep=''),
                       B=paste('b',
                               1:2,
                               sep=''),
                       C=paste('c',
                               1:3,
                               sep=''),
                       D=paste('d',
                               1:2,
                               sep=''),
                       E=paste('e',
                               1:5,
                               sep='')),
               round=0)
summary(crd_fe)

#! Factorial Experiment (FE) - Multivariated - CRD - QL
# Error = 0 - Nonsense (you can easily undertand the effects)
crd_femn <- gexp(mu=c(0, 0),
                 r=1,
                 err=mvtnorm::rmvnorm(n=3^1 * 2^1 * 1,
                                      sigma=matrix(c(0, 0,
                                                     0, 0),
                                                   ncol=2)),
                                   #Y1 Y2
                 fe=list(f1=matrix(c(0, 3,  #X1  X1
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
summary(crd_femn)

#! Factorial Experiment (FE) - Multivariated - CRD - QL
# Using default error
set.seed(30)
crd_femd <- gexp(mu=c(0, 2),
                 r=3,
                 fe=list(f1=matrix(c(1, 1,
                                     5, 1,
                                     1, 1),
                                   ncol=2,
                                   byrow=TRUE),
                         f2=matrix(c(1, 3,
                                     2, 2),
                                   ncol=2,
                                   byrow=TRUE)),
                 round=1)
summary(crd_femd)

set.seed(30)
crd_femi <- gexp(mu=c(0, 2),
                 err=mvtnorm::rmvnorm(n=3^1 * 2^1 * 3,
                                      sigma=matrix(c(1, 0,  # The same that the default error
                                                     0, 1),
                                                   ncol=2)),
                 r=3,
                 fe=list(f1=matrix(c(1, 1,
                                     5, 1, 
                                     1, 1),
                                   ncol=2,
                                   byrow=TRUE),
                         f2=matrix(c(1, 3, 
                                     2, 2),
                                   ncol=2,
                                   byrow=TRUE)),
                 round=1)
summary(crd_femi)

crd_femd$dfm[, 4:5] # Use of the dafault error
crd_femi$dfm[, 4:5] # Use of the user error (same as the default!)
crd_femd$dfm[, 4:5] == crd_femi$dfm[, 4:5]

#! Factorial Experiment (FE) - With interaction - CRD - QL
fe_crd <- gexp(mu=30,
               fe=list(f1=c(1, 1, 3),
                       f2=c(1, 1)),
               fl=list(A=paste('a',
                               1:3,
                               sep=''),
                       B=paste('b',
                               1:2,
                               sep='')),
               inte=c(3, 1, 1, 1, 1, 5),
               round=1,
               type='FE')
summary(fe_crd)

#! Split-plot Experiment (SPE) - CRD - QL
split_crd <- gexp(mu=30,
                  fe=list(f1=c(1, 1),
                          f2=c(2, 3)),
                  fl=list(P=paste('p',
                                  1:2,
                                  sep=''),
                          SP=paste('sp',
                                   1:2,
                                   sep='')),
                  inte=c(1, 15, 1, 1),
                  round=1,
                  type='SPE')
summary(split_crd)  #TODO: Não está errado mas é melhor se o resultado for orddnado pela parcela, não pela subparcela!

#! Randomized Complete Block Design (RCBD) - QL
# 1 factor, 3 blocks
rcbd <- gexp(mu=0,
             fe=list(f1=c(5, 1, 1)),
             fl=list(TR=LETTERS[1:3]),
             blke=c(1, 2, 3),
             blkl=list(BLK=paste('B',
                                 1:3,
                                 sep='')),
             round=1,
             type='RCBD')
summary(rcbd) # TODO: não dá o plot

#! Factorial Experiment (FE) - RCBD - QL
fe_rcbd <- gexp(mu=30,
                fe=list(f1=c(1, 1, 1),
                        f2=c(2, 3)),
                blke=c(1, 3),
                inte=c(1, 15, 1, 1, 5, 1),
                round=1,
                type='FE')
summary(fe_rcbd)

#! Multivariated - RCBD - QL
rcbd_m <- gexp(mu=c(0, 2),
               fe=list(f1= matrix(c(1, 1,
                                    5, 1, 
                                    1, 1),
                                  ncol=2,
                                  byrow=TRUE)),
               blke=matrix(c(2, 1,
                             1, 2,
                             1, 1),
                           ncol=2,
                           byrow=TRUE),
               round=1,
               type='RCBD')
summary(rcbd_m) 

#! Split-plot Experiment (SPE) - RCBD - QL
split_rcbd <- gexp(mu=30,
                   fe=list(f1=c(1, 1),
                           f2=c(2, 3),
                           f3=c(1, 1, 1)),
                   fl=list(A=paste('a',
                                   1:2,
                                   sep=''),
                           B=paste('b',
                                   1:2,
                                   sep=''),
                           C=paste('c',
                                   1:3,
                                   sep='')),
                   blke=c(1, 2),
                   blkl=list(BLK=paste('B',
                                       1:2,
                                       sep='')),
                   inte=c(1, 15, 1, 1, 1, 3, 4, 2, 1, 1, 4, 1,
                          1, 2, 1, 1,
                          1, 1, 1, 1, 1, 1,
                          1, 1, 3, 3, 3, 3),
                   round=1,
                   type='SPE')
summary(split_rcbd)  #TODO: é melhor sair o label real do bloco, não o nome da lista!

#! Latin Square Design (LSD) - QL
lsd <- gexp(mu=30,
            fe=list(f1=c(1, 1, 10)),
            rowe=c(1, 1, 1),
            cole=c(1, 1, 1),
            rowl=list(Row=paste('r',
                                1:3,
                                sep='')),
            coll=list(Col=paste('c',
                                1:3,
                                sep='')),
            round=1,
            type='LSD')
summary(lsd)

#! Factorial Experiment (FE) - LSD - QL
fe_lsd <- gexp(mu=30,
               fe=list(f1=c(1, 1),
                       f2=c(2, 3)),
               rowe=c(1, 3, 2, 1),
               cole=c(2, 2, 1, 1),
#               rowl=list(Row=paste('r',
#                                   1:4,
#                                   sep='')),
#               coll=list(Col=paste('c',
#                                   1:4,
#                                   sep='')),
               inte=c(1, 15, 1, 1),
               round=1,
               type='FE')
summary(fe_lsd)  #TODO: não entendi mensagem de erro ao informar rowl e coll!

#! Split-plot Experiment (SPE) - LSD - QL
split_lsd <- gexp(mu=30,
                  fe=list(f1=c(1, 1, 2),
                          f2=c(2, 3, 1)),
                  fl=list(P=paste('p',
                                  1:3,
                                  sep=''),
                          SP=paste('sp',
                                   1:3,
                                   sep='')),
                  inte=c(1, 15, 1, 1, 1, 1, 1, 1, 1),
                  rowe = c(1, 1, 1),
                  cole = c(1, 1, 1),
                  rowl=list(Row=paste('r',
                                      1:3,
                                      sep='')),
                  coll=list(Col=paste('c',
                                      1:3,
                                      sep='')),
                  round=1,
                  type='SPE')
summary(split_lsd)

#!___________________________
#! Quantitative Factor(s) (QT)
#!___________________________

# CRD - Orthogonal polynomials
level <- c(0, 10, 20, 30)
cont_crd <- contr.poly(length(level))

# Linear effect
crd_lo <- gexp(mu=NULL,
               r=4,
                         # B0 B1 B2 B3  Linear only
               fe=list(f1=c(2, 5, 0, 0)),
               fl=list(Dose=ordered(level)),
               contrasts=list(f1=cont_crd))
summary(crd_lo)
plot(Y1 ~ Dose,
     crd_lo$dfm)

# Quadratic effect
crd_qo <- gexp(mu=NULL,
               r=4,
                         # B0 B1 B2 B3  quadratic
               fe=list(f1=c(2, 0, 5, 0)),
               fl=list(Dose=ordered(c(0, 10, 20, 30))),
               contrasts=list(f1=cont_crd))
summary(crd_qo)
plot(Y1 ~ Dose,
     crd_qo$dfm)

# Cubic effect
crd_co <- gexp(mu=NULL,
               r=4,
                         # B0 B1 B2 B3  cubic
               fe=list(f1=c(2, 0, 0, 5)),
               fl=list(Dose=ordered(c(0, 10, 20, 30))),
               contrasts=list(f1=cont_crd))
summary(crd_co)
plot(Y1 ~ Dose,
     crd_co$dfm)

# Not orthogonal polynomials
# Linear
cont_crd <- matrix(c(level,
                     level^2,
                     level^3),
                   ncol=3)

crd_l <- gexp(mu=NULL,
              r=4,
              fe=list(f1=c(2, 10, 0, 0)),
              fl=list(Dose=ordered(c(0, 10, 20, 30))),
              contrasts=list(f1=cont_crd))
summary(crd_l)
plot(Y1 ~ Dose,
     crd_l$dfm)

reg <- lm(Y1 ~ Dose + I(Dose^2) + I(Dose^3),
          data=crd_l$dfm)

summary(reg)

# Linear and quadratic
# No caso de dois ou mais fatores, basta informar o valor de Beta0 apenas para o primeiro fator.  #TODO: não entendi!
crd_lq <- gexp(mu=NULL,
               r=3,
               fe=list(f1=c(0, 10, 0, 0),  #linear
                       f2=c(0,  3, 4, 0)), #quadratic
               fl=list(P=ordered(level),
                       N=ordered(1:4)),
               contrasts=list(f1=cont_crd,
                              f2=cont_crd))
summary(crd_lq)

with(crd_lq$dfm,
     plot(Y1 ~ P))

with(crd_lq$dfm,
     plot(Y1 ~ N))

# Multivariated!
crd_m <- gexp(mu=NULL,
              r=4,                #L  Q
              fe=list(f1=matrix(c( 2,  2,
                                  10,  0,
                                   0, 10,
                                   0,  0),
                                ncol=2,
                                byrow=TRUE)),
              fl=list(Dose=ordered(level)),
              contrasts=list(f1=cont_crd))

with(crd_m$dfm,
     plot(Y1 ~ Dose))

with(crd_m$dfm,
     plot(Y2 ~ Dose))

# Orthogonal polynomios - RCBD
cont_rcbd <- contr.poly(4)

rcbd <- gexp(fe=list(f1=c(1, 3, 0, 0)),
             blke=c(1, 2, 3),
             r=2,
             fl=list(Dose=ordered(c(0, 2, 4, 6))),
             blkl=list(Blk=c('B1', 'B2', 'B3')),
             contrasts=list(f1=cont_rcbd,
                            Blk=diag(3)),
             type='RCBD')
summary(rcbd)

#! Hibrid: qualitative and quantitative factors in the same experiment - hb
# CRD
r <- 2
(error <- matrix(rep(0,
                     4^1*3^1*r),
                 ncol=1))

crd_hb <- gexp(mu=NULL,
               err=error,
               r=r,
               fe=list(f1=c(0, 1, 0),      # Qualitative
                       f2=c(2, 1, 0, 0)),  # Quantitative linear
               fl=list(Var=LETTERS[1:3],
                       Dose=ordered(level)),
               contrasts=list(f1=diag(3),
                              f2=cont_crd))
summary(crd_hb)

# RCBD - Tentativa - JCF
r <- 2
(error <- matrix(rep(0,
                     4^1*3^1*r),
                 ncol=1))

#TODO: Comentado para tentar R RCMD check/build
#TODO: Não doda em blocos?
#rcbd_hb <- gexp(mu=NULL,
#                err=error,
#                r=r,
#                fe=list(f1=c(0, 1, 0),      # Qualitative
#                        f2=c(2, 1, 0, 0)),  # Quantitative linear
#                fl=list(Var=LETTERS[1:3],
#                        Dose=ordered(level)),
#                blke=c(1, 2),
#                blkl=list(Blk=c('B1', 'B2')),
#                contrasts=list(f1=diag(3),
#                               f2=cont_crd,
#                               Blk=diag(2)),
#                type='RCBD')
#summary(rcbd_hb)

}
