\name{lba}
\alias{lba}
\alias{lba.matrix}
\alias{lba.table}
\alias{lba.formula}
\alias{lba.ls}
\alias{lba.mle}
\alias{lba.ls.fe}
\alias{lba.mle.fe}
\alias{lba.ls.logit}
\alias{lba.mle.logit}
\alias{lba.default}

\docType{package}

\title{Latent Budget Analysis (LBA) for Compositional Data}

\description{Latent budget analysis (LBA) is a method for the analysis of contingency tables, from where the compositional data is derived. It is used to understand the relationship between the table rows and columns, where the rows denote the categories of the explanatory variable and the columns denote the categories of the response variable.
}

\usage{
lba(obj, ...)

\method{lba}{matrix}(obj,
    A           = NULL,
    B           = NULL,
    K           = 1L,
    cA          = NULL,
    cB          = NULL,
    logitA      = NULL,
    logitB      = NULL,
    omsk        = NULL,
    psitk       = NULL,
    S           = NULL,
    T           = NULL,
    row.weights = NULL,
    col.weights = NULL,
    tolG        = 1e-10,
    tolA        = 1e-05,
    tolB        = 1e-05,
    itmax.unide = 1e3,
    itmax.ide   = 1e3,
    trace.lba   = TRUE,
    toltype     = "all",
    method      = c("ls", "mle"),
    what        = c("inner","outer"), ...)

\method{lba}{table}(obj,
    A           = NULL,
    B           = NULL,
    K           = 1L,
    cA          = NULL,
    cB          = NULL,
    logitA      = NULL,
    logitB      = NULL,
    omsk        = NULL,
    psitk       = NULL,
    S           = NULL,
    T           = NULL,
    row.weights = NULL,
    col.weights = NULL,
    tolG        = 1e-10,
    tolA        = 1e-05,
    tolB        = 1e-05,
    itmax.unide = 1e3,
    itmax.ide   = 1e3,
    trace.lba   = TRUE,
    toltype     = "all",
    method      = c("ls", "mle"),
    what        = c("inner","outer"), ...)

\method{lba}{formula}(formula, data,
    A           = NULL,
    B           = NULL,
    K           = 1L,
    cA          = NULL,
    cB          = NULL,
    logitA      = NULL,
    logitB      = NULL,
    omsk        = NULL,
    psitk       = NULL,
    S           = NULL,
    T           = NULL,
    row.weights = NULL,
    col.weights = NULL,
    tolG        = 1e-10,
    tolA        = 1e-05,
    tolB        = 1e-05,
    itmax.unide = 1e3,
    itmax.ide   = 1e3,
    trace.lba   = TRUE,
    toltype     = "all",
    method      = c("ls", "mle"),
    what        = c("inner","outer"), ...)

\method{lba}{ls}(obj,
    A           ,
    B           ,
    K           ,
    row.weights ,
    col.weights ,
    tolA        ,
    tolB        ,
    itmax.unide ,
    itmax.ide   ,
    trace.lba   ,
    what        , ...)

\method{lba}{mle}(obj,
    A           ,
    B           ,
    K           ,
    tolG        ,
    tolA        ,
    tolB        ,
    itmax.unide ,
    itmax.ide   ,
    trace.lba   ,
    toltype     ,
    what        , ...)

\method{lba}{ls.fe}(obj,
    A           ,
    B           ,
    K           ,
    cA          ,
    cB          ,
    row.weights ,
    col.weights ,
    itmax.ide   ,
    trace.lba   , ...)

\method{lba}{mle.fe}(obj,
    A          ,
    B          ,
    K          ,
    cA         ,
    cB         ,
    tolG       ,
    tolA       ,
    tolB       ,
    itmax.ide  ,
    trace.lba  ,
    toltype    , ...)

\method{lba}{ls.logit}(obj,
    A           ,
    B           ,
    K           ,
    cA          ,
    cB          ,
    logitA      ,
    logitB      ,
    omsk        ,
    psitk       ,
    S           ,
    T           ,
    row.weights ,
    col.weights ,
    itmax.ide   ,
    trace.lba   , ...)

\method{lba}{mle.logit}(obj,
    A          ,
    B          ,
    K          ,
    cA         ,
    cB         ,
    logitA     ,
    logitB     ,
    omsk       ,
    psitk      ,
    S          ,
    T          ,
    itmax.ide  ,
    trace.lba  , ...)
}

\arguments{
 \item{obj,formula}{The function is generic, accepting some forms of the principal argument for specifying a two-way frequency table. Currently accepted forms are matrix, data frame (coerced to frequency tables), objects of class \code{"xtabs"} or \code{"table"} and one-sided formulae of the form \code{Col1 + Col2 + \dots + Coln ~ Row1 + Row2 + \dots + Rown}, where \code{Rown} and \code{Coln} are nth row (the mixing parameters) and column variable (the latent components).}
 \item{data}{A data frame containing variables in \code{formula}.}
 \item{A}{The starting value of a (I x K) matrix containing the mixing parameters, if given. The default is \code{NULL}, producing random starting values.}
 \item{B}{The starting value of a (J x K) matrix containing the latent components, if given. The default is \code{NULL}, producing random starting values.}
 \item{K}{Integer giving the number of latent budgets chosen by the user. The default is \code{1}.}
 \item{cA}{The value of a (I x K) matrix containing the constraints on the mixing parameters. Fixed constraints are the values themselves which are numbers in the [0,1] interval. The optional equality constraints are indicated by an integer starting from 2, such that parameters that must be equal have the same integer. The default is \code{NULL}, indicating no constraints.}
 \item{cB}{The value of a (J x K) matrix containing the constraints on the latent components. Fixed constraints are the values themselves which are numbers in the [0,1] interval. The optional equality constraints are indicated by an integer starting from 2, such that parameters that must be equal have the same integer. The default is \code{NULL}, indicating no constraints.}
 \item{logitA}{Design (IxS) matrix for row-covariates. The first column contains 1´s, indicating a constant covariate. The entries may be continuous or dummy coded values.}
 \item{logitB}{Design (JxT) matrix for column-covariates. The entries may be continuous or dummy coded values.}
 \item{omsk}{A (SxK) matrix giving the starting values for the multinomial logit parameters of the row covariates. The default is \code{NULL}, producing random starting values.}
 \item{psitk}{A (TxK) matrix giving the starting values for the multinomial logit parameters of the column covariates. The default is \code{NULL}, producing random starting values.}
 \item{S}{Number of row-covariates. The default is \code{NULL}.}
 \item{T}{Number of column-covariates. The default is \code{NULL}.}
 \item{row.weights}{A vector with the same number of rows of the matrix of the weighted least squares method. If is \code{NULL} (default), the weights are \deqn{\sqrt{n_{i+}/n_{++}}}.}
 \item{col.weights}{A vector with the same number of columns of the matrix of the weighted least squares method. If is \code{NULL} (default), the weights are \deqn{1/\sqrt{n_{+j}/n_{++}}}.}
 \item{tolG}{A tolerance value for judging when convergence has been reached. It is based on the estimated likelihood ratio statistics G2. The default is \code{1e-10}.}
 \item{tolA}{A tolerance value for judging when convergence has been reached. When the one-iteration change in the maximum of the absolute value of the  element wise difference of the estimated matrices A is less than tolA. The default is \code{1e-05}.}
 \item{tolB}{A tolerance value for judging when convergence has been reached. When the one-iteration change in the maximum of the absolute value of the  element wise difference of the estimated matrices B is less than tolB. The default is \code{1e-05}.}
 \item{itmax.unide}{Maximum number of iterations performed by the mle or ls method, if convergence is not achieved, before identification parameters. The default is \code{1e3}.}
 \item{itmax.ide}{Maximum number of iterations performed by the mle or ls method in the identification process. Is used too when the constrained fixed, equality and logit are required. The default is \code{1e3}.}
 \item{trace.lba}{Logical, indicating whether the \code{base} function \code{optim} and \code{constrOptim.nl} from package \pkg{alabama}, will trace their results. The default is \code{TRUE}.}
 \item{toltype}{String indicating which kind of tolerance to be used. That is, the EM algorithm stops updating and considers the maximum log-likelihood  to have been found. Their types are: \code{"all"} when the one-iteration change in the estimated likelihood ratio statistics G2 is less than tolG, and the one-iteration change in the maximum of the absolute value of the element wise difference of the estimated matrices A is less than tolA and the same for estimated matrices B with respect to tolB; \code{"G2"} when the only one-iteration change in the estimated likelihood ratio statistics G2 is less than tolG; \code{"ab"} when only the one-iteration change in the maximum of the absolute value of the element wise difference of the estimated matrices A is less than tolA and the same for estimated matrices B with respect to tolB. \code{toltype} works only for method = \code{"mle"}. The default is \code{"all"}. The ls method uses only \code{"ab"} as tolerance limit.}
 \item{method}{String indicating which kind of estimating method. They are: \code{"ls"} when least squares, either weighted or ordinary, method is 
used; \code{"mle"} when maximum likelihood method is used. The default is \code{"ls"}.}
 \item{what}{String indicating which kind identified solutions for mixing parameters and latent budgets matrices. They are: the \code{"inner"} extreme solution and the \code{"outer"} extreme solution. The default is \code{"inner"}.}
 \item{\dots}{Further arguments (required by generic).}
}

\value{
 The method \code{lba.ls} and \code{lba.mle} returns a list of class \code{lba.ls} and \code{lba.mle} respectively with the slots:

 \item{P}{The compositional data matrix which is formed by dividing the raw data matrix by their corresponding total, its rows are called observed budgets.}
 \item{pij}{Matrix whose rows are the expected budgets.}
 \item{residual}{Residual matrix \code{P - pij.}}
 \item{A}{(I x K) matrix of the unidentified the mixing parameters.}
 \item{B}{(J x K) matrix of the unidentified the latent components.}
 \item{Aoi}{(I x K) matrix of the identified mixing parameters, they may be either the inner extreme values or the outer extreme values.}
 \item{Boi}{(J x K) matrix of the identified latent componentes, they may be either the inner extreme values or the outer extreme values.}
 \item{rescB}{(J x K) matrix of the rescaled latent components.}
 \item{pk}{Budget proportions.}
 \item{val_func}{Value of least squared or likelihood function achieved.}
 \item{iter_unide}{Number of unidentified iterations.}
 \item{iter_ide}{Number of identified iterations.}

 The method \code{lba.ls.fe} and \code{lba.mle.fe} returns a list of class \code{lba.ls.fe} and \code{lba.mle.fe} respectively with the slots:

 \item{P}{The compositional data matrix which is formed by dividing the raw data matrix by their corresponding row total, its rows are called observed budgets.}
 \item{pij}{Matrix whose rows are the expected budgets.}
 \item{residual}{Residual matrix \code{P - pij.}}
 \item{A}{(I x K) matrix of the unidentified the mixing parameters.}
 \item{B}{(J x K) matrix of the unidentified the latent components.}
 \item{rescB}{(J x K) matrix of the rescaled latent components.} 
 \item{pk}{Budget proportions.}
 \item{val_func}{Value of least squared or likelihood function achieved.}
 \item{iter_ide}{Number of identified iteractions.}

 The method \code{lba.ls.logit} and \code{lba.mle.logit} returns a list of class \code{lba.ls.logit} and \code{lba.mle.logit} respectively with the slots:  

 \item{P}{The compositional data matrix which is formed by dividing the raw data matrix by their corresponding total, its rows are called observed budgets.}
 \item{pij}{Matrix whose rows are the expected budgets.}
 \item{residual}{Residual matrix \code{P - pij.}}
 \item{A}{(I x K) matrix of the unidentified the mixing parameters.}
 \item{B}{(J x K) matrix of the unidentified the latent componentes.}
 \item{rescB}{(J x K) matrix of the rescaled latent components.} 
 \item{pk}{Budget proportions.}
 \item{val_func}{Value of least squared or likelihood function achieved.}
 \item{iter_ide}{Number of identified iterations.}
 \item{omsk}{A (SxK) matrix giving estimated values of the multinomial logit parameters of the row covariates.}
 \item{psitk}{A (TxK) matrix giving the estimated values for the multinomial logit parameters of the column covariates.}
}

\references{
 Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley \& Sons.

 de Leeuw, J., and van der Heijden, P.G.M. 1988. "The analysis of time-budgets with a latent time-budget model". In E. Diday (Ed.), \emph{Data Analysis and Informatics V}. pp. 159-166. Amsterdam: North-Holland.

 de Leeuw, J., van der Heijden, P.G.M., and Verboon, P. 1990. "A latent time budget model". \emph{Statistica Neerlandica}. 44, 1, 1-21.

 Dempster, A.P., Laird, N.M., and Rubin, D.B. 1977. "Maximum likelihood from incomplete data via the EM algorithm". \emph{Journal of the Royal Statistical Society, Series}. 39, 1-38.

 van der Ark, A.L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of comositional data}. Ph.D. Thesis University of Utrecht.

 van der Heijden, P.G.M., Mooijaart, A., and de Leeuw, J. 1992. "Constrained latent budget analysis". In P.V. Marsden (Ed.), \emph{Sociological Methodology} pp. 279-320. Cambridge: Blackwell Publishers.
}

\note{
 The user has two options to entry the data: the raw data and the tabulated data. If the raw data is imported, he may indicate which, among the variables, comprises the row and which the column variable and let the \code{lba.formula} function make the tabulation. The user may also tabulate the data with the available functions in R. Recalling that if this second option is used, the object must be of the class \code{xtabs}, \code{table} or \code{matrix}. If the user imports the tabulated data, the class is, in general, \code{data.frame} and so, it is necessary to transform the object data into a \code{matrix}.

 The function \code{lba} uses EM algorithm to maximise the latent budget model log-likelihood function; the Active Constraints Methods (ACM) to minimise either the weighted least squares (wls), or ordinary least squares (ols) functions; and "BFGS" variable metric method in \code{constrOptim.nl} function of \pkg{alabama} package and in \code{optim} function of \pkg{stats} package used in identification for K >= 3, in constraint algorithm for ls method, in multinomial logit constraints and in some parts of constraining for mle method. Depending on the starting parameters, those algorithms may only locate a local, rather than global, maximum.  This becomes more and more of a problem as K, the number of latent budgets, increases. It is therefore highly advisable to run \code{lba} multiple times until you are relatively certain that you have located the global maximum log-likelihood or the global minimum least squares.
}

\seealso{
 \code{\link{goodnessfit}}, \code{\link{summary.lba.ls}},\code{\link{summary.lba.mle}},\code{\link{plotlba}},\code{\link{plotcorr}}
}

\examples{
data('votB')

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(parties ~ city,
           votB,
           K = 2)
ex1 

# Already tabulated data? Ok!
data('PerfMark') 

\dontrun{
ex2 <- lba(as.matrix(PerfMark),
           K = 2,
           what='outer')
ex2
}
# Using LS method (default) with constraint
# Fixed constraint to mixing parameters
cakiF1 <- matrix(c(0.2, NA, NA,
                   NA , NA,0.2,
                   NA , NA,0.2,
                   0.3, NA, NA,
                   0.2, NA, NA,
                   NA , NA, NA),
                 byrow = TRUE,
                 ncol  = 3)  

# K = 3
\dontrun{
exf1 <- lba(parties ~ city,
            votB,
            cA = cakiF1,
            K = 3)
exf1 
}
# Using LS method (default) with LOGIT constrain
data('housing')

# Make cross-table to matrix design.
tbh <- xtabs(value ~ Influence + Housing, housing)

Xis <- model.matrix(~ Housing*Influence,
                    tbh,
                    contrasts=list(Housing='contr.sum',
                                   Influence='contr.sum'))

tby <- xtabs(value ~ Satisfaction + Contact, housing)

Yis <- model.matrix(~ Satisfaction*Contact,
                    tby,
                    contrasts=list(Satisfaction='contr.sum',
                                   Contact='contr.sum'))[,-1]

S <- 12
T <- 5

tabs <- xtabs(value ~ interaction(Housing,
                                  Influence) + interaction(Satisfaction,
                                  Contact),
              housing)
\dontrun{
exlogit2 <- lba(tabs,
           K = 2,
           logitA = Xis,
           logitB = Yis,
           S = S,
           T = T,
           trace.lba=FALSE) 
exlogit2
}
}
