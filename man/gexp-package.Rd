\name{Generator of Experiments}
\alias{gexp-package}

\docType{package}

\title{gexp: Generator of Experiments}
\description{The package \code{gexp} was created to make it possible to plan, create and to model structured experiments, that is, under a experimental design. In the modeling it is possible to simulate results of experiments with possibility of user to report effects and random error(s). The designs are: Completely Randomized Design (CRD), Randomized Complete Block Design (RCBD) and Latin Squares Design (LSD). The types of experiments are: Simple (SIMPLE), Factorial Experiment (FE) and Split-plot Experiment (SPE).

  The experiments can be generated with one or more response variables, in the latter case, a strict covariance structure can be imposed. It is also possible to plan experiments using the graphic functions for use in planning from pictures or pictures of the experimental area.

  The possible uses are multiple: in the planning it makes possible to distribution and randomization of treatments and experimental units; in the data analysis allows to generate experiments for application in evaluations individual and can also be used to generate experiments for validations of new computational resources in the area of structured data analysis.

  In summary, the package provides computational resources useful in planning and modeling of structured experiments in the \code{R}.
}

\details{In some situations, we are interested in simulating a variable randomized according to the experimental procedure where the differences between treatments are predetermined. In a completely randomized design with two treatments for example, we may have an interest in simulating a variable random whose treatment A will have a 1-deviation effect and treatment B a effect of 3 deviations from a given overall average. In addition, may be interested in imposing a pre-established error structure for purposes evaluation in the various analysis strategies.
}

\author{
Ivan Bezerra Allaman <ivanalaman@gmail.com> and Jose Claudio Faria <joseclaudio.faria@gmail.com>.

Maintainer: Ivan Bezerra Allaman <ivanalaman@gmail.com>
}
