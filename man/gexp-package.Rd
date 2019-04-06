\name{Generator of Experiments}
\alias{gerexp-package}

\docType{package}

\title{gexp: Generator of Experiments}
\description{The package \code{gexp} was created to make it possible to plan, create
  and to model structured experiments, that is, under a experimental design. In the modeling 
  it is possible to simulate results of experiments with
  possibility of user to report effects and random error(s).
  The designs are: Completely Randomized Design (CRD),
  Randomized Complete Block Design (RCBD) and Latin Squares Design (LSD). 
  The types of experiments are: Factorial Experiment (FE) and Split-plot Experiment (SPE).

  The experiments can be generated with one or more response variables,
  In the latter case, a strict covariance structure can be imposed.
  It is also possible to plan experiments using the graphic functions
  for use in planning from pictures or pictures of the experimental area.

  The possible uses are multiple: in the planning it makes possible to
  distribution and randomization of treatments and experimental units;
  in the data analysis allows to generate experiments for application in evaluations
  individual and can also be used to generate experiments for validations
  of new computational resources in the area of structured data analysis.

  In summary, the package provides computational resources useful in planning and
  modeling of structured experiments in the R.
}

\details{In some situations, we are interested in simulating a variable
  randomized according to the experimental procedure where the differences between
  treatments are predetermined. In a completely randomized design
  with two treatments for example, we may have an interest in simulating a variable
  random whose treatment A will have a 1-deviation effect and treatment B a
  effect of 3 deviations from a given overall average. In addition,
  may be interested in imposing a pre-established error structure for
  purposes evaluation in the various analysis strategies.
}
%\description{O pacote \code{gexp} foi criado para possibilitar planejar, criar
% e modelar experimentos estruturados, ou seja, sob delineamento
% experimental. Na modelagem permite simular resultados de experimentos com
% possibilidade do usuário informar os efeitos e o(s) erro(s) aleatório(s).
% Os delineamentos são: Completely Randomized Design (CRD),
% Randomized Complete Block Design (RCBD) e
% Latin Squares Design (LSD). Os tipos de experimentos são: Factorial Experiment
% (FE) e Split-plot Experiment (SPE).
%
% Os experimentos podem ser gerados com uma ou mais variável de resposta,
% no último caso, pode se impor uma estrutura de covariância estre elas.
% É também possível planejar experimentos com o auxílio das funções gráficas
% para uso em planejamento a partir de figuras ou fotos da área experimental.
%
% Os possíveis usos são múltiplos: no planejamento possibilita efetuar a
% distribuição e aleatorização dos tratamentos e unidades experimentais;
% na análise de dados permite gerar experimentos para aplicação em avaliações
% individuais e pode ser usado também para gerar experimentos para validações
% de novos recursos computacionais na área de análise de dados estruturados.
%
% Em síntese, o pacote provê recursos computacionais úteis no planejamento e
% modelagem de experimentos estruturados no ambiente R.
%}
%
%\details{Em algumas situações, temos o interesse em simular uma variável
%aleatória de acordo com o procedimento experimental onde as diferenças entre
%tratamentos sejam pré-determinadas. Em um delineamento inteiramente ao acaso
%com dois tratamentos por exemplo, podemos ter interesse em simular uma variável
%aleatória cujo o tratamento A terá um efeito de 1 desvio e o tratamento B um
%efeito de 3 desvios em relação a uma determinada média geral. Adicionalmente,
%pode-se estar interessado em impor uma estrutura de erros pré-estabelecida para
%fins avaliação nas diversas estratégias de análise.
%}
%
\author{
Ivan Bezerra Allaman <ivanalaman@gmail.com> and José Cláudio Faria <joseclaudio.faria@gmail.com>.

Maintainer: Ivan Bezerra Allaman <ivanalaman@gmail.com>
}

\keyword{package}
