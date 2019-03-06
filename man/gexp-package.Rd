\name{Generator of Experiments}
\alias{gerexp-package}

\docType{package}

\title{gexp: Generator of Experiments}

\description{O pacote \code{gexp} foi criado para possibilitar planejar, criar
 e modelar experimentos estruturados, ou seja, sob delineamento
 experimental. Na modelagem permite simular resultados de experimentos com
 possibilidade do usuário informar os efeitos e o(s) erro(s) aleatório(s).
 Os delineamentos são: Completely Randomized Design (CRD),
 Randomized Complete Block Design (RCBD) e
 Latin Squares Design (LSD). Os tipos de experimentos são: Factorial Experiment
 (FE) e Split-plot Experiment (SPE).

 Os experimentos podem ser gerados com uma ou mais variável de resposta,
 no último caso, pode se impor uma estrutura de covariância estre elas.
 É também possível planejar experimentos com o auxílio das funções gráficas
 para uso em planejamento a partir de figuras ou fotos da área experimental.

 Os possíveis usos são múltiplos: no planejamento possibilita efetuar a
 distribuição e aleatorização dos tratamentos e unidades experimentais;
 na análise de dados permite gerar experimentos para aplicação em avaliações
 individuais e pode ser usado também para gerar experimentos para validações
 de novos recursos computacionais na área de análise de dados estruturados.

 Em síntese, o pacote provê recursos computacionais úteis no planejamento e
 modelagem de experimentos estruturados no ambiente R.
}

\details{Em algumas situações, temos o interesse em simular uma variável
aleatória de acordo com o procedimento experimental onde as diferenças entre
tratamentos sejam pré-determinadas. Em um delineamento inteiramente ao acaso
com dois tratamentos por exemplo, podemos ter interesse em simular uma variável
aleatória cujo o tratamento A terá um efeito de 1 desvio e o tratamento B um
efeito de 3 desvios em relação a uma determinada média geral. Adicionalmente,
pode-se estar interessado em impor uma estrutura de erros pré-estabelecida para
fins avaliação nas diversas estrat´gis de análise.
}

\author{
Ivan Bezerra Allaman <ivanalaman@gmail.com> and José Cláudio Faria <joseclaudio.faria@gmail.com>.

Maintainer: Ivan Bezerra Allaman <ivanalaman@gmail.com>
}

\keyword{package}
