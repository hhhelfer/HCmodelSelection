\name{Model.Selection.Phase}
\alias{Model.Selection.Phase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct sets of well-fitting models as proposed by Cox, D. R. & Battey, H. S. (2017)}
\description{
This function tests low dimensional subsests of the set of retained variables from the reduction phase and any squared or interaction terms suggested at the exploratory phase. Lists of well-fitting models of each dimension are returned.
}
\usage{
Model.Selection.Phase(X,Y, list.reduction, family=binomial, signif=0.05,
                      square.terms=NULL, interaction.terms=NULL, modelSize=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Design matrix.}
  \item{Y}{Response vector.}
  \item{list.reduction}{Indices of variables that where chosen at the reduction phase.}
  \item{family}{A description of the error distribution and link function to be used in the model. For glm this can be a character string naming a family function, a family function or the result of a call to a family function. See \code{\link{family}} for more details.}
  \item{signif}{Significance level of the likelihood ratio test against the comprehensive model. The default is 0.01.}
  \item{square.terms}{Indices of squared terms suggested at the exploratory phase (See \code{\link{Exploratory.Phase}})}
  \item{interaction.terms}{Indices of pairs of variables suggested at the exploratory phase (See \code{\link{Exploratory.Phase}})}
  \item{modelSize}{Maximum size of the models to be tested. Curently the maximum is 7. If not provided a default is used.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{goodModels}{List of models that are in the confidence set of size 1 to modelSize.}
%% ...

}
\references{
Cox, D. R. and Battey, H. S. (2017). Large numbers of explanatory variables, a semi-descriptive analysis. \emph{Proceedings of the National Academy of Sciences}, 114(32), 8592-8595.

Battey, H. S. and Cox, D. R. (2018). Large numbers of explanatory variables: A probabilistic assessment. \emph{Proceedings of the Royal Society of London, A.}, to appear.
}
\author{
Hoeltgebaum, H. H.
}
\note{
%%  ~~further notes~~
}

\section{Acknowledgement}{The work was supported by the UK Engineering and Physical Sciences Research Council under grant number EP/P002757/1}

\seealso{
\code{\link{Reduction.Phase}}, \code{\link{Exploratory.Phase}}
}
\examples{
## Generates a random DGP using the DGP from the HCmodelSelection package
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=100, intercept=5, noise=1,
          var=1, d=1000, DGP.seed = 2019)

## Design matrix
dgp$X

## Response variable
dgp$Y

#Reduction Phase using only the first 70 observations
outcome.Reduction.Phase =  Reduction.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
                                           dmHC = 5, family=gaussian, seed.HC = 1093)

# Exploratory Phase using only the first 70 observations, choosing the variables which
# were selected at least two times in the third dimension reduction

idxs = outcome.Reduction.Phase$List.Selection$`Hypercube with dim 3`$numSelected2
outcome.Exploratory.Phase =  Exploratory.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
                                               list.reduction = idxs,
                                               family=gaussian, signif=0.01)


# Model Selection Phase using only the remainer observations
sq.terms = outcome.Exploratory.Phase$mat.select.SQ
it.terms = outcome.Exploratory.Phase$mat.select.INTER

outcome.Model.Selection.Phase = Model.Selection.Phase(X=dgp$X[71:100,],
                                                      Y=dgp$Y[71:100],  family = gaussian,
                                                      list.reduction = idxs,
                                                      square.terms = sq.terms,
                                                      interaction.terms = it.terms, signif=0.01)




}
