\name{DGP}
\alias{DGP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data generating process used by Battey, H. S. & Cox, D. R. (2018). }
\description{
This function generates realizations of random variables as described in the simple example of Battey, H. S. & Cox, D. R. (2018).
}
\usage{
DGP(s,a,sigStrength,rho,n,noise,var,d,intercept,DGP.seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{Number of signal variables.}
  \item{a}{Number of noise variables correlated with signal variables.}
  \item{sigStrength}{Signal strength}
  \item{rho}{Correlation among signal variables and noise variables correlated with signal variables.}
  \item{n}{Sample size.}
  \item{noise}{Variance of the observations around the true regression line.}
  \item{var}{Variance of the potential explanatory variables.}
  \item{d}{Number of potential explanatory variables.}
  \item{intercept}{Expected value of the response variable when all potential explanatory variables are at zero.}
  \item{DGP.seed}{Seed for the random number generator.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X}{The simulated design matrix.}
\item{Y}{The simulated response variable.}
%% ...
}
\references{
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generates DGP
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=100, intercept=5, noise=1,
          var=1, d=1000, DGP.seed = 2019)

## Design matrix
dgp$X

## Response variable
dgp$Y
}
