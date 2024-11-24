\name{chirp_model_fit}
\alias{chirp_model_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{chirp_model_fit
%%  ~~function to do ... ~~
}
\description{
For real-life data, the exact number of Chirp components (p) needed for a good fit
is unknown. To estimate p, we use a sequential procedure: starting with a
one-component Chirp model, we fit the data using the "chirp_est" function, obtain
the residuals, and check their stationarity. We iteratively add components until the residuals become stationary, with k as the maximum number of components based on criteria like sample size. The p-component model is then fitted using "chirp_est" to obtain 4p parameter estimates. ~~
}
\usage{
chirp_model_fit(x)
}
%- To use this function, we need to provide the data and the maximum number of components (k) to consider, typically set to 15. This choice balances complexity and data size, as fitting a 15-component Chirp model requires estimating 60 parameters. For datasets with 100-200 observations, starting with k=15 or 16 is often optimal.


\arguments{
  \item{d}{A vector of time series observations}.
%%     ~~Describe \code{x} here~~
}
\item{k}{k is the maximum number of components you want to fit to the model, ideally p < k}.
%%     ~~Describe \code{x} here~~
}
}

\value{\item{rss_mean}{rss_mean is the mean of the squared residuals.}
       \item{kpss_test }{kpss_test returns the output of the KPSS test to check
       whether the residuals are stationary.}
       \item{adf_test}{adf_test returns the output of the Augmented Dickey-Fuller
       test to check whether the residuals are stationary.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{An article "Instructions for Creating Your Own R Package" by In Song
Kim, Phil Martin, Nina McMurry, Andy Halterman.
%% ~put references to the literature/web site here ~Dutta
}
\author{Tiyasa
%%  ~~who you are~~
}
\note{The same function can be applied to fit a multicomponent Elementary Chirp
model to a real life data taking α=0.
%%  ~~further notes~~
}

\examples{
d<-rnorm(100)
k=8
chirp_model_fit(d,k)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

