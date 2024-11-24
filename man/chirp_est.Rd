\name{chirp_est}
\alias{chirp_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{chirp_est
%%  ~~function to do ... ~~
}
\description{
The Chirp model consists of two non-linear parameters: α (non zero unknown
frequency present in the model) and β (unknown chirp rate present in the model)
and two linear parameters A and B which are the amplitudes. Now this function takes the input data and generates grid for the non-linear parameters over the interval (0,pi) X (0,pi). The Approximate Least Squares approach is being used, in this approach, the periodogram function is maximized iteratively to update the estimates of α and𝛽.The process continues until a convergence criterion is satisfied, yielding the Approximate Least Squares estimates.
}
\usage{
chirp_est(x)
}
%- For using this function we need to give the data or a vector of time series observations, the function will return the estimates of the four parameters of the Chirp model and also the estimate of the standard error.
\arguments{
  \item{y_t}{A vector of time series observations}.
%%     ~~Describe \code{x} here~~
}
}
\value{\item{est}{est is a vector which returns the initial choices of α,β, the
     approximate least square estimates of α, β, A, B and also the estimate of
     the standard error.}
     \item{res}{res is the residuals obtained by fitting the model.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{An article "Instructions for Creating Your Own R Package" by In Song
Kim, Phil Martin, Nina McMurry, Andy Halterman.
%% ~put references to the literature/web site here ~
}
\author{Tiyasa Dutta
%%  ~~who you are~~
}
\note{The above discussion is applicable for Elementary Chirp model taking α=0.
%%  ~~further notes~~
}



\examples{
y_t<-rnorm(100)
chirp_est(y_t)
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
