\name{BMI}
\alias{BMI}
\docType{data}
\title{Data on BMI of Dutch citizens}
\description{
  The \code{BMI} data contains information on Dutch individuals' BMI, in addition to select socio-demographic variables.}

\format{
  A data frame with 3323 observations on the following 6 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{male} and \code{female}}
    \item{\code{education}}{a factor with levels \code{lowest}, \code{middle}, and \code{highest}}
    \item{\code{year}}{a factor with levels \code{2000}, \code{2005}, and \code{2011}}
    \item{\code{BMI}}{interval variable representing respondents' Body Mass Index (BMI)}
    \item{\code{children}}{a factor with levels \code{Children} and \code{Childless}}
    \item{\code{log_age}}{interval variable representing the natural log of respondents'age}
  }
}

\source{
These data are a subset from three waves of the `Socio-Cultural Developments in the Netherlands' (SOCON) datasets, collected at the Radboud University in the Netherlands (see references for original codebooks).}

\references{
Eisinga, R., G., Kraaykamp, P. Scheepers, P. Thijs (2012). \emph{Religion in Dutch society 2011-2012. Documentation of a national survey on religious and secular attitudes and behaviour in 2011-2012}, DANS Data Guide 11, The Hague: DANS/Pallas Publications  Amsterdam University Press, 184p.

Eisinga, R., A. Need, M. Coenders, N.D. de Graaf, M. Lubbers, P. Scheepers, M. Levels, P. Thijs (2012).
\emph{Religion in Dutch society 2005. Documentation of a national survey on religious and secular attitudes and behaviour in 2005}, DANS Data Guide 10, The Hague: DANS/Pallas Publications  Amsterdam University Press, 246p.

Eisinga, R., M. Coenders, A. Felling, M. te Grotenhuis, S. Oomens, P. Scheepers (2002). \emph{Religion in Dutch society 2000. Documentation of a national survey on religious and secular attitudes in 2000}, Amsterdam: NIWI-Steinmetz Archive, 374p.

}
\examples{

data(BMI)

# Without Controls
model.dummy <- lm(BMI ~ education, data=BMI)
summary(model.dummy)

# With Controls
model.dummy.controls <- lm(BMI ~ education + sex + log_age + children + year, data=BMI)
summary(model.dummy.controls)

}
\keyword{datasets}
