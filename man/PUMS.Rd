\name{PUMS}
\alias{PUMS}
\docType{data}
\title{Public Use Microdata Sample files (PUMS) 2013}
\description{
  The ACS Public Use Microdata Sample files (\code{PUMS}) are a sample of the actual responses to the American Community Survey and include most population and housing characteristics. }

\format{
  A data frame with 10000 observations on the following 4 variables.
  \describe{
    \item{\code{wage}}{annual wages (binned to 1000s, top-coded, in US dollar)}
    \item{\code{race}}{a factor with levels \code{Hispanic}, \code{Black}, \code{Asian}, and \code{White}}
    \item{\code{education.int}}{level of education}
    \item{\code{education.cat}}{a factor variable with levels \code{High School}, and \code{Degree}}
  }
}

\source{
These data are a random subset of 10000 observations from working individuals aged over 25 in the 2013 ACS Public Use Microdata Sample files (\code{PUMS}).}

\examples{

data(PUMS)

PUMS$race.wec <- factor(PUMS$race)
contrasts(PUMS$race.wec) <- contr.wec(PUMS$race.wec, "White")
contrasts(PUMS$race.wec)

m.wec <- lm(wage ~ race.wec, data=PUMS)
summary(m.wec)

PUMS$race.educint <- wec.interact(PUMS$race.wec, PUMS$education.int)
m.wec.educ <- lm(wage ~ race.wec + education.int + race.educint, data=PUMS)
summary(m.wec.educ)

}
\keyword{datasets}
