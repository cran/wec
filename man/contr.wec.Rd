\name{contr.wec}
\alias{contr.wec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function calculates contrasts for a factor variable based on weighted effect coding.}

\description{This function calculates contrasts for a factor variable based on weighted effect coding. In weighted effect coding the estimates from a first order regression model show the deviations per group from the sample mean. This is especially useful when a researcher has no directional hypotheses and uses a sample from a population in which the number of observations per group is different.}

\usage{
contr.wec(x, omitted)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Factor variable}
  \item{omitted}{Label of the factor label that should be taken as the omitted category}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{Returns a contrast matrix based on weighted effect coding.}

\references{
Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). When size matters: advantages of weighted effect coding in observational studies. \emph{International Journal of Public Health}, online access: http://link.springer.com/article/10.1007/s00038-016-0901-1

Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). Weighted effect coded interactions: a novel moderation regression analysis for observational studies. \emph{International Journal of Public Health}, online access: http://link.springer.com/article/10.1007/s00038-016-0902-0

Sweeney, Robert E. and Ulveling, Edwin F. (1972) A Transformation for Simplifying the Interpretation of Coefficients of Binary Variables in Regression Analysis. \emph{The American Statistician}, 26(5): 30-32.
}

\author{Rense Nieuwenhuis, Manfred te Grotenhuis, Ben Pelzer, Alexander Schmidt, Ruben Konig, Rob Eisinga}

\seealso{
\code{\link{wec.interact}}}

\examples{

data(BMI)


# Without controls
BMI$educ.wec.lowest <- BMI$educ.wec.highest <- BMI$educ

contrasts(BMI$educ.wec.lowest) <- contr.wec(BMI$education, omitted="lowest")
contrasts(BMI$educ.wec.highest) <- contr.wec(BMI$education, omitted="highest")

model.wec.lowest <- lm(BMI ~ educ.wec.lowest, data=BMI)
model.wec.highest <- lm(BMI ~ educ.wec.highest, data=BMI)
summary(model.wec.lowest)
summary(model.wec.highest)

# With Controls
BMI$sex.wec.female <- BMI$sex.wec.male <- BMI$sex
contrasts(BMI$sex.wec.female) <- contr.wec(BMI$sex, omitted="female")
contrasts(BMI$sex.wec.male) <- contr.wec(BMI$sex, omitted="male")

BMI$year.wec.2000 <- BMI$year.wec.2011 <- BMI$year
contrasts(BMI$year.wec.2000) <- contr.wec(BMI$year, omitted="2000")
contrasts(BMI$year.wec.2011) <- contr.wec(BMI$year, omitted="2011")

model.wec.lowest.controls <- lm(BMI ~ educ.wec.lowest +
  sex.wec.female + log_age + year.wec.2000,
  data=BMI)
model.wec.highest.controls <- lm(BMI ~ educ.wec.highest +
  sex.wec.male + log_age + year.wec.2011,
  data=BMI)
summary(model.wec.lowest.controls)
summary(model.wec.highest.controls)
}

\keyword{classes}
\keyword{category}
\keyword{manip}

