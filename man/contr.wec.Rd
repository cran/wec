\name{contr.wec}
\alias{contr.wec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to create a factor variable with contrasts based on weighted effect coding.}
\description{This function provides contrasts to a factor variable that are based on weighted effect coding. In weighted effect coding the estimates from a first order regression model show the deviations per group from the sample mean. This is especially useful when a researcher has no directional hypotheses and uses a sample from a population in which the number of observations per group is different.}
\usage{
contr.wec(x, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Factor variable}
  \item{ref}{Label of the factor label that should be taken as the reference category}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{Returns a factor variable (with weighted effect coding).}

\references{
Sweeney, Robert E. and Ulveling, Edwin F. (1972) A Transformation for Simplifying the Interpretation of Coefficients of Binary Variables in Regression Analysis. \emph{The American Statistician}, 26(5): 30-32.
}

\author{Rense Nieuwenhuis, Manfred te Grotenhuis, Ben Pelzer, Alexanter Schmidt, Ruben Konig, Rob Eisinga}

\note{It should be noted that the procedure of applying weighted effect coding differs from the convential way to apply contrasts in R. This is because to apply weighted effect coding, unlike with for example treatment coding, information is required on the sample mean and the distribution of the factor categories. Hence, instead of specifying a type of contrasts when creating the factor variable (or in the specification of the regression analysis) a new factor variable with weighted effect coding is created. It is this new factor variable that is then to be used in the regression analyses.}

\seealso{
\code{\link{contr.wec.interact}}}

\examples{

data(BMI)

#
# Treatment Coding (Default coding in R)
#

# Without Controls
model.dummy <- lm(BMI ~ education, data=BMI)
summary(model.dummy)

# With Controls
model.dummy.controls <- lm(BMI ~ education + sex + log_age + children + year, data=BMI)
summary(model.dummy.controls)

#
# Weighted Effect Coding
#

# Without controls
BMI$educ.wec.lowest <- contr.wec(BMI$education, ref="lowest")
BMI$educ.wec.highest <- contr.wec(BMI$education, ref="highest")
model.wec.lowest <- lm(BMI ~ educ.wec.lowest, data=BMI)
summary(model.wec.lowest)
model.wec.highest <- lm(BMI ~ educ.wec.highest, data=BMI)
summary(model.wec.highest)

# With Controls
BMI$sex.wec.female <- contr.wec(BMI$sex, ref="female")
BMI$sex.wec.male <- contr.wec(BMI$sex, ref="male")
BMI$year.wec.2000 <- contr.wec(BMI$year, ref="2000")
BMI$year.wec.2011 <- contr.wec(BMI$year, ref="2011")
BMI$children.wec <- contr.wec(BMI$children, ref="Children")
BMI$children.wec.childless <- contr.wec(BMI$children, ref="Childless")

model.wec.lowest.controls <- lm(BMI ~
  educ.wec.lowest +  sex.wec.female + log_age + children.wec + year.wec.2000,
  data=BMI)
summary(model.wec.lowest.controls)
model.wec.highest.controls <- lm(BMI ~
  educ.wec.highest +  sex.wec.male + log_age + children.wec.childless + year.wec.2011,
  data=BMI)
summary(model.wec.highest.controls)
}

\keyword{classes}
\keyword{category}
\keyword{manip}

