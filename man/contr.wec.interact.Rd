\name{contr.wec.interact}
\alias{contr.wec.interact}

\title{Function to create an interaction between two factor variables based on weighted effect coding.}
\description{This function facilitates the estimation of an interaction between two factor variables that are based on weighted effect coding. To that end, it creates a third variable that, together with the two original factor variables, forms the complete interaction. In second order (interaction) regression models the weighted effect coded interactions reflect the additional deviation from the sample mean on top of the main effects. This is useful in case one has an unbalanced design, i.e., the number of observations varies across the interaction terms.}
\usage{
contr.wec.interact(x1, x2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{Factor variable 1 (with contrasts based on weighted effect coding)}
  \item{x2}{Factor variable 2 (with contrasts based on weighted effect coding)}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{Returns a factor variable (with weighted effect coding).}

\references{
Sweeney, Robert E. and Ulveling, Edwin F. (1972) A Transformation for Simplifying the Interpretation of Coefficients of Binary Variables in Regression Analysis. \emph{The American Statistician}, 26(5): 30-32.
}

\author{Rense Nieuwenhuis, Manfred te Grotenhuis, Ben Pelzer, Alexanter Schmidt, Ruben Konig, Rob Eisinga}

\note{
It should be noted that the procedure of applying weighted effect coding differs from the convential way to apply contrasts in R. This is because to apply weighted effect coding, unlike with for example treatment coding, information is required on the sample mean and the distribution of the factor categories. This also applies to the \emph{interaction} between two factors with weighted effect coding.}

\seealso{
\code{\link{contr.wec}}
}
\examples{
data(BMI)

###
# Model I: No Interaction
###

BMI$educ.wec.lowest <- contr.wec(BMI$education, ref="lowest")
BMI$sex.wec.female <- contr.wec(BMI$sex, ref="female")
model1a <- lm(BMI ~ educ.wec.lowest +  sex.wec.female, data=BMI)
summary(model1a)

# To obtain estimates for all categories, also run:
#BMI$educ.wec.highest <- contr.wec(BMI$education, ref="highest")
#BMI$sex.wec.male <- contr.wec(BMI$sex, ref="male")
#model1b <- lm(BMI ~ educ.wec.highest +  sex.wec.male, data=BMI)
#summary(model1b)

###
# Model II: Interaction sex * education
###

BMI$int.2a <- contr.wec.interact(BMI$educ.wec.lowest, BMI$sex.wec.female)
model2a <- lm(BMI ~ educ.wec.lowest +  sex.wec.female + int.2a, data=BMI)
summary(model2a)

# To obtain estimates for all categories, also run:
#BMI$int.2b <- contr.wec.interact(BMI$educ.wec.highest, BMI$sex.wec.female)
#BMI$int.2c <- contr.wec.interact(BMI$educ.wec.lowest, BMI$sex.wec.male)
#BMI$int.2d <- contr.wec.interact(BMI$educ.wec.highest, BMI$sex.wec.male)

#model2b <- lm(BMI ~ educ.wec.highest +  sex.wec.female + int.2b, data=BMI)
#model2c <- lm(BMI ~ educ.wec.lowest +  sex.wec.male + int.2c, data=BMI)
#model2d <- lm(BMI ~ educ.wec.highest +  sex.wec.male + int.2d, data=BMI)

#summary(model2b)
#summary(model2c)
#summary(model2d)
}

\keyword{classes}
\keyword{category}
\keyword{manip}

