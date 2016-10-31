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
Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). When size matters: advantages of weighted effect coding in observational studies. \emph{International Journal of Public Health}, online access:http://link.springer.com/article/10.1007/s00038-016-0901-1

Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). Weighted effect coded interactions: a novel moderation regression analysis for observational studies. \emph{International Journal of Public Health}, online access: http://link.springer.com/article/10.1007/s00038-016-0902-0

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

BMI$childless.wec.yes <- contr.wec(BMI$childless, ref="yes")
BMI$childless.wec.no <- contr.wec(BMI$childless, ref="no")

BMI$age.wec.young <- contr.wec(BMI$age_categorical, ref="Young (18-30)")
BMI$age.wec.older <- contr.wec(BMI$age_categorical, ref="Older (60-70)")


model3a <- lm(BMI ~ childless.wec.yes + age.wec.young, data=BMI)
model3b <- lm(BMI ~ childless.wec.no + age.wec.older, data=BMI)
summary(model3a)
summary(model3b)


# Interaction
BMI$interact_c <- contr.wec.interact(BMI$childless.wec.yes, BMI$age.wec.young)
BMI$interact_d <- contr.wec.interact(BMI$childless.wec.yes, BMI$age.wec.older)
#BMI$interact_e <- contr.wec.interact(BMI$childless.wec.no, BMI$age.wec.young)
#BMI$interact_f <- contr.wec.interact(BMI$childless.wec.no, BMI$age.wec.older)

model3c <- lm(BMI ~ childless.wec.yes + age.wec.young + interact_c, data=BMI)
model3d <- lm(BMI ~ childless.wec.yes + age.wec.older + interact_d, data=BMI)
#model3e <- lm(BMI ~ childless.wec.no + age.wec.young + interact_e, data=BMI)
#model3f <- lm(BMI ~ childless.wec.no + age.wec.older + interact_f, data=BMI)

summary(model3c)
summary(model3d)
#summary(model3e)
#summary(model3f)

}

\keyword{classes}
\keyword{category}
\keyword{manip}

