\name{wec.interact}
\alias{wec.interact}

\title{Function to create an interaction between two variables based on weighted effect coding.}
\description{This function facilitates the estimation of an interaction between two factor variables that are based on weighted effect coding. To that end, it creates a third variable that, together with the two original factor variables, forms the complete interaction. In interaction models, weighted effect coding displays the extra effect on top of the main effects found in a model without the interaction effect(s).}

\usage{
wec.interact(x1, x2, output.contrasts)
}


\arguments{
  \item{x1}{Factor variable (with contrasts based on weighted effect coding)}
  \item{x2}{Factor variable (with contrasts based on weighted effect coding) \emph{or} interval or ratio variable.}
  \item{output.contrasts}{Specifies whether the contrast matrix of the interaction should be returned. Defaults to FALSE, returning the model matrix. Option currently only implemented for interactions between one weighted effect coded and one interval or ratio variable.}}

\value{Returns a model matrix or contrast matrix for the interaction terms of (a.) two weighted effect coded variables, or (b.) one weighted effect coded and one interval or ratio variable.}

\references{
Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). When size matters: advantages of weighted effect coding in observational studies. \emph{International Journal of Public Health}, online access:http://link.springer.com/article/10.1007/s00038-016-0901-1

Grotenhuis, M. Te, Pelzer, B., Schmidt-Catran, A., Nieuwenhuis, R., Konig, R., and Eisinga, R. (2016). Weighted effect coded interactions: a novel moderation regression analysis for observational studies. \emph{International Journal of Public Health}, online access: http://link.springer.com/article/10.1007/s00038-016-0902-0

Sweeney, Robert E. and Ulveling, Edwin F. (1972) A Transformation for Simplifying the Interpretation of Coefficients of Binary Variables in Regression Analysis. \emph{The American Statistician}, 26(5): 30-32.
}

\author{Rense Nieuwenhuis, Manfred te Grotenhuis, Ben Pelzer, Alexander Schmidt, Ruben Konig, Rob Eisinga}

\note{
It should be noted that the procedure of applying weighted effect coding with interactions differs from the convential way to apply contrasts in R. This is becasue the contrast matrix of the interaction differs from the multiplication of the contrast matrix/matrices of the interacted variables.}

\seealso{
\code{\link{contr.wec}}
}
\examples{
data(BMI)

# Interaction two weighted effect coded categorical variables
BMI$childless.wec.yes <- BMI$childless.wec.no <- BMI$childless
contrasts(BMI$childless.wec.yes) <- contr.wec(BMI$childless, omitted="yes")
contrasts(BMI$childless.wec.no) <- contr.wec(BMI$childless, omitted="no")

BMI$age.wec.young <- BMI$age.wec.older <- BMI$age
contrasts(BMI$age.wec.young) <- contr.wec(BMI$age_categorical, omitted="Young (18-30)")
contrasts(BMI$age.wec.older) <- contr.wec(BMI$age_categorical, omitted="Older (60-70)")

model3a <- lm(BMI ~ childless.wec.yes + age.wec.young, data=BMI)
model3b <- lm(BMI ~ childless.wec.no + age.wec.older, data=BMI)
summary(model3a)
summary(model3b)


# Interaction
BMI$interact_c <- wec.interact(BMI$childless.wec.yes, BMI$age.wec.young)
BMI$interact_d <- wec.interact(BMI$childless.wec.yes, BMI$age.wec.older)
#BMI$interact_e <- wec.interact(BMI$childless.wec.no, BMI$age.wec.young)
#BMI$interact_f <- wec.interact(BMI$childless.wec.no, BMI$age.wec.older)

model3c <- lm(BMI ~ childless.wec.yes + age.wec.young + interact_c, data=BMI)
model3d <- lm(BMI ~ childless.wec.yes + age.wec.older + interact_d, data=BMI)
#model3e <- lm(BMI ~ childless.wec.no + age.wec.young + interact_e, data=BMI)
#model3f <- lm(BMI ~ childless.wec.no + age.wec.older + interact_f, data=BMI)

summary(model3c)
summary(model3d)
#summary(model3e)
#summary(model3f)


# Interaction weighted effect coded categorical variable and ratio/interval variable

data(PUMS)
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

\keyword{classes}
\keyword{category}
\keyword{manip}

