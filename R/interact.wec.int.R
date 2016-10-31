interact.wec.int <- function(model, x.int, x.wec, ref)
{

  updated.data <- model.frame(model)
  # TAke a model and two specified variables (1 weighted effect coding, 1 interval)

  #get predicted scores of the model

  # Set correct reference
  #update the model
  updated.data$x.wec <- relevel(updated.data$x.wec, ref=ref)
  model.update1 <- update(model, data=updated.data)
  #update the model to include the interaction


  # Center the interval variable
  updated.data$x.int <- updated.data$x.int - mean(updated.data$x.int)
  # apply WEC on categorical variable

  updated.data$x.wec <- contr.wec(updated.data$x.wec, ref)

  update(model, data=updated.data)
}
