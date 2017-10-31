wec.interact <- function (x1, x2, output.contrasts = FALSE)
{
  if (class(x1) != "factor") {
    stop("x1 needs to be a factor variable")
  }
  if (class(x2) == "factor") {
    mm <- model.matrix(~x1 * x2)
    mm <- mm[, -1]
    mm.full <- mm
    mm <- unique(mm)
    n.int <- (length(levels(x1)) - 1) * (length(levels(x2)) -
                                           1)
    n.main <- (length(levels(x1)) - 1) + (length(levels(x2)) -
                                            1)

    mm.int <- matrix(mm[, (dim(mm)[2] - n.int + 1):dim(mm)[2]], ncol=n.int)
    colnames(mm.int) <- colnames(mm)[-(1:n.main)]
    rownames(mm.int) <- rownames(mm)

    mm.main <- mm[, 1:n.main]
    ref.x1 <- paste(names(attr(table(x1, x2), "dimnames"))[1],
                    names(table(x1)), sep = "")
    ref.x2 <- paste(names(attr(table(x1, x2), "dimnames"))[2],
                    names(table(x2)), sep = "")
    reftable <- matrix(data = NA, nrow = length(table(x1)),
                       ncol = length(table(x2)))
    for (i in 1:dim(reftable)[1]) {
      for (j in 1:dim(reftable)[2]) {
        reftable[i, j] <- paste(ref.x1[i], ref.x2[j],
                                sep = ":")
      }
    }
    refcat.x1 <- setdiff(rownames(contrasts(x1)), colnames(contrasts(x1)))
    refcat.x2 <- setdiff(rownames(contrasts(x2)), colnames(contrasts(x2)))
    refcat.n <- table(x1, x2)[refcat.x1, refcat.x2]
    cond1 <- apply(mm.main < 0, 1, sum) == n.main
    values <- NA
    for (i in 1:ncol(mm.int)) {
      values[i] <- table(x1, x2)[which(reftable == colnames(mm.int)[i],
                                       arr.ind = TRUE)]/refcat.n
    }

    mm.int[cond1, ] <- rep(values, each = sum(cond1))
    cond2 <- which(mm.int < 0, arr.ind = TRUE)

    for (i in 1:nrow(cond2)) {
      num <- table(x1, x2)[which(reftable == colnames(mm.int)[cond2[i,
                                                                    2]])]
      denom <- table(x1, x2)[x1[as.numeric(rownames(cond2)[i])],
                             x2[as.numeric(rownames(cond2)[i])]]
      mm.int[cond2[i, 1], cond2[i, 2]] <- -1 * num/denom
    }

    mm[, -1:-n.main] <- mm.int
    if (output.contrasts) {
      warning("When interacting to weighted effect coded variables, coding matrix cannot be used to set contrasts.")
      ind <- as.numeric(row.names(mm.int))
      rownames(mm.int) <- paste(x1[ind], x2[ind], sep = ":")
      return(mm.int)
    }

    mm.full <- mm.full[, 1:n.main]
    mm.full <- as.data.frame(mm.full)
    mm <- as.data.frame(mm)
    output.data <- left_join(mm.full, mm, by = names(mm.full))
    output.data <- output.data[, -1:-n.main]
    output.data <- as.matrix(output.data)
    return(output.data)
  }
  if (class(x2) != "factor") {
    ref <- which(apply(contrasts(x1), 1, sum) != 1)
    weights <- tapply(x2, x1, function(x) {
      sum((x - mean(x))^2)
    })
    weights <- -1 * weights/weights[ref]
    effect.matrix <- contr.sum(length(levels(x1)))
    effect.matrix[-ref, ] <- effect.matrix[1:length(levels(x1)) -
                                             1, ]
    effect.matrix[ref, ] <- weights[-ref]
    if (output.contrasts) {
      return(effect.matrix)
    }
    colnames(effect.matrix) <- names(weights)[-ref]
    rownames(effect.matrix) <- names(weights)
    interact <- x1
    contrasts(interact) <- effect.matrix
    output <- model.matrix(~interact)[, -1] * x2
    output <- apply(output, 2, function(x) x - ave(x, x1))
    return(output)
  }
}

