contr.wec <-
function(x, omitted)
{
    frequencies <- table(x)
    n.cat       <- length(table(x))
    omitted         <- which(levels(x) == omitted)

    new.contrasts <- contr.treatment(n.cat, base=omitted)
    new.contrasts[omitted,] <- -1 * frequencies[-omitted] / frequencies[omitted]

    colnames(new.contrasts) <- names(frequencies[-omitted])

    return(new.contrasts)
}
