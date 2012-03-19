v <-
function(x,u, ...){
    # u
    if (!is.factor(u))
       stop("effect modifier must be nominal or ordinal \n")
    int <- rep(1, times = length(u))
    options(contrasts = c("contr.treatment", "contr.treatment"))
    dummies <- as.matrix(model.matrix(~ u)[,-1])
    # x
    if (is.factor(x)&&nlevels(x)>2)
       stop("varying coefficient not well defined. \n")
    options(na.action=na.pass)
    if (is.factor(x) && nlevels(x)==2)
       x <- -1*model.matrix(~ x, contrasts = list(x="contr.sum"))[,-1]
    # design
    design <- cbind( int-rowSums(dummies) , dummies) * as.vector(x)
    colnames(design) <- paste(".",levels(u), sep="")
    return(design)
}

