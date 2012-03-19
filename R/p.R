p <-
function(u, ...){
    # u
    if (!is.factor(u))
       stop("effect 'u' in 'p(u)' must be nominal or ordinal \n")
    options(contrasts = c("contr.treatment", "contr.treatment"))
    design <- as.matrix(model.matrix(~ u)[,-1])
    # design
    colnames(design) <- paste(".",levels(u)[-1], sep="")
    return(design)
}
