p <-
function(u, ...){
    old.contr <- getOption("contrasts")
    # u
    if (is.factor(u)){
     #  stop("effect 'u' in 'p(u)' must be nominal or ordinal \n")
    options(contrasts = c("contr.treatment", "contr.treatment")) # dummy
    design <- as.matrix(model.matrix(~ u)[,-1])
    # design
    colnames(design) <- paste(".",levels(u)[-1], sep="")
    } else {
    design <- matrix(model.matrix(~ u)[,-1],ncol=1)    
    }
    options(contrasts = old.contr)
    return(design)
}
