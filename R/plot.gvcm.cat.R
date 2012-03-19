plot.gvcm.cat <-
function(

x,
accuracy=2,
type="path",

individual.paths=FALSE, # for type=="path"
xlim,
ylim,
main = NULL,

indent = 0,
color = TRUE,
...
)


{

# check 
  # x
  if (!("gvcm.cat" %in% is(x)))
       stop ("x must be a 'gvcm.cat' object. \n")
  # type
  if (!(type %in% c("path","score"))) 
       stop ("type is incorrect. \n")
     
# penalties...
if (x$method %in% c("nlm", "lqa")){
if (type=="path"){

    # check 
    if (is.matrix(x$plot[[1]])==FALSE)
    {stop ("Input argument 'plot' must be 'TRUE' for plotting coefficient paths. \n")}
    
    if (!is.logical(individual.paths) && !is.character(individual.paths) && !is.vector(individual.paths))
    {stop ("Error in input argument 'individual.paths'. \n")}

    if (length(as.vector(accuracy))!=1 || accuracy < 1 || accuracy > 4)
         stop ("accuracy must be a single integer > 0. \n")
         
    if (missing(xlim)) {xlim <- c(0,1)}
    if (!is.numeric(xlim) || !is.vector(xlim) || min(xlim)<0 || max(xlim)>1 || length(xlim)!=2)
         stop ("Error in arguemnt xlim. \n")  
         
    if (!is.logical(color))     
     {stop ("Error in input argument 'color'. \n")}

    if (indent<0 || !is.numeric(indent)) {indent <- 0}
    if (indent>0.8) {indent <- 0.8}

    
    # definitions
    lambda <- x$tuning[[1]]
    assured.intercept <- x$control$assured.intercept
    f <- x$formula
    number.selectable.parameters <- x$number.selectable.parameters
    index1 <- x$index1
    index2 <- x$index2 + x$index3
    control <- x$control
    n.p <- length(index1)

    # prepare path
    path <- x$plot[[1]]
    path[,-1] <- round(path[,-1], digits=accuracy)

    i <- 2  
    while (i<dim(path)[2] && dim(path)[1]-length(reduceBeta(path[,i],control$index1, control$index2, control$index3)) < number.selectable.parameters) {i<-i+1}
    path <- path[,1:i]
    path <- path[,order(as.numeric(colnames(path)),decreasing=TRUE)]
    
    # other definitions
    lambdas <- as.numeric(colnames(path))
    lambda.upper <- lambdas[1]
    x. <- 1- lambdas/lambda.upper
    s <- 1- (lambda/lambda.upper)
    if (color==TRUE){
    line.type <- - abs(rep(index2, times=index1)) + 2
    colour <- rep(1:n.p, times=index1)
    } else {
    line.type <- rep(rep(1:6, length.out=length(index1)),times=index1) 
    colour <- rep(gray( (0:(n.p-1)) / n.p),times=index1)  
    }
    
    # zoooom
    path <- path[,min(which(x.>=min(xlim))):max(which(x.<=max(xlim)))]
    x. <- x.[min(which(x.>=min(xlim))):max(which(x.<=max(xlim)))]

    # plot options
    lambda.max.print = TRUE
    lambda.opt.plot = TRUE
    lambda.opt.print = TRUE
    
    anteil.b <- 1/(3*22)
    range.x <- max(xlim)-min(xlim)
    coefs.left <- 4*anteil.b*range.x
    
    # remove NA's (e.g. if ML-estimate does not exisit)
    if (sum(is.na(path))>0) 
       warning("Some estimates contain NAs; set to zero for plotting. \n")
    for (k in 1:ncol(path)) path[which(is.na(path[,k])),k] <- 0 
     
    # individual.paths=TRUE    
    if (!is.logical(individual.paths) || individual.paths==TRUE) {
        lambda.max.print <- FALSE
        lambda.opt.print <- FALSE
         e <- c(1,cumsum(index1)[1:(n.p-1)]+1)
        index <- matrix(c(e, e+1, cumsum(index1)),nrow=n.p,ncol=3,byrow=FALSE)
        if (!is.logical(individual.paths)) {
            coefs.names <- c("Intercept")  
            tf <- terms(f,specials="v", "p")
            if(grepl("v\\(", as.character(attr(tf,"variables")[3]))) w <- 1 else w <- 0            
            for (i in (3+w):length(attr(tf,"variables"))) {
               if( grepl("p\\(", as.character(attr(tf,"variables")[i])) || grepl("v\\(", as.character(attr(tf,"variables")[i]))) {
                 coefs.names <- c(coefs.names,as.character(attr(tf, "variables")[[i]][2]))
               } else {
                 coefs.names <- c(coefs.names,as.character(attr(tf, "variables")[[i]]))
               }}
            which.coefs <- c()
            for (i in 1:length(coefs.names)) {
                   if (coefs.names[i] %in% individual.paths) {
                   which.coefs <- c(which.coefs, i)}
                 }
             } else {which.coefs <- 1:nrow(index)}
        if(length(which.coefs)<4){mfrows<-c(1,length(which.coefs))}else{mfrows<-c(2,ceiling(length(which.coefs)/2))}
        par(mfrow=mfrows)
        main <- ""# unlist(strsplit(as.character(f[3]),"\\+"))
        } else {
        index <- matrix(c(1, 2, sum(index1)),ncol=3, nrow=1)
        which.coefs <- 1:nrow(index)
        }
        
       
    for (j in which.coefs){# in 1:nrow(index)){
    
        if (!missing(ylim)) {
           if (!is.numeric(ylim) || !is.vector(ylim) || length(xlim)!=2)
               stop ("Error in arguemnt ylim. \n")  
           ylimj <- c(min(ylim), max(ylim))
           } else {
           ylimj <- c(floor(min(path[index[j,1]:index[j,3],])), 
           ceiling(max(path[index[j,1]:index[j,3],])))
           } 
        heigth <- ylimj[2] - ylimj[1]
    
        # legend
          legend.y <- path[index[j,1],dim(path)[2]]+0.05*heigth
          if (index[j,1] < index[j,3]) {
          for (i in index[j,2]:index[j,3]) {
              legend.y <- c(legend.y, path[i,dim(path)[2]]+0.05*heigth)
              }
          }
          # legend.x
            legend.x <- rep(1-coefs.left, times=length(legend.y))
            if (length(legend.y)>1) {
            legend.y.ordered <- legend.y[order(legend.y,decreasing=TRUE)]
            diff.legend <- (legend.y.ordered[1:(length(legend.y)-1)] - legend.y.ordered[2:(length(legend.y))])/heigth
            for (i in 1:length(diff.legend)){
            if (diff.legend[i] <0.02180){
            if(legend.x[order(legend.y,decreasing=TRUE)[i]] < 1-coefs.left+2*indent){
            legend.x[order(legend.y,decreasing=TRUE)[i+1]] <- legend.x[order(legend.y,decreasing=TRUE)[i]] + indent
            } else {
            legend.x[order(legend.y,decreasing=TRUE)[i+1]] <- legend.x[order(legend.y,decreasing=TRUE)[i]] - 2*indent            
            }
            }
            }
            }
          # correct margin
          anteil <- anteil.b * max(nchar(dimnames(x$coefficients)[[1]])) +  (indent!=0)*.9*(max(legend.x) - 1)*(1-indent) + (1-indent!=0)*.05
          right.margin <- anteil*range.x/(1-anteil)
          const. <- 0
          
        # plot
        matplot(x.,path[index[j,1],], type="l", lty=line.type[index[j,1]],    # lty=1 => solid
                col = colour[index[j,1]], lwd=2,                              # lty=2 => dashed
                main = main[j], axes = FALSE,                                 # lty=3 => punkte
                xlim = c(min(xlim),max(xlim) + right.margin),
                ylim = ylimj,
                xlab= expression(paste(1-lambda/lambda[max])), ylab= "estimated coefficient"
                )
        legend(x=legend.x[1]+const., y=legend.y[1],
               legend=rownames(path)[index[j,1]],box.lty=0, text.col=colour[index[j,1]])
        if (index[j,1] < index[j,3]) {
        for (i in index[j,2]:index[j,3]) {
                matlines(x.,path[i,], type="l", lty= line.type[i] ,col= colour[i],  lwd=2)
                legend(x=legend.x[i-index[j,2]+2]+const., y=legend.y[i-index[j,2]+2],legend=rownames(path)[i],box.lty=0, text.col=colour[i]) } }

        
        # x$lambda => dotted line
        if (lambda.opt.plot==TRUE){
            if(s>0){
            matlines(x=c(s, s), y=ylimj, type="l", lty=3, lwd=1)} else {s<-1} }
        if (lambda.opt.print==TRUE){
            legend("bottomleft", legend=c(expression(paste(lambda[CV], " = ")),(round(lambda, digits=2))) , box.lty=0, horiz=TRUE)
            }
        
        # max.lambda
        if (lambda.max.print==TRUE){
            legend("topleft", legend=c(expression(paste(lambda[max], " = ")),(round(lambda.upper, digits=2))) , box.lty=0, horiz=TRUE)
            }
            
        # the good look
        box()
        axis(1, at=seq(from=min(xlim),to=max(xlim),length.out=6), labels=as.character(seq(from=min(xlim),
                to=max(xlim),length.out=6)), col = "black", lty = 1, lwd = 1 )          
        axis(2, #at=seq(from = ylimj[1], to = ylimj[2], length.out=heigth),
                #labels=as.character(seq(from = ylimj[1], to = ylimj[2], length.out=heigth)),
                col = "black", lty = 1, lwd = 1 )    
        
        }

}
if (type=="score"){
    # check x$plot
    if (is.numeric(x$plot[[2]])==FALSE){stop ("type='score' requires cross-validation of lambda. \n")}
    
    # defintions
    score <- x$plot[[2]]
    ph <- as.numeric(rownames(score))
    phi <- x$tuning[[2]]
    lambda <- x$tuning[[1]]
    
    # phi fix
    if (is.matrix(score)==FALSE) {
        z <- score
        l <- as.numeric(names(score))
       } else {
        best <- which(rownames(score)==as.character(phi))
        rest <- which(rownames(score)!=as.character(phi))
        z <- score[best,]
        l <- as.numeric(colnames(score))
       }
       
    if(missing(xlim)){xlim <- c(floor(min(l)),ceiling(max(l)))}
    if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2)
         stop ("Error in arguemnt xlim. \n")  
    if(missing(ylim)){ylim <- c(floor(min(z)),ceiling(max(z)))}
    if (!is.numeric(ylim) || !is.vector(ylim) || length(ylim)!=2)
         stop ("Error in arguemnt ylim. \n")  
    
    matplot(l, z, type="l", lty=1, 
            col = 1,                         
            main = main, axes = FALSE, 
            xlim = xlim, 
            ylim = ylim,
            xlab= expression(paste(lambda)), ylab= "score"
            )
    
    matlines(x=c(lambda, lambda), y=ylim, type="l", lty=3, lwd=1)
    
    # phi flexible => add rest    
    if (length(ph)>1) {
        j <- 2
        for (i in rest) {
            z <- score[i,]
            matlines(l, z, type="l", lty= 1 ,col=j)
            j <- j+1
           }
        label <- factor(1:length(ph), levels= 1:length(ph), labels = as.character(ph[c(best,rest)]) )
        legend("bottomright", levels(label), fill=1:length(ph))    
       }
       
    # good look
    box()
    axis(1, at=seq(from = floor(min(l)), to = ceiling(max(l)), by = 1),
            labels=as.character(seq(from = min(xlim), to = max(xlim), by = 1)),
            col = "black", lty = 1, lwd = 1 )    
    axis(2, #at=seq(from = floor(min(z)), to = ceiling(max(z)), by = 1),
            #labels=as.character(seq(from = min(ylim), to = max(ylim), by = 1)),
            col = "black", lty = 1, lwd = 1 )    
    
}
}

# forward selection
if (x$method %in% c("AIC", "BIC")){}


}

