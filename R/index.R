index <-
function(formula, data)
{
# definitions
    special <- c("v", "p")
    int <- if (grepl("v\\(", strsplit(deparse(formula[3]), "\\+")[[1]][1]) ) 
        0 else 1
    m <- model.frame(formula=terms(formula, specials=special, data=data), data)
    tf <- attr(m, "terms")
    attr(tf, "intercept") <- 1
#    options(contrasts = c("contr.effect", "contr.effect"))
    X <- model.matrix(tf, m)
    if (int==0) X <- X[,-1]

# index1
  # intercept
  if (int==0) {index1 <- ncol(as.matrix(m[2])); j <- 3} else 
     {index1 <- 1; j <- 2}
  
  # others
  if (j <= ncol(m))  {
  for (i in j:ncol(m)){
    if (i %in% attr(tf,"specials")$v || i %in% attr(tf,"specials")$p) 
        {index1 <- c(index1, ncol(as.matrix(m[i])))}

    if (attr(tf,"dataClasses")[[i]] == "numeric")
        {index1 <- c(index1, 1)}

    if (attr(tf,"dataClasses")[[i]] %in% c("factor", "ordered"))
        {index1 <- c(index1, nlevels(as.factor(as.matrix(m[i])))-1)}
  }
  }
  
  # consistency
  if (dim(X)[2] != sum(index1))
    stop ("inconsistency concerning 'index1'! \n")

# index2
  # intercept
  if (j == 2) {index2 <- c(0)}
  if (j == 3) {
    u <- model.frame(~eval(parse(text=gsub(")", "", gsub(" ", "", 
         strsplit(names(m)[2], "\\,")[[1]][2])))), data)
    erms <- attr(u, "terms")
    if (attr(erms, "dataClasses")[[1]]=="factor") {index2 <- c(-1)}
    if (attr(erms, "dataClasses")[[1]]=="ordered"){index2 <- c(+1)}
    }
    
  # others
  if (j <= ncol(m))  {
  for (i in j:ncol(m)){
    if (i %in% attr(tf,"specials")$v){ 
    u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "", 
         strsplit(names(m)[i], "\\,")[[1]][2])))), data)
    erms <- attr(u, "terms")
    if (attr(erms, "dataClasses")[[1]]=="factor") {index2 <- c(index2, -1)}
    if (attr(erms, "dataClasses")[[1]]=="ordered"){index2 <- c(index2, +1)}
    } else {index2 <- c(index2, 0)}
  }
  }
  if (length(index1) != length(index2))
    stop ("inconsistency concerning 'index1' and 'index2'! \n")
    
# index3
  # intercept
  index3 <- 0
  
  # others
  if (2 <= ncol(m))  {
  for (i in j:ncol(m)){
    if (i %in% attr(tf,"specials")$p){ 
    u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "", 
         strsplit(names(m)[i], "p\\(")[[1]][2])))), data)
    erms <- attr(u, "terms")
    if (attr(erms, "dataClasses")[[1]]=="factor") {index3 <- c(index3, -1)}
    if (attr(erms, "dataClasses")[[1]]=="ordered"){index3 <- c(index3, +1)}
    if (attr(erms, "dataClasses")[[1]]=="numeric"){index3 <- c(index3, +1)}
    } else {index3 <- c(index3, 0)}
  }
  }
  if (length(index1) != length(index3))
    stop ("inconsistency concerning 'index1' and 'index3'! \n")
    

# return
return(list(index1, index2, index3))

}

