index <-
function(
dsgn, data=data
)

{
# definitions
  Terms <- dsgn$Terms
  label <- attr(Terms, "term.labels")
  variables <- as.character(attr(Terms,"variables"))[-c(1)]
  m <- dsgn$m
  int <- dsgn$int
  X <- model.matrix(Terms, m)
  ass <- attr(X,"assign")

# index1
  index1 <- summary(as.factor(ass))
    if (int==0) index1 <- index1[-1]
  names(index1) <- NULL

# index2-4
  index2 <- index3 <- rep(0,length(index1))
  index4 <- rep(1,length(index1))
  j <- 1 # label
  r <- 1

  # index2
  if(!is.null(attr(Terms,"specials")$v)){
      # intercept
      if (index1[1]!=1) {
        u <- model.frame(~eval(parse(text=gsub(")", "", gsub(" ", "",
             strsplit(names(m)[2], "\\,")[[1]][2])))), data)
        erms <- attr(u, "terms")
        if (attr(erms, "dataClasses")[[1]]=="factor") {index2[1] <- c(-1)}
        if (attr(erms, "dataClasses")[[1]]=="ordered"){index2[1] <- c(+1)}
        j <- 2
        r <- 0
      }

      # others
      if (length(label)>=j) {
      for (i in j:length(label)) {
        variable <- which(variables==label[i])
        if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$v){
          u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
               strsplit(label[i], "\\,")[[1]][2])))), data)
          erms <- attr(u, "terms")
          if (attr(erms, "dataClasses")[[1]]=="factor") {index2[i+r] <- -1}
          if (attr(erms, "dataClasses")[[1]]=="ordered"){index2[i+r] <- +1}
          }}
        } # for
      } # uberhaupt others
  } # uberhaupt


# index3 + 4
  if(!is.null(attr(Terms,"specials")$p)){
      if (length(label)>=j) {
      for (i in j:length(label)) {
        variable <- which(variables==label[i])
        if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$p){
          u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
               strsplit(label[i], "p\\(")[[1]][2])))), data)
          erms <- attr(u, "terms")
          if (attr(erms, "dataClasses")[[1]]=="factor") {index3[i+r] <- -1}
          if (attr(erms, "dataClasses")[[1]]=="ordered"){index3[i+r] <- +1}
          if (attr(erms, "dataClasses")[[1]]=="numeric"){index3[i+r] <- +1}
          }}

       if (length(variable)==0 && grepl(":", label[i]) && !grepl("v\\(", label[i])) { # interaktion ohne v
          eins <- strsplit(label[i], ":")[[1]][1]
          zwei <- strsplit(label[i], ":")[[1]][2]
          for (m in c(eins, zwei)){
          welche <- which(variables==m)
            if (length(welche)>0) {if (welche %in% attr(Terms,"specials")$p){
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(variables[welche], "p\\(")[[1]][2])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index3[i+r] <- -1}
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index3[i+r] <- +1}
              if (attr(erms, "dataClasses")[[1]]=="numeric"){index3[i+r] <- +1}

              if (variables[welche] %in% label && attr(erms, "dataClasses")[[1]]!="numeric") # nur int mit p(cat) in index4 schreiben!!
                 {index4[i+r] <- which(label==variables[welche])+r}
              }}
          }}

        } # for
      } # uberhaupt others
  } # uberhaupt
  

# return
# return(list(index1, index2, index3, index4))
return(list(index1, index2, index3))

}
