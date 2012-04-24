weight.function <-
function(
phi,
indices,
oml, 
control,
acoefs,
...
)

{
if (missing(acoefs))
    acoefs <- a.coefs(indices, control, oml)

weight <- c()

index1 <- indices[[1]]
index2 <- indices[[2]]
index3 <- indices[[3]]

assured.intercept <- control$assured.intercept
adapted.weights <- control$adapted.weights

if (assured.intercept==TRUE){
  if (index2[1]<0){
      weight <- c(weight,rep(phi, times=(0.5*index1[1]*(index1[1]-1))))
      } 
  if (index2[1]>0){weight <- c(weight,rep(phi, times=(index1[1]-1)))}  
  b <- 2 
} else {b <- 1} 

# phi <- .5 # Sarah!!
if(b<=length(index1)){
  for (i in b:length(index1)) { 

  # index2
  if (index2[i]<0) {
    weight <- c(weight,rep(c(1-phi,phi), times=c(index1[i], 0.5*index1[i]*
             (index1[i]-1)))) 
    }  
  if (index2[i]>0) {
    weight <- c(weight,rep(c(1-phi,phi), times=c(index1[i], (index1[i]-1)))) 
    }  
  
  # index3
  if (index3[i]<0) {
    weight <- c(weight,rep(c(1,1), times=c(index1[i], 0.5*index1[i]*(index1[i]-1)))) 
    }  
  
  if (index3[i]>0) {
    weight <- c(weight,rep(c(1,1), times=c(1, (index1[i]-1))))
    }  
  }
}

if (adapted.weights == TRUE)  {  
weight <- weight * (1/abs(t(acoefs)%*%oml))  
}

return(weight)

}

