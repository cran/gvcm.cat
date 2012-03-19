weight.function <-
function(
phi, 
control,
...
)

{
weight <- c()

if (control$assured.intercept==TRUE){
  if (control$index2[1]<0){
      weight <- c(weight,rep(phi, times=(0.5*control$index1[1]*(control$index1[1]-1))))
      } 
  if (control$index2[1]>0){weight <- c(weight,rep(phi, times=(control$index1[1]-1)))}  
  b <- 2 
} else {b <- 1} 

if(b<=length(control$index1)){
  for (i in b:length(control$index1)) {
  if (control$index2[i]<0 || control$index3[i]<0) {
    weight <- c(weight,rep(c(1-phi,phi), times=c(control$index1[i], 0.5*control$index1[i]*
             (control$index1[i]-1)))) 
    }  
  if (control$index2[i]>0) {
    weight <- c(weight,rep(c(1-phi,phi), times=c(control$index1[i], (control$index1[i]-1)))) 
    }  
  if (control$index3[i]>0) {
    if (control$p.ord.abs == TRUE){
    weight <- c(weight,rep(c(1-phi,phi), times=c(control$index1[i], (control$index1[i]-1)))) 
    } else {
    weight <- c(weight,rep(c(1-phi,phi), times=c(1, (control$index1[i]-1))))
    }
    }  
  }
}

if (control$adapted.weights == TRUE)  {  
weight <- weight * (1/abs(t(a.coefs(control$index1,control$index2,control$index3,
               control$assured.intercept,control$p.ord.abs))%*%control$oml))  
}

return(weight)

}

