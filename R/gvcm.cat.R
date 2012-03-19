gvcm.cat <-
function(
formula,
data,
family = gaussian,
method = "lqa", 
tuning = list(lambda=TRUE, phi=0.5),
weights,
control,
model = FALSE,
x = FALSE,
y = FALSE,
plot=FALSE,
...
)
{
UseMethod("gvcm.cat")
}

