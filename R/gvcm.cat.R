gvcm.cat <-
function(
formula,
data,
family = gaussian,
method = c("lqa", "AIC", "BIC"),
tuning = list(lambda=TRUE, specific=FALSE, phi=0.5, grouped.fused=0.5, elastic=0.5, vs=0.5, spl=0.5),
weights,
offset,
start,
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

