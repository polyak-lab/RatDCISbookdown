# Function to bootstrap shannon idx values.

BootstrapShannonIdx=function(X.counts, N){
  xsum=sum(X.counts)
  XFrac=X.counts/xsum
  replVal=sapply(1:N, function(x) sample(c(1:length(X.counts)), xsum, replace = T, prob=XFrac))
  replFrac=lapply(1:N, function(x) table(replVal[ ,x])/xsum)
  Sidx=sapply(1:N, function(x) -sum(replFrac[[x]]*log(replFrac[[x]]), na.rm = T))
  CI=quantile(Sidx, c(0.025, 0.975))
  CI
}


PermuteGini=function(X.counts, N){
  xsum=sum(X.counts)
  XFrac=X.counts/xsum
  replVal=sapply(1:N, function(x) sample(c(1:length(X.counts)), xsum, replace = T, prob=XFrac))
  replFrac=lapply(1:N, function(x) table(replVal[ ,x])/xsum)
  Gini=sapply(replFrac, function(x) gini(x))
  CI =quantile(Gini, c(0.025, 0.975))
  CI
}