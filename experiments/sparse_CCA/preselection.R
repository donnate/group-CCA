

C = cor(example$Data)
diag(example$Sigma)
C = example$Sigma
color  =rep(0, nrow(example$Sigma))
color[which(apply(example$u,1, norm)>0)]=1
color[nrow(example$u) + which(apply(example$v,1, norm)>0)]=1
t = apply(C-diag(diag(C)), 1, function(x){max(x^2)})
ggplot(data.frame(y = t,
                  x = 1:length(t),
                  c = color), aes(x=x, y=y,colour=c)) + 
  geom_point()
indices <- order(-t)

J = order(-t)[1: ceiling(0.8  * n/4)]
set_u = J[which(J <= p1)]
set_v = J[which(J > p1)]
#### We take these and we hope that they are enough essentially

r=2
t=CCA::cc(as.matrix(example$Data[,set_u]), as.matrix(example$Data[, set_v]))
Uhat = matrix(0, p, r)
Uhat[set_u, ] =  t$xcoef[,1:r]
Uhat[set_v, ] =  t$ycoef[,1:r]
subdistance(Uhat, example$a)



preselection <-function(Data, CorrelationMat, p1, r, alpha){
  p = ncol(Data)
  p2 =  p-p1
  
  t = apply(CorrelationMat -diag(diag(CorrelationMat)), 1, 
            function(x){max(x^2)})
  J = order(-t)[1: ceiling(alpha  * n/(2*r))]
  set_u = J[which(J <= p1)]
  set_v = J[which(J > p1)]
  t=CCA::cc(as.matrix(Data[,set_u]), as.matrix(Data[, set_v]))
  Uhat = matrix(0, p, r)
  Uhat[set_u, ] =  t$xcoef[,1:r]
  Uhat[set_v, ] =  t$ycoef[,1:r]
  subdistance(Uhat, example$a)
  
}