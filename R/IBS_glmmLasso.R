library(glmmLasso)

load("./data/ibs_all.RData")
# log-transform
for (i in 1:nrow(ibs)) {
  for (j in 7:ncol(ibs)) {
    if(ibs[i,j]!=0){
      ibs[i,j]=log(ibs[i,j], base = 10)
    }else(ibs[i,j]=1)
  }
}
# include necessary variables and rename
excl<-c("time","Group")
ibs<-ibs[,!colnames(ibs) %in% excl]
colnames(ibs)[3:28]<-c(paste0("x", seq(1,26)))
xnam <- paste("x", 1:26, sep="")

# find best lambda parameter
res<-as.data.frame(matrix(NA,nrow=length(lseq), ncol=3))
colnames(res)<-c("l","bic","iter")
start<-Sys.time()
lseq<-seq(0,300); i<-1
for (l in lseq) {
  lm2 <- glmmLasso(reformulate(xnam, "status"),
                   family = binomial(link="logit"), rnd = list(subject=~1+x1), lambda=l,
                   data = ibs, control=list(print.iter=TRUE))
  res[i,]<-c(l,lm2$bic,lm2$conv.step)
  #save(res,file = paste0("restmp_i_",i,"_lambda_",l,".Rdata"))
  i<-i+1
  print(l)
  print(lm2$coefficients)
}
end<-Sys.time()
runtime<-end-start

# BIC window between 0-30
plot(res$l[res$l<=30],res$bic[res$l<=30])
res2<-res[res$l<=30,]
min.b<-min(res2$bic, na.rm = T)
l=res[res$bic==min.b,]$l # final l=29
lm2 <- glmmLasso(reformulate(xnam, "status"), 
                 family = binomial(link="logit"), rnd = list(subject=~1+x1), lambda=l, 
                 data = ibs, control=list(print.iter=TRUE))
reslist<-list(
  lm2$call,
  summary(lm2),
  lm2$aic,
  lm2$bic,
  lm2$loglik,
  lm2$coefficients,
  lm2$StdDev,
  lm2$lambda.max
)