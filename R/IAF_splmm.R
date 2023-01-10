library(splmm)
library(tidyverse)
library(readr)

##########################################################################
# SP-LMM in IAF Metabolomics, CVD risk scores as response
# load metabolomics data
load("./data/iaf_metabolomics-cvd.RData")
x<-as.matrix(cbind(1,metabo[,c("Age","sex")],metabo[,5:253]))
colnames(x)[1]<-"Intercept"
y<-metabo$ASCVD
z<-as.matrix(x[,1:3])
grp<-as.factor(parse_number(metabo$person_id))

# find best lambda parameters
lam1 = seq(0.1,1,0.1)
lam2 = seq(0.1,1,0.1)
start<-Sys.time()
set.seed(17638191)
tuning<-splmmTuning(x=x,y=y,z=z,grp=grp,lam1=lam1,lam2=lam2,
                    nonpen.b = c(1,2,3), nonpen.L = c(1,2),
                    penalty.b="scad", penalty.L="scad",
                    control = splmmControl(maxArmijo = 50))
end<-Sys.time()
end-start

# model fit
set.seed(17638191)
fit<-splmm(x=x,y=y,z=z,grp=grp,nonpen.b = c(1,2,3), nonpen.L = c(1,2),
           lam1=tuning$min.lam1,lam2=tuning$min.lam2,
           penalty.b = "scad", penalty.L = 'scad',
           control = splmmControl(maxArmijo = 50))
summary(fit)


##########################################################################
# SP-LMM in IAF Proteomics, CVD risk scores as response
# load proteomics data
proteom<-read.csv("./data/iaf_proteomics-cvd.csv", sep=",", row.names = 1)
x<-as.matrix(cbind(1,proteom[,c("Age","sex")],proteom[10:275]))
colnames(x)[1]<-"Intercept"
y<-proteom$ASCVD
z<-as.matrix(x[,1:3,drop=FALSE])
grp<-as.factor(parse_number(proteom$person_id))

# find best lambda parameters
lam1 = seq(0.1,1,0.1)
lam2 = seq(0.1,1,0.1)
start<-Sys.time()
set.seed(17638191)
tuning<-splmmTuning(x=x,y=y,z=z,grp=grp,lam1=lam1,lam2=lam2,
                    nonpen.b = c(1,2,3), nonpen.L = c(1,2),
                    penalty.b="scad", penalty.L="scad",
                    control = splmmControl(maxArmijo = 50))
end<-Sys.time()
end-start

# model fit
set.seed(17638191)
fit <- splmm(x=x,y=y,z=z,grp=grp,lam1=tuning$min.lam1,lam2=tuning$min.lam2,
             nonpen.b = c(1,2,3), nonpen.L = c(1,2),
             penalty.b="scad", penalty.L="scad",
             control = splmmControl(maxArmijo = 50))
summary(fit)
