library(dplyr)
source("./main/fun_volcano.R")

load("./data/ibs_all_h_conf.RData") #healthy
load("./data/ibs_all_c_conf.RData") #ibs-c
load("./data/ibs_all_d_conf.RData") #ibs-d

catch.err.lqmm<-function(db,alpha,lambda.u.seq,lambda.z.seq){
  out<-tryCatch(res<-jqm(db=db,
                         alpha=alpha,
                         lambda.u.seq=lambda.u.seq,
                         lambda.z.seq=lambda.z.seq), 
                error=function(e){
                  message("catch error")
                  print(e)
                })
  out<-inherits(out, what="error")
  if(out==FALSE){
    res2<-res
  }else{res2<-NA}
  return(list(out=out,
              res=res2))
}

iri_fun<-function(db,df2,alpha){
  d1 <- trendsub(db=db, lim=2.3)
  d2 <- varboot(db=db, lim=2.3)
  Subject<-unique(c(d1$exc_sub, d2$out.mad))
  df2<-df2[!df2$subject %in% Subject,]
  colnames(df2)<-c("subject","time","age","sex","y")
  
  exmad<-d1$d_mad
  inc<-exmad[!exmad$subject %in% Subject,]
  df2<-df2 %>% left_join(., inc, by="subject")
  # define another outlier: 5MAD and exclude them
  df2$outlier<-ifelse(df2$y<df2$mad_loww | df2$y>df2$mad_upp,3,
                      ifelse(df2$y<df2$mad_low | df2$y>df2$mad_up,1,0))
  df2$outlier<-ifelse(is.na(df2$outlier)==T,0,df2$outlier)
  df2<-df2[df2$outlier!=3,]
  db<-df2
  subjects<-unique(db$subject)
  cnt<-1; db2<-0; db3<-0
  for(s in subjects) {
    db.y<-db[(db$subject==s),]
    #db.y$outlier<-ifelse(db.y$time==max(db.y$time),2,db.y$outlier) # uncoment this if we want to exclude the last observation
    dba<-db.y
    db.y<-db.y[(db.y$time!=max(db.y$time)),] # exclude the last observation in the IRI estimation
    
    db2<-rbind(db2,dba)
    db3<-rbind(db3,db.y)
  }
  df2<-db2[-1,]; df3<-db3[-1,]
  age<-NULL;sex<-NULL;cnt<-1
  for (s in subjects) {
    df.s<-df2[df2$subject==s,]
    age[cnt]<-df.s$age[1]
    sex[cnt]<-df.s$sex[1]
    cnt<-cnt+1
  }
  source("./main/JQM_Function_FE.R")
  err<-catch.err.lqmm(db=df2,
                      alpha=alpha,
                      lambda.u.seq = seq(0.02,0.1,0.02),
                      lambda.z.seq = seq(0.5,5,0.5))
  res<-err$res
  if(err$out==FALSE){
    uz <- cbind.data.frame(res$u, res$z)
    uz$low <- res$beta0 + res$u + res$z*res$beta1 + res$beta3*age + res$beta4*sex
    uz$up  <- res$beta0 + res$u + res$z*res$beta2 + res$beta3*age + res$beta4*sex
    uz$id <- unique(df2$subject)
    colnames(uz) <- c("u", "z", "low", "up","id")
  }else{
    source("./main/JQM_Function_age.R")
    err2<-catch.err.lqmm(db=df2,
                        alpha=alpha,
                        lambda.u.seq = seq(0.02,0.1,0.02),
                        lambda.z.seq = seq(0.5,5,0.5))
    res<-err2$res
    uz <- cbind.data.frame(res$u, res$z)
    uz$low <- res$beta0 + res$u + res$z*res$beta1 + res$beta3*age
    uz$up  <- res$beta0 + res$u + res$z*res$beta2 + res$beta3*age
    uz$id <- unique(df2$subject)
    colnames(uz) <- c("u", "z", "low", "up","id")
    }

  return(list(res=res,
         uz=uz, df2=df2,
         out=err$out))
}

ibs<-ibs_all_h
ibs<-ibs[,c(1:2,8:9,3:7,10:28)]
res.est<-matrix(0,nrow=ncol(ibs_all_h)-4,ncol=14);
uz_h<-0;df_h<-0;cnt<-1#ncol(ibs)
for (i in 5:ncol(ibs)) {
  db<-ibs[,c(1:2,i)]
  df2<-ibs[,c(1:4,i)]
  iri<-iri_fun(db=db,df2=df2,alpha=0.05)
  res<-iri$res
  uz<-iri$uz
  uz$group<-"Healthy"
  uz$label<-colnames(df2)[5]
  df<-iri$df2
  df$group<-"Healthy"
  df$label<-colnames(df2)[5]
  save(uz, file =paste0("./processed_data/uz_Healthy_",colnames(df2)[5],".RData"))
  save(df, file =paste0("./processed_data/df_Healthy_",colnames(df2)[5],".RData"))
  
  if(iri$out==FALSE){
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,res$beta4,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
    
  }else{
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,NA,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
  }
  
  cnt<-cnt+1
  uz_h<-rbind(uz_h,uz)
  df_h<-rbind(df_h,df)
}
est_h<-res.est
save(est_h, file =paste0("./processed_data/est_Healthy.RData"))

ibs<-ibs_all_c
ibs<-ibs[,c(1:2,8:9,3:7,10:28)]
res.est<-matrix(0,nrow=ncol(ibs_all_c)-4,ncol=14);
uz_c<-0;df_c<-0;cnt<-1#ncol(ibs)
for (i in 5:ncol(ibs)) {
  db<-ibs[,c(1:2,i)]
  df2<-ibs[,c(1:4,i)]
  iri<-iri_fun(db=db,df2=df2,alpha=0.05)
  res<-iri$res
  uz<-iri$uz
  uz$group<-"IBS-C"
  uz$label<-colnames(df2)[5]
  df<-iri$df2
  df$group<-"IBS-C"
  df$label<-colnames(df2)[5]
  save(uz, file =paste0("./processed_data/uz_IBS-C_",colnames(df2)[5],".RData"))
  save(df, file =paste0("./processed_data/df_IBS-C_",colnames(df2)[5],".RData"))
  
  if(iri$out==FALSE){
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,res$beta4,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
    
  }else{
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,NA,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
    }
  
  cnt<-cnt+1
  uz_c<-rbind(uz_c,uz)
  df_c<-rbind(df_c,df)
}
est_c<-res.est
save(est_c, file =paste0("./processed_data/est_IBS-C"))

ibs<-ibs_all_d
ibs<-ibs[,c(1:2,8:9,3:7,10:28)]
res.est<-matrix(0,nrow=ncol(ibs_all_d)-4,ncol=14);
uz_d<-0;df_d<-0;cnt<-1#ncol(ibs)
for (i in 5:ncol(ibs)) {
  db<-ibs[,c(1:2,i)]
  df2<-ibs[,c(1:4,i)]
  iri<-iri_fun(db=db,df2=df2,alpha=0.05)
  res<-iri$res
  uz<-iri$uz
  uz$group<-"IBS-D"
  uz$label<-colnames(df2)[5]
  df<-iri$df2
  df$group<-"IBS-D"
  df$label<-colnames(df2)[5]
  save(uz, file =paste0("./processed_data/uz_IBS-D_",colnames(df2)[5],".RData"))
  save(df, file =paste0("./processed_data/df_IBS-D_",colnames(df2)[5],".RData"))
  
  if(iri$out==FALSE){
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,res$beta4,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
    
  }else{
    res.est[cnt,]<-c(df$group[1],colnames(df2)[5],length(unique(df$subject)),iri$out,
                     res$beta0,res$beta1,res$beta2,res$beta3,NA,
                     res$lambda.u,res$lambda.z,res$cov.subj,res$cov.time,res$cov.tot)
    colnames(res.est)<-c("group","label","nr.subject","only.age",
                         "beta0","beta1","beta2","beta3","beta4",
                         "l.u","l.z","SEC","TEC","EC")
  }
  
  cnt<-cnt+1
  uz_d<-rbind(uz_d,uz)
  df_d<-rbind(df_d,df)
}
est_d<-res.est
save(est_d, file =paste0("./processed_data/est_IBS-D.RData"))


uz_all<-rbind(uz_h[-1,],uz_c[-1,],uz_d[-1,])
df_all<-rbind(df_h[-1,],df_c[-1,],df_d[-1,])
est_all<-rbind(est_h,est_c,est_d)

save(uz_all, file =paste0("./processed_data/uz_all.RData"))
save(df_all, file =paste0("./processed_data/df_all.RData"))
save(est_all, file =paste0("./processed_data/est_all.RData"))
