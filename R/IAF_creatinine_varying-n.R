library(tidyverse)
# load data of several iaf metabolomics features 
load("./data/iaf_metabolomics-cvd_difftime.RData")
pars<-unique(all$var)

start<-Sys.time()
datres<-0; iter<-c(4:11) # only for data with all 11 observation (10 as training)
for(i in pars){
  for (j in iter) {
    # prepare data of each feature in different n, separate the last observation as training
    db<-filter(all, var==i & time %in% 1:j)
    subject<-unique(all$subject)
    
    dbu<-0
    for(s in subject){
      dt<-db[db$subject==s,]
      dt$obs2[nrow(dt)]<-"last"
      dbu<-rbind(dbu,dt)
    }
    dbu<-dbu[-1,]
    df2<-dbu[dbu$obs2=="prev",]
    
    # function of to catch errors while computing IRI
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
    
    source("./R/main/JQM_Function.R") # change the source function following the model specification (what covariates to include?)
    alpha=0.05
    subject<-unique(df2$subject)
    len<-length(unique(df2$subject))
    err<-catch.err.lqmm(db=df2,
                        alpha=alpha,
                        lambda.u.seq = seq(0.02,0.1,0.02),
                        lambda.z.seq = seq(0.5,5,0.5))
    # organise and save the results
    res<-err$res
    uz <- cbind.data.frame(res$u, res$z)
    uz$id <- unique(df2$subject)
    colnames(uz) <- c("u", "z","id")

    df1<-dbu[dbu$obs2=="last",]
    uz<-uz %>% left_join(df1[,c(1,3:5)], by=c("id"="subject"))
    uz$low <- res$beta0 + res$u + res$z*res$beta1
    uz$up  <- res$beta0 + res$u + res$z*res$beta2
    uz$len<-uz$up-uz$low
    uz$cov<-ifelse(uz$y>=uz$low & uz$y<=uz$up,1,0)
    uz$lenprop<-NA
    for(k in 1:nrow(uz)){
      avglen<-mean(uz[!rownames(uz) %in% k,]$len)
      if(uz$len[k]>avglen){
        uz$lenprop[k]<-2-(uz$len[k]/avglen)
      }else{uz$lenprop[k]<-uz$len[k]/avglen}
    }
    uz$var<-i; uz$n<-j-1
    save(uz, file=paste0("processed/results0612/iri_sex_",i,"_n",j-1,".RData"))
    
    var<-i; n<-j-1
    TEC<-mean(uz$cov)
    avglen<-mean(uz$len); lenprop<-mean(uz$lenprop)
    uzsum<-data.frame(var,n,TEC,avglen,lenprop)
    datres<-rbind(datres, uzsum)
  }
}
end<-Sys.time()
end-start

datres<-datres[-1,]
save(datres, file="processed/results/datres_n3-10.RData")

