# function to calculate coefficient of variation (temporal or spatial)
source("scripts/install_packages_function.R")
lp("tidyverse")
lp("vegan")

# this function's input is two data sets the environment dataset and 
# the community dataset
# site.id is the column name of the unique site identifier

# compositional turnover
cturn.jadist<-function(ds.env,ds.com,site.id){
  ds.com$dummy<-1
  t.jac<-as.matrix(vegdist(decostand(ds.com,method="pa"),method = "jaccard"))
  env<-ds.env
  env$nr<-rownames(ds.com)

# write a function to extract dissimilarities
compturn<-function(site.id,site,distmat,env){
  t1<-env[,grep(x=colnames(env),pattern=site.id)]
  t2<-env$nr[t1==site]
  td<-distmat[t2,t2]
  t3<-1:nrow(td)
  t4<-NA
  for(i in 1:(length(t3)-1))t4[i]<-paste0("s",t3[i],".s",t3[i+1])
  return(data.frame(site.id=site,
                    dsts=td[row(td)==(col(td)-1)],
                    tms=t4))
}

sts<-data.frame(unique(env[,grep(x=colnames(env),pattern=site.id)]))
cturn<-compturn(site=as.vector(sts[1,1]),
           site.id=site.id,
           distmat = t.jac,
           env=env)
for(i in 2:nrow(sts))cturn<-bind_rows(cturn,
                                        compturn(site=as.vector(sts[i,1]),
                                           site.id=site.id,
                                           distmat = t.jac,
                                           env=env))
return(cturn)
}


