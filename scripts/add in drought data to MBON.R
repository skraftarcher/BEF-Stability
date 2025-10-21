# script to organize drought data
# drought data is downloaded from here https://droughtmonitor.unl.edu/DmData/DataDownload/ComprehensiveStatistics.aspx


source("scripts/install_packages_function.R")
lp("tidyverse")

envdat<-read.csv("odata/MBON full environmental dataset for community.csv")
# add in dates of drought - sample will either be before, during, or after
drought.t<-read.csv("odata/drought_October25.csv")%>%
  filter(County %in% c("Terrebonne Parish"))%>%
  mutate(is.drought=ifelse(None==100,0,1),
         sev=case_when(
           D4!=0~5,
           D3!=0 & D4==0 ~4,
           D2!=0 &D3==0 &D4==0~3,
           D1!=0 & D2==0 &D3==0 & D4==0~2,
           D0!=0 & D1==0 & D2==0 &D3==0 & D4==0~1,
           TRUE~0),
         ValidStart=ymd(ValidStart),
         ValidEnd=ymd(ValidEnd))%>%
  arrange(ValidStart)

drought.b<-read.csv("odata/drought_October25.csv")%>%
  filter(County %in% c("Lafourche Parish"))%>%
  mutate(is.drought=ifelse(None==100,0,1),
         sev=case_when(
           D4!=0~5,
           D3!=0 & D4==0 ~4,
           D2!=0 &D3==0 &D4==0~3,
           D1!=0 & D2==0 &D3==0 & D4==0~2,
           D0!=0 & D1==0 & D2==0 &D3==0 & D4==0~1,
           TRUE~0),
         ValidStart=ymd(ValidStart),
         ValidEnd=ymd(ValidEnd))%>%
  arrange(ValidStart)

envdat$max6week.sev<-NA
envdat$mean6week.sev<-NA
envdat$max2week.sev<-NA
envdat$mean2week.sev<-NA

envdat$date.retrieved<-ymd(envdat$date.retrieved)

for(i in 1:nrow(envdat)){
  t1<-envdat$date.retrieved[i]
  t6w<-t1-weeks(6)
  t2w<-t1-weeks(2)
  if(envdat$site[i] %in% c("BB1","BB2","BB3"))td<-drought.b
  if(!envdat$site[i] %in% c("BB1","BB2","BB3"))td<-drought.t 
  t4<-filter(td,ValidStart<=t1)
  t42w<-filter(t4,ValidStart>=t2w)
  t46w<-filter(t4,ValidStart>=t6w)
  envdat$max6week.sev[i]<-max(t46w$sev)
  envdat$mean6week.sev[i]<-mean(t46w$sev)
  envdat$max2week.sev[i]<-max(t42w$sev)
  envdat$mean2week.sev[i]<-mean(t42w$sev)
}


envdat$hur<-"before.francine"
envdat$hur[envdat$date.retrieved>="2024-09-11"]<-"after.francine"

write.csv(envdat,"wdata/MBON environmental data WITH DROUGHT.csv")

