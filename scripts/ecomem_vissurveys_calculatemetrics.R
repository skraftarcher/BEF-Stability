# This script will  calculate stability metrics on
# visual survey data

# load packages and download up to date data
source("scripts/install_packages_function.R")
source("scripts/metrics/cv.R")
source("scripts/metrics/compositionalturnover.R")
lp("readxl")
lp("codyn")

# bring in data----
viswide<-read.csv("odata/EcoMEM_wide_vissurvey_datawithdummy.csv")%>%
  mutate(Scar=ifelse(treat %in% c("SU","SG"),"scar","no.scar"),
         Graze=ifelse(treat %in% c("SG","UG"),"graze","no.graze"),
         Bay=ifelse(blockID %in% c("SA1","SA2","SA3","SA4","SA5","SA6"),"StAndrew","StJoe"))
viswide<-viswide[,c(1:4,72:74,5:71)]

vislong<-viswide%>%
  pivot_longer(-1:-7,names_to="taxa",values_to="abund")%>%
  separate(sampling,into=c("s","time.var"),sep=-1,remove=FALSE,convert = T)%>%
  select(-s)%>%
  filter(taxa!="dummy")

# time covariance----
tvar<-cvhome(ds=viswide,
         nc=7)
# space covariance----
svar<-cvhome(ds=viswide,
             nc=7,
             tors="spatial")

# compositional turnover----
cturn<-cturn.jadist(ds.env=viswide[,1:7],ds.com=viswide[,-1:-7],site.id="plotID")

# extinctions----

dis<-turnover(
  df=vislong,
  time.var="time.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="plotID",
  metric="disappearance")

#invasions----
app<-turnover(
  df=vislong,
  time.var="time.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="plotID",
  metric="appearance")

# calculate resistance----
# think through this one more
viswide2<-viswide%>%
  filter(sampling=="s5")

vis.hel<-decostand(viswide2[,-1:-7],method="hellinger")
vis.rda<-rda(vis.hel)
plot(vis.rda)
vis.ris<-data.frame(viswide2[,1:9],
                    scores(vis.rda,choices=c(1,2),display="sites"))

vis.cent<-vis.ris%>%
  filter(Graze=="No.graze" & Scar=="No.scar")%>%
  group_by(Bay)%>%
  summarize(cx=mean(PC1),
            cy=mean(PC2))

vis.cent2<-vis.ris%>%
  filter(Graze=="No.graze" & Scar=="No.scar")


vis.ris2<-vis.ris%>%
  left_join(vis.cent)%>%
  mutate(resistance=1/sqrt((PC1-cx)^2+(PC2-cy)^2))%>%
  select(blockID,plotID,sampling,Bay,Scar,Graze,resistance)


