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

vis2<-data.frame(viswide[,1:7],spr=specnumber(viswide[,c(-1:-7)]))
vis3<-vislong%>%
  group_by(blockID,plotID,sampling,Scar,Graze)%>%
  summarize(abund=sum(abund))%>%
  ungroup()%>%
  left_join(vis2)

eco.sg<-read.csv("odata/EcoMEM_thalassia p cover shoot density and canopy by quadrat long.csv")


#######################
## OVERALL QUESTION ##
######################

## Response variables we will use
# MBON 
#   Oyster volume
#   Total Abundance/density
#   animal species richness
# ECOMEM
#   seagrass shoot density
#   visual survey species richness
#   visual survey total abundance


# time covariance----
############################################
## QUESTION - should we detrend the data? ##
############################################
# look at temporal variation in species richness and species abundance over time
ggplot(data=vis3%>%filter(sampling!="s1"),aes(x=sampling,y=spr,group=plotID,color=Bay))+geom_line()+facet_grid(Scar~Graze)
ggplot(data=vis3,aes(x=sampling,y=abund,group=plotID,color=Bay))+geom_line()+facet_grid(Scar~Graze)
# looks like there might be a trend in species richness across treatments

tvar.ecomem<-cvhome(ds=viswide,
         nc=7)
tvar.sg<-cvhome(ds=eco.sg[,c(-8,-10,-12)],
                nc=8,
                tors = "temporal",
                g.vars=c("blockID","plotID","sampling","mnths","bay","scar","graze"),
                t.vars=c("blockID","plotID","scar","graze"))%>%
  select(-temp.cvs,-temp.dt.cvs)%>%
  rename(temp.cvsg=temp.cva,
         temp.dt.cvsg=temp.dt.cva)

# spatial covariance----
############################################
## QUESTION - should we detrend the data? ##
############################################
svar<-cvhome(ds=viswide,
             nc=7,
             tors="spatial")

# compositional turnover----

cturn<-cturn.jadist(ds.env=viswide[,1:7],ds.com=viswide[,-1:-7],site.id="plotID")

# extinctions----
#################################################
## from codyn package, pretty straight forward ##
#################################################
dis<-turnover(
  df=vislong,
  time.var="time.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="plotID",
  metric="disappearance")

#invasions----
#################################################
## from codyn package, pretty straight forward ##
#################################################
app<-turnover(
  df=vislong,
  time.var="time.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="plotID",
  metric="appearance")

# calculate resistance----
# think through this one more

## From Donohue et al. 2013 
## Measured as the inverse of the Euclidian distance from each experimental plot to the 
## centroid of the unmanipulated uncaged treatment at the end of the experiment 
## [based  on Bray-Curtis similarity matrices calculated from log (x + 1)-transformed algal cover  data]. 
## This provides a holistic measure of the extent of change in algal community  structure over the duration 
## of the experiment, over and above natural background  dynamics. The resistance of the caged plots with no 
## experimental consumer  removals did not differ statistically from the inverse of the distance of the uncaged 
## plots to their treatment centroid (t3 = 0.24, P = 0.82).

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


