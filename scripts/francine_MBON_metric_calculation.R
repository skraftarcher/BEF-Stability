# calculate metrics for francine effects on stability for MBON data

# load packages and download up to date data
source("scripts/install_packages_function.R")
lp("tidyverse")
lp("patchwork")
lp("vegan")

# load data 
mba<-read.csv("odata/MBON wide abundance community dataset.csv")
mbe<-read.csv("odata/MBON full environmental dataset for community.csv")

nts<-mbe%>%
  group_by(site,season,yr)%>%
  summarize(nt=n())

mbep<-mbe%>%
  bind_cols(mba)%>%#add on the community data before changing dataset structure
  group_by(site)%>%
  mutate(sampling=paste(season,yr))%>%
  # filter out the data we want from prior to Francine
  filter(sampling=="Summer 2024" & site %in% c("BB1","LUMO3","LUMO6", "TB1"))

mbea<-mbe%>%
  bind_cols(mba)%>%#add on the community data before changing dataset structure
  group_by(site)%>%
  mutate(sampling=paste(season,yr))%>%
  # filter out the data we want from after Francine
  filter(sampling=="Fall 2024" & site %in% c("BB1","LUMO3","LUMO6", "TB1"))

# put all non-community variables "up front"
mbep<-mbep[,c(1:28,149,29:148)]
mbea<-mbea[,c(1:28,149,29:148)]

# spatial variability----
source("scripts/metrics/cv.R")
source("scripts/metrics/compositionalturnover.R")

prior.space.within.cv<-cvhome(ds=mbep,
                       nc=29,
                       tors="spatial",
                       g.vars = c("site","sampling","tray"),
                       s.vars = c("site","sampling"),
                       dset="m")

after.space.within.cv<-cvhome(ds=mbea,
                         nc=29,
                         tors="spatial",
                         g.vars = c("site","sampling","tray"),
                         s.vars = c("site","sampling"),
                         dset="m")

prior.space.between.cv<-cvhome(ds=mbep,
                                nc=29,
                                tors="spatial",
                                g.vars = c("basin","site","sampling"),
                                s.vars = c("basin","sampling"),
                                dset="m")

after.space.between.cv<-cvhome(ds=mbea,
                                  nc=29,
                                  tors="spatial",
                                  g.vars = c("basin","site","sampling"),
                                  s.vars = c("basin","sampling"),
                                  dset="m")

# resistance----

mbea.env<-mbea[,1:29]
mbea.com<-mbea[,-1:-29]

mbea.com.hel<-decostand(mbea.com,"hellinger")

mbea.rda<-rda(mbea.com.hel)

mbea.scores<-data.frame(scores(mbea.rda,choices=c(1,2),display="sites"))

mbea.env<-bind_cols(mbea.env,mbea.scores)


mbep.env<-mbep[,1:29]
mbep.com<-mbep[,-1:-29]

mbep.com.hel<-decostand(mbep.com,"hellinger")

mbep.rda<-rda(mbep.com.hel)

mbep.scores<-data.frame(scores(mbep.rda,choices=c(1,2),display="sites"))

mbep.env<-bind_cols(mbep.env,mbep.scores)%>%
  group_by(site)%>%
  mutate(cent.x=mean(PC1),
         cent.y=mean(PC2),
         dists=1/sqrt((PC1-cent.x)^2+(PC2-cent.y)^2))

centroids<-mbep.env%>%
  select(site,cent.x,cent.y)%>%
  distinct()

mbea.env<-left_join(mbea.env,centroids)%>%
  mutate(dists=1/sqrt((PC1-cent.x)^2+(PC2-cent.y)^2))
