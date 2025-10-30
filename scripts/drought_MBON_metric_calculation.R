# calculate metrics for drought effects on stability for MBON data

# load packages and download up to date data
source("scripts/install_packages_function.R")
lp("tidyverse")
lp("patchwork")

# load data
mba<-read.csv("odata/MBON wide abundance community dataset.csv")
mbe<-read.csv("odata/MBON environmental data WITH DROUGHT.csv")

nts<-mbe%>%
  group_by(site,season,yr)%>%
  summarize(nt=n())
# identify drought periods
mbe.dr<-mbe%>%
  filter(season=="Winter")%>%
  filter(yr==2023)

mbe2<-mbe%>%
  bind_cols(mba)%>%#add on the community data before changing dataset structure
  anti_join(mbe.dr)%>% # pull out middle summer data
  filter(!site %in% c("LUMO2","LUMO5","LUMO4","SL1"))%>%
  mutate(d6.max=ifelse(max6week.sev>=3,1,0),#assign if in drought or not
         d2.max=ifelse(max2week.sev>=3,1,0))%>%
  filter(ymd(date.retrieved)<="2024-09-11")# filter out data after Francine to avoid hurricane conflicts


# create distinct datasets
# drought and not based on 6 week max severity
d6<-mbe2%>%
  filter(d6.max==1)
nd6<-mbe2%>%
  filter(d6.max==0)


# now assign sampling # based on group (drought v no drought, site)

d6<-d6%>%
  group_by(site)%>%
  mutate(sampling=paste(season,yr))

unique(d6$sampling)

d6$sampling<-factor(d6$sampling,levels=c("Spring 2022",
                                         "Summer 2022",
                                         "Summer 2023",
                                         "Fall 2023",
                                         "Winter 2024"))
nd6<-nd6%>%
  group_by(site)%>%
  mutate(sampling=paste(season,yr))

unique(nd6$sampling)

nd6$sampling<-factor(nd6$sampling,levels=c("Fall 2022",
                                         "Summer 2023",
                                         "Spring 2024",
                                         "Summer 2024"))

# put all non-community variables "up front"
d6<-d6[,c(1:34,155:157,35:154)]
nd6<-nd6[,c(1:34,155:157,35:154)]


# temporal variability
source("scripts/metrics/cv.R")
source("scripts/metrics/compositionalturnover.R")

drought.tempcv<-cvhome(ds=d6,
                       nc=37,
                       tors="temporal",
                       g.vars = c("site","sampling"),
                       t.vars = c("site"),
                       dset="m")

nodrought.tempcv<-cvhome(ds=nd6,
                       nc=37,
                       tors="temporal",
                       g.vars = c("site","sampling"),
                       t.vars = c("site"),
                       dset="m")


drought.space.within.cv<-cvhome(ds=d6,
                       nc=37,
                       tors="spatial",
                       g.vars = c("site","sampling","tray"),
                       s.vars = c("site","sampling"),
                       dset="m")

nodrought.space.within.cv<-cvhome(ds=nd6,
                         nc=37,
                         tors="spatial",
                         g.vars = c("site","sampling","tray"),
                         s.vars = c("site","sampling"),
                         dset="m")

drought.space.between.cv<-cvhome(ds=d6,
                                nc=37,
                                tors="spatial",
                                g.vars = c("basin","site","sampling"),
                                s.vars = c("basin","sampling"),
                                dset="m")

nodrought.space.between.cv<-cvhome(ds=nd6,
                                  nc=37,
                                  tors="spatial",
                                  g.vars = c("basin","site","sampling"),
                                  s.vars = c("basin","sampling"),
                                  dset="m")

# compositional turnover
d6.site<-d6%>%
  group_by(site,sampling)%>%
  mutate(nt=n())%>%
  ungroup()%>%
  pivot_longer(38:157,names_to = "taxa",values_to="abund")%>%
  group_by(site,sampling,taxa)%>%
  mutate(abund=sum(abund)/nt)%>%
  select(site,sampling,taxa,abund)%>%
  distinct()%>%
  pivot_wider(names_from=taxa,values_from=abund)%>%
  filter(site %in% c("BB1","LUMO6","TB1"))

drought.compt<-cturn.jadist(ds.env=d6.site[,1:2],ds.com = d6.site[,-1:-2],site.id = "site")


nd6.site<-nd6%>%
  group_by(site,sampling)%>%
  mutate(nt=n())%>%
  ungroup()%>%
  pivot_longer(38:157,names_to = "taxa",values_to="abund")%>%
  group_by(site,sampling,taxa)%>%
  mutate(abund=sum(abund)/nt)%>%
  select(site,sampling,taxa,abund)%>%
  distinct()%>%
  pivot_wider(names_from=taxa,values_from=abund)%>%
  filter(site %in% c("BB1","LUMO6","TB1"))

nodrought.compt<-cturn.jadist(ds.env=nd6.site[,1:2],ds.com = nd6.site[,-1:-2],site.id = "site")


# 