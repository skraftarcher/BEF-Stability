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


# temporal variability----
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

# compositional turnover ----
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


# extinctions----
lp("codyn")

d6.long<-d6%>%
  filter(site %in% c("BB1","LUMO6","TB1"))%>%
  pivot_longer(38:157,names_to = "taxa",values_to="abund")%>%
  group_by(site,sampling,taxa)%>%
  mutate(abund=sum(abund))%>%
  select(site,sampling,taxa,abund)%>%
  distinct()%>%
  mutate(samp.var=as.numeric(sampling))
  


d6dis<-turnover(
  df=d6.long,
  time.var="samp.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="site",
  metric="disappearance")


nd6.long<-nd6%>%
  filter(site %in% c("BB1","LUMO6","TB1"))%>%
  mutate(keep=case_when(
    site %in% c("BB1","TB1")~1,
    site=="LUMO6" & sampling %in% c("Fall 2022","Summer 2024")~1,
    TRUE~0))%>%
  filter(keep==1)%>%
  pivot_longer(38:157,names_to = "taxa",values_to="abund")%>%
  group_by(site,sampling,taxa)%>%
  mutate(abund=sum(abund))%>%
  select(site,sampling,taxa,abund)%>%
  distinct()%>%
  mutate(samp.var=as.numeric(sampling))



nd6dis<-turnover(
  df=nd6.long,
  time.var="samp.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="site",
  metric="disappearance")


#invasions ----
d6inv<-turnover(
  df=d6.long,
  time.var="samp.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="site",
  metric="appearance")

nd6inv<-turnover(
  df=nd6.long,
  time.var="samp.var",
  species.var="taxa",
  abundance.var="abund",
  replicate.var="site",
  metric="appearance")

# resistance----

d6.env<-d6[,1:37]
d6.com<-d6[,-1:-37]

lp("vegan")

d6.com.hel<-decostand(d6.com,"hellinger")

d6.rda<-rda(d6.com.hel)

d6.scores<-data.frame(scores(d6.rda,choices=c(1,2),display="sites"))

d6.env<-bind_cols(d6.env,d6.scores)


nd6.env<-nd6[,1:37]
nd6.com<-nd6[,-1:-37]

nd6.com.hel<-decostand(nd6.com,"hellinger")

nd6.rda<-rda(nd6.com.hel)

nd6.scores<-data.frame(scores(nd6.rda,choices=c(1,2),display="sites"))

nd6.env<-bind_cols(nd6.env,nd6.scores)%>%
  group_by(site)%>%
  mutate(cent.x=mean(PC1),
         cent.y=mean(PC2),
         dists=1/sqrt((PC1-cent.x)^2+(PC2-cent.y)^2))

centroids<-nd6.env%>%
  select(site,cent.x,cent.y)%>%
  distinct()

d6.env<-left_join(d6.env,centroids)%>%
  mutate(dists=1/sqrt((PC1-cent.x)^2+(PC2-cent.y)^2))


# organize, visualize, and analyze data ----
theme_set(theme_bw()+theme(panel.grid = element_blank()))

# look at temporal variability
temp.cv<-bind_rows(bind_cols(drought.tempcv,data.frame(dgroup="drought")),
                   bind_cols(nodrought.tempcv,data.frame(dgroup="nodrought")))

(temp.cva.t<-t.test(temp.cv$temp.cva~temp.cv$dgroup))
(temp.cvspr.t<-t.test(temp.cv$temp.cvs~temp.cv$dgroup))

(temp.cva.kw<-kruskal.test(temp.cv$temp.cva~temp.cv$dgroup))
(temp.cvspr.kw<-kruskal.test(temp.cv$temp.cvs~temp.cv$dgroup))

ggplot(data=temp.cv)+
  geom_boxplot(aes(x=dgroup,y=temp.cva))

ggplot(data=temp.cv)+
  geom_boxplot(aes(x=dgroup,y=temp.cvs))

# look at within site spatial variability
site.cv<-temp.cv<-bind_rows(bind_cols(drought.space.within.cv,data.frame(dgroup="drought")),
                            bind_cols(nodrought.space.within.cv,data.frame(dgroup="nodrought")))

(site.cva.t<-t.test(site.cv$space.cva~site.cv$dgroup))
(site.cvspr.t<-t.test(site.cv$space.cvs~site.cv$dgroup))

(site.cva.kw<-kruskal.test(site.cv$space.cva~site.cv$dgroup))
(site.cvspr.kw<-kruskal.test(site.cv$space.cvs~site.cv$dgroup))

ggplot(data=site.cv)+
  geom_boxplot(aes(x=dgroup,y=space.cva))

ggplot(data=site.cv)+
  geom_boxplot(aes(x=dgroup,y=space.cvs))

# look at between site spatial variability
site2.cv<-bind_rows(bind_cols(drought.space.between.cv,data.frame(dgroup="drought")),
                            bind_cols(nodrought.space.between.cv,data.frame(dgroup="nodrought")))

(site2.cva.t<-t.test(site2.cv$space.cva~site2.cv$dgroup))
(site2.cvspr.t<-t.test(site2.cv$space.cvs~site2.cv$dgroup))

(site2.cva.kw<-kruskal.test(site2.cv$space.cva~site2.cv$dgroup))
(site2.cvspr.kw<-kruskal.test(site2.cv$space.cvs~site2.cv$dgroup))

ggplot(data=site2.cv)+
  geom_boxplot(aes(x=dgroup,y=space.cva))

ggplot(data=site2.cv)+
  geom_boxplot(aes(x=dgroup,y=space.cvs))

# comp.turnover----
turnover<-bind_rows(bind_cols(drought.compt,data.frame(dgroup="drought")),
                    bind_cols(nodrought.compt,data.frame(dgroup="nodrought")))

(turn.t<-t.test(turnover$dsts~turnover$dgroup))

ggplot(data=turnover)+
  geom_boxplot(aes(x=dgroup,y=dsts))
# no difference

# extinctions----
ext<-bind_rows(bind_cols(d6dis,data.frame(dgroup="drought")),
                 bind_cols(nd6dis,data.frame(dgroup="nodrought")))

(ext.t<-t.test(ext$disappearance~ext$dgroup))


ggplot(data=ext)+
  geom_boxplot(aes(x=dgroup,y=disappearance))

# invasions----
inv<-bind_rows(bind_cols(d6inv,data.frame(dgroup="drought")),
               bind_cols(nd6inv,data.frame(dgroup="nodrought")))

(inv.t<-t.test(inv$appearance~inv$dgroup))

ggplot(data=inv)+
  geom_boxplot(aes(x=dgroup,y=appearance))

# resistance
d6.res<-d6.env%>%
  select(site,sampling,basin,dists,PC1,PC2)%>%
  mutate(dgroup="drought")


nd6.res<-nd6.env%>%
  select(site,sampling,basin,dists,PC1,PC2)%>%
  mutate(dgroup="nodrought")

res<-bind_rows(d6.res,nd6.res)

(res.t<-t.test(res$dists~res$dgroup))
(res.kw<-kruskal.test(res$dists~res$dgroup))

ggplot(data=res)+
  geom_boxplot(aes(x=dgroup,y=dists))

# make plot of resistance
ggplot()+
  stat_ellipse(aes(x=PC1,y=PC2,color=site),data=res%>%filter(dgroup=="nodrought"))+
  # geom_point(data=centroids,aes(x=cent.x,y=cent.y,color=site),shape=15,size=6)+
  # geom_point(aes(x=PC1,y=PC2,color=site),size=4,alpha=.5,data=res%>%filter(dgroup=="drought"))+
  stat_ellipse(aes(x=PC1,y=PC2,color=site,fill=site),alpha=.5,geom="polygon",data=res%>%filter(dgroup=="drought"))+
  facet_wrap(~site)
  

ggplot()+
  stat_ellipse(aes(x=PC1,y=PC2,color=site,fill=site),geom="polygon",alpha=.5,data=res%>%filter(dgroup=="nodrought"))+
  # geom_point(data=centroids,aes(x=cent.x,y=cent.y,color=site),shape=15,size=6)+
  # geom_point(aes(x=PC1,y=PC2,color=site),size=4,alpha=.5,data=res%>%filter(dgroup=="drought"))+
  stat_ellipse(aes(x=PC1,y=PC2,color=site,fill=site),alpha=.5,geom="polygon",data=res%>%filter(dgroup=="drought"))+
  facet_grid(basin~dgroup)

