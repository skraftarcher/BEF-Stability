# simple visualizations of MBON data

# load packages and download up to date data
source("scripts/install_packages_function.R")
lp("tidyverse")
lp("patchwork")
theme_set(theme_bw()+theme(panel.grid=element_blank(),
                           axis.title = element_text(size=16),
                           axis.text = element_text(size=12),
                           axis.text.x = element_text(angle=45,hjust=1),
                           legend.title=element_text(size=16)))

# load data
mba<-read.csv("odata/MBON wide abundance community dataset.csv")
mbe<-read.csv("odata/MBON full environmental dataset for community.csv")


# plot species richness over time
mbe.sum<-mbe%>%
  group_by(date.retrieved,site)%>%
  summarize(oys.m=mean(oys.count/(.25*.25)),
            tax.rich.m=mean(taxa.rich),
            tax.rich.sd=sd(taxa.rich))
ggplot(mbe.sum%>%
         filter(!site%in%c("SL1","LUMO2","LUMO4","LUMO5")))+
  geom_line(aes(x=date.retrieved,y=tax.rich.m,group=site,color=site),linewidth=2)+
  ylab("Number of taxa")+
  xlab("Date")+
  scale_color_viridis_d(name="Site",option="B",begin=.2,end=.8)

ggsave("figures/change in taxa richness over time.jpg",width=10,height=6)

# bring in larger drought dataset
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
  arrange(ValidStart)%>%
  filter(ValidStart>=min(ymd(mbe$date.retrieved)))%>%
  filter(ValidEnd<=max(ymd(mbe$date.retrieved)))

ggplot(data=drought.t)+
  geom_line(aes(x=ValidStart,y=sev))

# oysters and biodiversity
ggplot(data=mbe%>%
         filter(site%in%c("BB2","BB3","LUMO3","LUMO6")),
       aes(x=oys.count/(.25*.25),y=taxa.rich,color=season))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~site,scales="free")+
  ylab("Number of taxa")+
  xlab("Live oysters per m2")

