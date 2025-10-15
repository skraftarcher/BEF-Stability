# This script will  calculate stability metrics on
# visual survey data

# load packages and download up to date data
source("scripts/install_packages_function.R")
# source("scripts/download_all_experiment_data-EX.R")
# 2
source("scripts/metrics/cv.R")
source("scripts/metrics/compositionalturnover.R")
lp("readxl")
lp("codyn")

# bring in data----
samplings<-read_xlsx(paste0("odata/downloaded_2025-07-17_EXPERIMENT - Sampling Dates.xlsx"),sheet=1)

vis<-read_xlsx(paste0("odata/downloaded_2025-07-17_EXPERIMENT - Visual Survey Data.xlsx"),sheet=2)%>%
  filter(!is.na(plotID))

# size<-read_xlsx(paste0("odata/downloaded_2025-07-17_EXPERIMENT - Visual Survey Data.xlsx"),sheet=3)

tax<-read_xlsx(paste0("odata/downloaded_2025-07-17_EXPERIMENT - LAB - Dipnet.xlsx"),sheet=4)

# organize data----
viswide<-vis[-grep("egg",vis$taxa),]%>%
  select(-scientific.name,-QAQC,-common.name,-Notes)%>%
  filter(taxa!="animal-0")%>%
  rename(taxaID=taxa)%>%
  left_join(tax)%>%
  group_by(date,blockID,plotID,UpdateTaxaID,Surveyor)%>%
  summarize(abund=sum(total.n))%>%
  mutate(sampling=case_when(
    date >= samplings$Start.date[1] & date <=samplings$End.date[1]~samplings$sampling[1],
    date >= samplings$Start.date[2] & date <=samplings$End.date[2]~samplings$sampling[2],
    date >= samplings$Start.date[3] & date <=samplings$End.date[3]~samplings$sampling[3],
    date >= samplings$Start.date[4] & date <=samplings$End.date[4]~samplings$sampling[4],
    date >= samplings$Start.date[5] & date <=samplings$End.date[5]~samplings$sampling[5]),
    time.var=case_when(
      sampling=="s1"~1,
      sampling=="s2"~2,
      sampling=="s3"~3,
      sampling=="s4"~4,
      sampling=="s5"~5,
    ))%>%
  pivot_wider(names_from = UpdateTaxaID,values_from = abund,values_fill = 0)

plot.info<-viswide%>%
  select(blockID,plotID,date)%>%
  distinct()
plot.info$Bay<-"SA"
plot.info$Bay[grep("SJ",plot.info$plotID)]<-"SJ"
plot.info$Scar<-"Scar"
plot.info$Scar[grep("UU",plot.info$plotID)]<-"No.scar"
plot.info$Scar[grep("UG",plot.info$plotID)]<-"No.scar"
plot.info$Graze<-"Graze"
plot.info$Graze[grep("SU",plot.info$plotID)]<-"No.graze"
plot.info$Graze[grep("UU",plot.info$plotID)]<-"No.graze"

viswide<-left_join(plot.info,viswide)

vislong<-viswide%>%
  pivot_longer(-1:-9,names_to="taxa",values_to="abund")

# time covariance----
tvar<-cvhome(ds=viswide,
         nc=9)
# space covariance----
svar<-cvhome(ds=viswide,
             tors="spatial")
# compositional turnover----
cturn<-cturn.jadist(ds.env=viswide[,1:9],ds.com=viswide[,-1:-9],site.id="plotID")

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
  mutate(dummy=1)%>%
  filter(sampling=="s5")

vis.hel<-decostand(viswide2[,-1:-9],method="hellinger")
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


