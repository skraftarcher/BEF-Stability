# function to calculate coefficient of variation (temporal or spatial)
source("scripts/install_packages_function.R")
lp("tidyverse")
lp("vegan")
lp("EnvStats")
# need a dataset that is wide (species as columns) and 
# has environmental variables as the first few columns (default is 9)
# g.vars is the first grouping, default is t.vars=c("blockID","plotID","Scar","Graze","sampling")
# t.vars is the temporal grouping and the default is t.vars=c("blockID","plot,ID","Scar","Graze")
# spatial grouping vars default is s.vars=c("Scar","Graze","Bay","sampling")
#tors default is "temporal" option is "spatial"
# baseval default is "abund" option is "spr"

cvhome<-function(ds,
             nc=9,
             tors="temporal",
             g.vars=c("blockID","plotID","sampling","Scar","Graze"),
             t.vars=c("blockID","plotID","Scar","Graze"),
             s.vars=c("Scar","Graze","Bay","sampling"),
             baseval="abund"){
  t1<-data.frame(ds[,1:nc],spr=specnumber(ds[,c(-1:-nc)]))
  t2<-ds%>%
    pivot_longer(-1:-all_of(nc),names_to = "taxa",values_to="abund")%>%
    group_by(across(all_of(g.vars)))%>%
    summarize(abund=sum(abund))%>%
    ungroup()%>%
    left_join(t1)
  if(tors=="temporal"){
    return(t2%>%
             ungroup()%>%
             group_by(across(all_of(t.vars)))%>%
             summarize(temp.cva=cv(abund),
                temp.cvs=cv(spr)))
  }
  if(tors=="spatial"){
    return(t2%>%
             ungroup()%>%
             group_by(across(all_of(s.vars)))%>%
             summarize(space.cva=cv(abund),
                       space.cvs=cv(spr)))
  }
}
