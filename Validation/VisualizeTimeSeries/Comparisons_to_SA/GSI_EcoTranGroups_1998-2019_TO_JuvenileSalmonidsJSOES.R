JS<-read.csv("Comparisons_to_SA/GSI_EcoTranGroups_1998-2019.csv")

head(JS)
library(dplyr)
JS<-JS %>% select(group,Fish.Wt..gm.,Year,Month) %>% rename(Wt=Fish.Wt..gm.)
JS<-JS[-which(JS$group==1&JS$Month=="September"),]
JS$group<-ifelse(JS$group=="3&4"&JS$Month=="September","4",
       ifelse(JS$group=="3&4"&JS$Month!="September","3",JS$group))
JS<-JS %>% group_by(group,Year,Month) %>% summarise(Biomass=sum(Wt,na.rm=T)/1000)


write.csv(JS,"Comparisons_to_SA/JuvenileSalmonidsJSOES.csv",row.names=F)
