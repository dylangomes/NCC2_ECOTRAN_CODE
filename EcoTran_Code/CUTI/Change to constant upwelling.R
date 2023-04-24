library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(dplyr) # package for plotting


## Choose dates for individual years
CD<-read.csv("CUTI_daily.csv")
# head(CD)
# length(CD$year)
# 
# Year<-2008
# Index<-which(CD$year==Year)

### DONT ACCIDENTLY OVERWRITE FILES!!!! CREATE A COPY FIRST

# nc_data <- nc_open('C:/Users/dgome/Documents/ECOTRAN_NCC/code_files/CUTI/CUTI_daily.nc',write=T)
nc_data <- nc_open('CUTI_daily - Copy.nc',write=T)
# print(nc_data)

CUTI<-ncvar_get(nc_data,"CUTI")
dim(CUTI) # 17 latitude values by number of days in file
# mean(CUTI)
CUTI[1:4,1:3] ## just look at a sample of huge file

CUTI<-CUTI[,c(Index)]
dim(CUTI) # 17 latitude values by number of days in file
# CUTI[,]<-0 ## zeros don't work
# CUTI[,]<-0.01 ## this is approximately mean(CUTI)/40
# CUTI[,]<-mean(CUTI)/10 ## scale mean down arbitrarily

ncvar_put(nc_data,"CUTI",vals=CUTI)

nc_close(nc_data)

Vis<-read.csv("Upwelling_CUTI_daily.csv")
AVG<-read.csv("CUTI_AVERAGE_daily.csv")
head(Vis)
head(AVG)
Vis<-Vis %>% filter(year<2022)

Vis$yday<-as.POSIXlt(paste0(Vis$year,"-",Vis$month,"-",Vis$day),format="%Y-%m-%d")$yday
# library(plotly)
# ggplot(Vis,aes(x=yday,y=X45N,group=year,color=year))+
#   geom_line()
# 
# p<-ggplot(Vis,aes(x=yday,y=X45N,group=year,color=year))+
#   geom_smooth(se=F)
# p
# ggplotly(p)

Years<-c(1993,1998,2008,2019)
Vis$col<-ifelse(Vis$year%in%Years,Vis$year,"other")
AVG$col<-"AVG"
AVG<-AVG %>% filter(year==1988)
AVG$year <- 2080
Vis<-bind_rows(Vis,AVG)

ggplot(Vis,aes(x=yday,y=X45N,group=year,color=col))+
  geom_line()+
  geom_line(data=Vis[Vis$year%in%Years,])+
  geom_line(data=Vis[Vis$col=="AVG",],size=1)+
  scale_color_manual(values=c("red","black","blue","orange","green","grey"))+
  theme_minimal()
  # coord_cartesian(xlim=c(0,100))

ggplot(Vis,aes(x=yday,y=X45N,group=year,color=col))+
  geom_smooth(se=F)+
  geom_smooth(se=F,data=Vis[Vis$year%in%Years,])+
  scale_color_manual(values=c("red","black","blue","orange","grey"))+
  theme_minimal()
  

AVG<-Vis  %>% group_by(yday) %>% 
  summarise(meanCUTI=mean(X45N))

ggplot(Vis)+
  geom_smooth(se=F,aes(x=yday,y=X45N,group=year,color=year))+
  geom_line(data=AVG,aes(x=yday,y=meanCUTI),color="red",se=F)+
  theme_minimal()+
  theme(axis.line=element_line(),
        legend.position=c(0.7,.3))+
  scale_x_continuous(expand = c(0, 0))+
  labs(y="CUTI",
       x="Day of year",
       color="Year")


## create file of average upwelling across years
AvgSeries<-Vis  %>% dplyr::group_by(yday) %>% 
  dplyr::summarise_at(vars("X31N":"X47N"),.funs=mean)
head(Vis)
head(AvgSeries)
Vismeta<-Vis %>% select(year,month,day,yday)

Avg<-left_join(Vismeta,AvgSeries,by="yday")

## check it

ggplot(Avg,aes(x=yday,y=X45N,group=year,color=year))+
  geom_smooth(se=F)+
  geom_line(se=F,color="red")

write.csv(Avg,"CUTI_AVERAGE_daily.csv",row.names = F)
