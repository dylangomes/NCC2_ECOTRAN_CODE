## load libraries ####
library(tidyr)
library(dplyr)
library(ggplot2)
# library(R.matlab) ## no longer load directly from matlab output

## set paths ####
setwd("/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE/Validation/VisualizeTimeSeries")
#Test<-read.csv("NCC2_09032022_001_1_0004_28-Mar-2023_SR1.csv",header=F) ## ecotran model run path
Test <- read.csv("/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/ECOTRANprojects/2D_upwelling_WA_12jul23/Output/temp/NCC2_09032022_008_1_0001_17-Jul-2023_SR1.csv", header=F)
SA<-read.csv("Comparisons_to_SA/SA_quick.csv")  ## stock assessments (or other) file path
Names<-read.csv("../../NCC2_09032022.csv",header = F) ## food web model file path

## set paths ####

NameKey<-data.frame(Name=Names$V4,FG=paste0("X",c(1:102)))

# subregion=c(1:5)
# i=1

# for(i in 1:length(subregion)){
# Test<-data.frame(DR$re.Y[,,subregion[i]])
Test$run<-c(1:nrow(Test))
TestL<-pivot_longer(Test,names_to = "FG",cols = 1:102)
TestL$FG<-gsub("V","X",TestL$FG)
TestL<-merge(TestL,NameKey,by="FG")

## remove large object if don't need other subregions
# rm(DR)

QB<-read.csv("QB.csv")
## go from yearly to daily:
QB$QB<-QB$QB/365

Sim<-merge(TestL,QB,by="Name")
# head(Sim)
# rm(Test,TestL,Names,NameKey)

## Get units back to original
## not sure which to use

CalcBD<-function(method="phytoplankton"){
    if(method=="phytoplankton"){
        ## phytoplankton conversions
        Sim$BiomassDensity = Sim$value/Sim$QB*15*(1/1e-6)*(1/1000)*(106/16)*(12/1)*(1/0.5)*(1/0.2)*(1/1e6)  
    }else{
        ## fish conversions
        Sim$BiomassDensity = Sim$value/Sim$QB*15*(1/1e-6)*(1/1000)*(106/16)*(12/1)*(1/0.65)*(1/0.3)*(1/1e6)  
        # Sim$BiomassDensity = Sim$value/Sim$QB*15*(1/1e-6)*(1/1000)*(106/16)*     (12.0107/1)*(1/0.036572)*(1/0.080178)* (1/1e6) *0.13
    }
    return(Sim)
}

Sim.p<-CalcBD()
Sim.f<-CalcBD(method="fish")

## create plotting function
GGplot<-function(Dat,Keep=unique(Dat$Name),max=50){
    p<-ggplot(Dat[Dat$Name%in%Keep,],
              aes(x=run/365,y=BiomassDensity,group=Name,color=Name))+
        geom_line()+
        labs(title=paste("subregion =",subregion[i]),
             x="Years (starts 1998)",
             y=expression("Biomass density mt/km"^2))+
        geom_point(aes(x=0, y=Biomass),size=2)+
        coord_cartesian(xlim=c(0,max))
    
    return(p)
}

# Select some groups for plotting

#### Primary Production ####

Phyto<-TestL[grep("phyto",TestL$Name),]
# head(Phyto)
temp<-data.frame(Name=unique(Phyto$Name),Biomass=c(45.75326,6.230052))
Phyto<-merge(Phyto,temp,by="Name")

## add original data on top:
PP<-read.csv("Comparisons_to_SA/Phytoplankton_AveragedByRegion_2002-2021.csv")

# max(PP$run)/365 ## check how years the run is
## create a "run" column to merge by, which is by day, starting in 1998 
PP$run<-(PP$Year-1998)*365+((PP$Month-1)*31+31)
PP<-PP[PP$run>=0,]
# head(PP)

scale1=10 # to get the first axis in correct place (relative to starting biomass)
scale2=1 # to get the second axis in correct place

## group PP groups
PhytoAgg<-Phyto %>% group_by(run) %>% summarise(Biomass=sum(Biomass),value=sum(value))
#png("Figures/PP_timeseriesFigure.png",height=9,width=12,units="in",res=600)
p <- ggplot(data=PhytoAgg)+
    geom_line(aes(x=run/365,y=value*scale1,color="EcoTran"))+
    geom_point(aes(x=-0.5, y=Biomass),size=3,shape=17,show.legend = F)+
    # geom_errorbar(aes(x=-0.5, ymin=Biomass-(0.069*Biomass),ymax=Biomass+(0.069*Biomass)),width=0.2)+
    labs( x="Year",
          # y=expression("Ecosystem model biomass density mt/km"^2))+
          y=expression("Phytoplankton biomass density mt/km"^2))+
    
    geom_line(data=PP,aes(x=run/365,y=PP.mtkm2/scale2,color="VGPM"),alpha=0.5)+
    
    scale_x_continuous(breaks=c(seq(0,20,by=2)),labels=c(seq(1998,2018,by=2)))+
    theme_minimal()+
    theme(
        axis.line =  element_line(),
        legend.position = c(0.92,0.95),
        legend.background = element_rect(fill="white",color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        axis.text.x = element_text(angle=90,size=20,vjust = 0.5),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=24),
        plot.title = element_text(size=24,hjust = 0.5),
        axis.ticks.x = element_line(size=1), 
        axis.ticks.length = unit(5, "pt"),
        plot.margin = margin(.1,.1,1,1, "cm")
    )+
    coord_cartesian(xlim=c(0,20))+
    scale_color_manual(values = c("VGPM"="blue","EcoTran"="black"),
                       breaks=c("VGPM","EcoTran"))+ 
    guides(colour = guide_legend(override.aes = list(alpha = c(.5,1),linewidth=c(1))))
ggsave("Figures/PP_timeseriesFigure.pdf", width=12, height=9)

## Plot other groups

# select=0 is even years, select=1 is odd years for easier visualization
PlotVal<-function(ValData,ETData,ValCol,TS=2,select=0,SPAN=.3){
    ETData$BiomassDensity<-scale(ETData$BiomassDensity)
    ValData$BIOM.s<-scale(ValData[,which(names(ValData)==ValCol)])
    # max(Hake$run)/365
    ## create a "run" column to merge by, which is by day, starting in 1998 
    ValData$run<-(ValData$Year-start)*365
    ValData<-ValData[ValData$run>=0,]
    ## which repetition in the timeseries to plot
    
    ggplot(data=ETData)+
        geom_line(aes(x=run/365,y=BiomassDensity))+
        labs(title=paste("group =",unique(ETData$Name)),
             x=paste0("Years (",start,"-",end,") X3"),
             y=expression("Standardized relative biomass"))+
        geom_point(data=ValData,aes(x=(run/365)+repeats*(TS-1),y=BIOM.s),
                   size=2,color="blue")+
        geom_smooth(data=ValData,aes(x=(run/365)+repeats*(TS-1),y=BIOM.s),
                    method="loess",se=F,span=SPAN)+
        geom_vline(xintercept = c(repeats,repeats*2))+
        scale_x_continuous(breaks=c(0:(length(Label)-1))[which(Label%%2==select)],labels=Label[which(Label%%2==select)])+
        theme(axis.text.x = element_text(angle=90))+
        annotate("rect", xmin=0, xmax=repeats*(TS-1), ymin=-Inf, ymax=Inf, alpha=0.4, fill="black")+
        annotate("text",label = "Burn-in", x = repeats*(TS-1)/2, y = 1.25, size = unit(8, "pt"),color="black")+
        annotate("rect", xmin=repeats*(TS), xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.4, fill="green")+
        annotate("text",label = "Stability", x = (repeats*(TS)+repeats*(TS+1))/2, y = 1.25, size = unit(8, "pt"),color="forestgreen")
}

PlotView<-function(ValData,ETData,ValCol,TS=2,SPAN=.3,NAME="default"){
    ETData<-ETData[which((ETData$run/365)>=0 & (ETData$run/365)<=33),]
    ETData$BiomassDensity<-scale(ETData$BiomassDensity)
    ValData$BIOM.s<-scale(ValData[,which(names(ValData)==ValCol)])
    # max(Hake$run)/365
    ## create a "run" column to merge by, which is by day, starting in 1998 
    ValData$run<-(ValData$Year-start)*365
    ValData<-ValData[ValData$run>=0,]
    ## which repetition in the timeseries to plot
    NAME<-ifelse(NAME=="default",paste(gsub(".*: ","",unique(ETData$Name))),
                 NAME)
    
    ggplot(data=ETData)+
        geom_line(aes(x=run/365,y=BiomassDensity,color="EcoTran"))+
        labs(title=NAME,
             x=paste0("\nYears (",start,"-",end,")"),
             y=expression("Standardized relative biomass\n"))+
        geom_point(data=ValData,aes(x=(run/365),y=BIOM.s,color="Independent"),
                   size=2)+
        geom_smooth(data=ValData,aes(x=(run/365),y=BIOM.s),
                    method="loess",se=F,span=SPAN)+
        scale_x_continuous(breaks=c(0:(length(Label)/3-1)),labels=Label[1:(length(Label)/3)])+
        theme_minimal()+
        theme(
            legend.title = element_blank(),
            axis.line =  element_line(),
            legend.text = element_text(size=18),
            axis.text.x = element_text(angle=90,size=20,vjust = 0.5),
            axis.text.y = element_text(size=20),
            axis.title = element_text(size=24),
            plot.title = element_text(size=24,hjust = 0.5),
            axis.ticks.x = element_line(size=1), 
            axis.ticks.length = unit(5, "pt"),
            plot.margin = margin(.1,.1,1,1, "cm")
        )+  
        scale_color_manual(values = c("Independent"="blue","EcoTran"="black"),
                           breaks=c("Independent","EcoTran"))+
        guides(colour = guide_legend(override.aes = list(size = c(2,0),
                                                         linewidth=c(1))))
}
## set up timeseries years
start=1988
end=2021
repeats=(end+1)-start

Label=rep(c(start:(end)),3)
Label[which(Label%%2==0)]<-""


## FIG 11 ####

#### SeaNettle ##

SeaNettle<-Sim.p[grep("large jell",Sim.p$Name),]
SN<-SA[which(SA$Name=="SeaNettle"),] 


# pdf("SeaNettle_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=SN,ETData=SeaNettle,ValCol="Biomass")
# dev.off()

p <- PlotView(ValData=SN,ETData=SeaNettle,ValCol="Biomass",NAME="Sea nettles")+
    theme(legend.position = c(0.9,0.95),
          legend.background = element_rect(fill="white",color="black"),
    )
ggsave("Figures/SeaNettle_timeseriesFigure.png", width=12, height=9, dpi=600)

#### squid ##

squid<-Sim.p[grep("small cephalopod",Sim.p$Name),]
MS<-SA[which(SA$Name=="Market Squid"),] 


# pdf("Market Squid_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=MS,ETData=squid,ValCol="Biomass")
# dev.off()

p <- PlotView(ValData=MS,ETData=squid,ValCol="Biomass",NAME="Market squid")+
    theme(legend.position = "none")
ggsave("Figures/Market Squid_timeseriesFigure.png", width=12, height=9, dpi=600)


#### Dungeness ##

Dungy<-Sim.p[grep("Dungeness",Sim.p$Name),]
DN<-read.csv("Comparisons_to_SA/crab_model_results_2020422.csv")
DN<-DN[which(DN$area!="Central CA"),] # remove Cen CA
DN$Year<-DN$season

# pdf("Figures/Dungeness_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=DN,ETData=Dungy,ValCol="mean_est_thousands_mt")
# dev.off()

p <- PlotView(ValData=DN,ETData=Dungy,
         ValCol="mean_est_thousands_mt")+
    theme(legend.position = "none")
ggsave("Figures/Dungeness_timeseriesFigure.png", width=12, height=9, dpi=600)

### FIG 12 ####

#### sardine ##

sardine<-Sim.p[grep("sardine",Sim.p$Name),]
Sar<-SA[which(SA$Name=="sardine"),] 

# pdf("sardine_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=Sar,ETData=sardine,ValCol="Biomass")
# dev.off()

p <- PlotView(ValData=Sar,ETData=sardine,ValCol="Biomass")+
    theme(legend.position = c(0.9,0.95),
          legend.background = element_rect(fill="white",color="black"))
ggsave("Figures/sardine_timeseriesFigure.png", width=12, height=9, dpi=600)

#### Anchovy ##

anchovy<-Sim.p[grep("anchovy",Sim.p$Name),]
AN<-SA[which(SA$Name=="anchovy"),] 

# 
# pdf("anchovy_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=AN,ETData=anchovy,ValCol="Biomass")
# dev.off()

p <- PlotView(ValData=AN,ETData=anchovy,ValCol="Biomass")+
    theme(legend.position = "none")
ggsave("Figures/anchovy_timeseriesFigure.png", width=12, height=9, dpi=600)

#### jack ##

jack<-Sim.p[grep("jack",Sim.p$Name),]
JM<-SA[which(SA$Name=="jack mackerel"),] 


# pdf("jack mackerel_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=JM,ETData=jack,ValCol="Biomass",SPAN=1)
# dev.off()

p <- PlotView(ValData=JM,ETData=jack,ValCol="Biomass",SPAN=1)+
    theme(legend.position = "none")
ggsave("Figures/jack mackerel_timeseriesFigure.png", width=12, height=9, dpi=600)

#### chub ##

chub<-Sim.p[grep("Pacific mackerel",Sim.p$Name),]
CM<-SA[which(SA$Name=="chub mackerel"),] 

# 
# pdf("chub mackerel_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=CM,ETData=chub,ValCol="Biomass",SPAN=.75)
# dev.off()

p <- PlotView(ValData=CM,ETData=chub,ValCol="Biomass",SPAN=.75)+
    theme(legend.position = "none")
ggsave("Figures/chub mackerel_timeseriesFigure.png", width=12, height=9, dpi=600)


#### SRKW ####

SRKW<-Sim.p[grep("Southern resident killer whales",Sim.p$Name),]
SR<-SA[which(SA$Name=="SRKW"),] 


# pdf("SRKW_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=SR,ETData=SRKW,ValCol="Biomass")
# dev.off()

p <- PlotView(ValData=SR,ETData=SRKW,ValCol="Biomass")+
    theme(legend.position = "none")
ggsave("Figures/SRKW_timeseriesFigure.png", width=12, height=9, dpi=600)

#### Murre ####

Murre<-Sim.p[grep("murre",Sim.p$Name),]
CM<-SA[which(SA$Name=="Murre"),] 

# pdf("Murre_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=CM,ETData=Murre,ValCol="Biomass",SPAN=.7)
# dev.off()

p <- PlotView(ValData=CM,ETData=Murre,ValCol="Biomass",SPAN=.7)+
    theme(legend.position = "none")
ggsave("Figures/Murre_timeseriesFigure.png", width=12, height=9, dpi=600)


#### Shearwater ####

Shearwater<-Sim.p[grep("shearwater",Sim.p$Name),]
SS<-SA[which(SA$Name=="Shearwater"),] 


# pdf("Shearwater_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=SS,ETData=Shearwater,ValCol="Biomass",SPAN=1)
# dev.off()

p <- PlotView(ValData=SS,ETData=Shearwater,ValCol="Biomass",SPAN=1)+
    theme(legend.position = "none")
ggsave("Figures/Shearwater_timeseriesFigure.png", width=12, height=9, dpi=600)

#### Humpback whale ####

Humpback<-Sim.p[grep("baleen",Sim.p$Name),]
HW<-SA[which(SA$Name=="Humpback whale"),] 


# pdf("Humpback_timeseriesValidation.pdf",height=9,width=12)
# PlotVal(ValData=HW,ETData=Humpback,ValCol="Biomass",SPAN=.7)
# dev.off()

p <- PlotView(ValData=HW,ETData=Humpback,ValCol="Biomass",SPAN=.7)+
    theme(legend.position = "none")
ggsave("Figures/Humpback_timeseriesFigure.png", width=12, height=9, dpi=600)

#### Juvenile salmonids ####
JS<-read.csv("Comparisons_to_SA/JuvenileSalmonidsJSOES.csv")
Groups<-c("PlnkF2bi: Chinook yearling Sp",
          "PlnkF2bii: Chinook yearling Fa",
          "PlnkF2biii: Chinook subyearling Fa early",
          "PlnkF2biv: Chinook subyearling Fa late")

Names<-c( "yearling spring Chinook",
          "yearling fall Chinook",
          "subyearling fall Chinook (early ocean entry)",
          "subyearling fall Chinook (late ocean entry)")



for(i in 1:length(Groups)){
    JS.ET<-Sim.p[grep(Groups[i],Sim.p$Name),]
    jsoes<-JS[which(JS$group==i),] 
    
    
    # pdf(paste0(Names[i],"_timeseriesValidation.pdf"),height=9,width=12)
    # print(PlotVal(ValData=jsoes,ETData=JS.ET,ValCol="Biomass",SPAN=.3))
    # dev.off()
    
    p <- PlotView(ValData=jsoes,ETData=JS.ET,ValCol="Biomass",SPAN=.5,NAME=Names[i])
    ggsave(paste0(Names[i],"_timeseriesFigure.png"), width=12, height=9, dpi=600)
    
}
