
## load sim with best tuning parameters
Names<-read.csv("../../../NCC2_09032022.csv",header = F)
NameKey<-data.frame(Name=Names$V4,FG=paste0("X",c(1:102)))


Test<-read.csv("NCC2_09032022_001_1_0003_28-Mar-2023_SR1.csv",header=F)
Test$run<-c(1:nrow(Test))
TestL<-pivot_longer(Test,names_to = "FG",cols = 1:102)
TestL$FG<-gsub("V","X",TestL$FG)
TestL<-merge(TestL,NameKey,by="FG")

# DR<-readMat(FileName)
# Test<-data.frame(DR$re.Y[,,subregion[j]])
# Test$run<-c(1:nrow(Test))
# TestL<-pivot_longer(Test,names_to = "FG",cols = 1:102)
# TestL<-merge(TestL,NameKey,by="FG")
QB<-read.csv("../QB.csv")
## go from yearly to daily:
QB$QB<-QB$QB/365

Sim<-merge(TestL,QB,by="Name")
# head(Sim)
# rm(Test,TestL,Names,NameKey)

SimSum<-Sim %>% group_by(Name) %>% summarise(
  start=value[run=1],
  end=mean(value[run=c(max(run):(max(run)-365))]), ## last year
  # end=value[run=max(run)], ## last event
  prop=end/start
)

SimSum[SimSum$prop>1.5,]

mean(Sim$value[run=c(max(Sim$run):(max(Sim$run)-365))])

# ggplot(data=Sim,aes(x=run,y=value,group=Name))+
#   geom_line()

SimProp<-Sim %>% group_by(Name) %>% summarise(
  start=value[run==1],
  prop=value/start,
  run=run
)

mean<-data.frame(Name="Average",start=NA,SimProp %>% group_by(run) %>% summarise(prop=mean(prop)))

Sim$value[Sim$Name=="PlnkF2biii: Chinook subyearling Fa early"&Sim$run==1]

Extinctions<-paste(unique(SimProp$Name[SimProp$prop<0.01]),collapse="\n")

# View(SimProp %>% select(Name,start) %>% unique())
ggplot(data=SimProp,aes(x=run/365,y=prop,group=Name))+
  geom_line()+
  geom_line(data=mean,aes(x=run/365,y=prop),color="blue",size=3)+
  # geom_line(data=SimProp[SimProp$prop<0.01,],
  #           aes(x=run,y=prop,group=Name),color="red")+
  geom_line(data=SimProp[SimProp$Name%in%SimProp$Name[SimProp$prop<0.01],],
            aes(x=run,y=prop,group=Name),color="red")+
  # geom_text(x=.9*max(SimProp$run),y=.9*max(SimProp$prop),label=Extinctions)
  labs(title=paste("Extinctions =",Extinctions))
# coord_cartesian(xlim=c(0,100))

SimProp$year<-SimProp$run/365
years=150
SimProp$prop[SimProp$year==years]

nonFlat<-SimProp$Name[SimProp$year==(years-20)][which(abs(SimProp$prop[SimProp$year==(years-20)]/SimProp$prop[SimProp$year==years]-1)>0.05)]
SimProp$flat=T
SimProp$flat[SimProp$Name%in%nonFlat]<-F

SimProp$Weight = ifelse(SimProp$run==1, 10000, 1)

SimPropTest<-SimProp[SimProp$run<=3650,]

pdf("150year_AVG_CUTI/Equilibrium_150Year.pdf",height=7,width=9)
ggplot(data=SimProp,aes(x=run/365,y=prop,group=Name,weight=Weight))+
  geom_smooth(se=F,aes(color=flat,group=Name,weight=Weight))+
  scale_color_manual(values=c("black","blue"))+
  labs(x="Burn-in time (years)\n",
       y="\nAbundance relative to starting")+
  theme_minimal()+
  theme(axis.line = element_line(),
        legend.position = "none")
dev.off()

write.csv(
data.frame(
  Name=SimProp$Name[SimProp$year==(years-20)],
  Change=(SimProp$prop[SimProp$year==(years-20)]/SimProp$prop[SimProp$year==years]-1)
),
"../../Tables/TableS6_150-year_simulation_percChange.csv",
row.names=F
)
