###################################
#
# R code for analysis of 
# Fleischer et al. Nature Geoscience 2019
# AmazonFACE model ensemble project

# AmazonFACE ecosystem model ensemble analysis
# K Fleischer
# 
# 8 May 2019
#
###################################
rm(list=ls())

#load required libraries
library(ggplot2)
library(reshape2)
library(lattice)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(grid)
library(dplyr)
library(extrafont)
#loadfonts(device = "pdf", quiet = FALSE)

#upload model output data
dfa<-read.table("amaface_CO2_biolyear.csv")

#ordering factors of model codes, name
dfa$MOD<-factor(dfa$MOD, levels = c("CAB","POP","ELM","ECA","GDP","ORC","LPJ","PON","OCN","GDA","JUL","ED2","FAT","INL"))
model_names<-c("CABLE","CABLE-POP","ELM-CTC","ELM-ECA","GDAY","ORCHIDEE","LPJ-GUESS","CABLE-POP(CN)","O-CN","GDAY(CN)","JULES","ED2","ELM-FATES","InLand")

#create colour scale for models
cols=data.frame(MOD=c(levels(dfa$MOD)),COL=c(brewer.pal(7,"Greens")[7:2],brewer.pal(6,"Blues")[2:6],brewer.pal(4,"Greys")[c(2:4)]),stringsAsFactors = F)

#function to calculate rel change of variables over initial, end time period (+-2 years) and complete period, 
relgr<-function(frame, init, end) {

ind<-frame$YEARbio%in%init
ind2<-frame$YEARbio%in%seq(end-2,end+2,1)
ind3<-frame$YEARbio%in%seq(init,end,1)
  
outinit<-aggregate(frame[ind,-c(1:4)],by=list(MOD=frame$MOD[ind]),FUN=mean,na.rm=T)    
outend<-aggregate(frame[ind2,-c(1:4)],by=list(MOD=frame$MOD[ind2]),FUN=mean,na.rm=T)
outall<-aggregate(frame[ind3,-c(1:4)],by=list(MOD=frame$MOD[ind3]),FUN=mean,na.rm=T)
outch<-data.frame(cbind("MOD"=outinit$MOD,(outend[,-1]-outinit[,-1])))
#outch<-data.frame(cbind("MOD"=outinit$MOD,(1-(outend[,-1]/outinit[,-1]))*-100))
#confusing here: what is the change? absolut %init-%end, or 1-(%end/%init) ?

out<-rbind(outinit,outend,outch,outall)
out$PRD<-c(rep("initial",nrow(outinit)),rep("final",nrow(outinit)),rep("change",nrow(outinit)),rep("wholetime",nrow(outinit)))
out<-out[,c("MOD","PRD",names(out)[!names(out) %in% c("PRD","MOD")])]
out$PRD<-as.factor(out$PRD)
out$PRD<-factor(out$PRD, levels = c("initial", "final", "change","wholetime"))
return(out)

}

#take out what we dont need...
#function for calculating additional variables
#addvars<-function(frame) {
  
  #add additonal variables
  frame$CBIO<-rowSums(frame[,c("CL","CFR","CCR","CW")],na.rm=T)
  frame$CLIT<-rowSums(frame[,c("CFLIT","CCLITB")],na.rm=T)
  frame$CROOT<-rowSums(frame[,c("CFR","CCR")],na.rm=T)
  frame$CAG<-rowSums(frame[,c("CW","CL")],na.rm=T)
  frame$CECO<-rowSums(frame[,c("CSOIL","CBIO","CLIT")],na.rm=T)
  #add ROOT? but comparing model which has CR with model that had only CFR does not make sense
  frame$NBIO<-rowSums(frame[,c("NL","NFR","NCR","NW")],na.rm=T)
  frame$NLIT<-rowSums(frame[,c("NFLIT","NCLITB")],na.rm=T)
  frame$NECO<-rowSums(frame[,c("NSOIL","NBIO","NLIT")],na.rm=T)
  
  frame$PBIO<-rowSums(frame[,c("PL","PFR","PCR","PW")],na.rm=T)
  frame$PLIT<-rowSums(frame[,c("PFLIT","PCLITB")],na.rm=T)
  frame$PECO<-rowSums(frame[,c("PSOIL","PBIO","PLIT")],na.rm=T)
  
  frame$CUE<-frame$NPP/frame$GPP
  frame$NUE<-frame$NPP/frame$NUP
  frame$PUE<-frame$NPP/frame$PUP
  
  #add sum of growth, allocation AG and BG
  frame$CG_TOT<-rowSums(frame[,c("CGL","CGW","CGCR","CGFR")],na.rm=T)
  
  frame$CALLOC_L<-frame$CGL/frame$CG_TOT
  frame$CALLOC_W<-frame$CGW/frame$CG_TOT
  frame$CALLOC_CR<-frame$CGCR/frame$CG_TOT
  frame$CALLOC_FR<-frame$CGFR/frame$CG_TOT
  
  frame$CG_AG<-rowSums(frame[,c("CGL","CGW")],na.rm=T)
  frame$CG_BG<-rowSums(frame[,c("CGCR","CGFR")],na.rm=T)
  
  frame$CALLOC_AG<-frame$CG_AG/frame$CG_TOT
  frame$CALLOC_BG<-frame$CG_BG/frame$CG_TOT
  
  #add sum of litter input, AG and BG
  frame$CLITIN_TOT<-rowSums(frame[,c("CLITIN","CWLIN","CCRLIN","CFRLIN")],na.rm=T)
  frame$CLITIN_AG<-rowSums(frame[,c("CLITIN","CWLIN")],na.rm=T)
  frame$CLITIN_BG<-rowSums(frame[,c("CCRLIN","CFRLIN")],na.rm=T)
  frame$CLITIN_AG_frac<-frame$CLITIN_AG/frame$CLITIN_TOT
  frame$CLITIN_BG_frac<-frame$CLITIN_BG/frame$CLITIN_TOT
  
  #calculate residence time in yrs, in DeKauwe its called lifespan and calculated as litter/pool (which is exactly the same...)
  frame$CBIO_RES<-1/(frame$CLITIN_TOT/frame$CBIO)  #residence time is the 1-fraction of litter versus pool size
  frame$CW_RES<-1/(frame$CWLIN/frame$CW)  #residence time is the 1-fraction of litter versus pool size
  frame$CL_RES<-1/(frame$CLITIN/frame$CL)  #residence time is the 1-fraction of litter versus pool size
  frame$CCR_RES<-1/(frame$CCRLIN/frame$CCR)  #residence time is the 1-fraction of litter versus pool size
  frame$CFR_RES<-1/(frame$CFRLIN/frame$CFR)  #residence time is the 1-fraction of litter versus pool size
  frame$CSOIL_RES<-1/(frame$RHET/frame$CSOIL)  #residence time is the 1-fraction of hetero. respiration versus pool size
  
  #calculate increment change of biomass sum(CG,CLITIN) what is the actual change in biomass from year to year)
  frame$CBIO_inc<-frame[,c("CG_TOT")]-frame[,c("CLITIN_TOT")]
  
  #calculate fluxes per biomass
  frame$RAU_CBIO<-frame$RAU/frame$CBIO
  frame$GPP_CBIO<-frame$GPP/frame$CBIO
  frame$GPP_LAI<-frame$GPP/frame$LAI
  frame$RHET_CSOIL<-frame$RHET/frame$CSOIL
  
  #caclulate water use efficiency
  frame$WUE<-frame$GPP/frame$T
  
  #calculate NPP_N
  frame$NPP_N<-frame$NPP/frame$NL 
  frame$NPP_P<-frame$NPP/frame$PL
  
  #calculate stoichiometry
  frame$LeafCN<-frame$CL/frame$NL
  frame$LeafNP<-frame$NL/frame$PL
  frame$LeafCP<-frame$CL/frame$PL
  
  frame$FrootCN<-frame$CFR/frame$NFR
  frame$FrootNP<-frame$NFR/frame$PFR
  frame$FrootCP<-frame$CFR/frame$PFR
  
  frame$CrootCN<-frame$CCR/frame$NCR
  frame$CrootNP<-frame$NCR/frame$PCR
  frame$CrootCP<-frame$CCR/frame$PCR
  
  frame$WoodCN<-frame$CW/frame$NW
  frame$WoodNP<-frame$NW/frame$PW
  frame$WoodCP<-frame$CW/frame$PW
  
  frame$BioCN<-frame$CBIO/frame$NBIO
  frame$BioNP<-frame$NBIO/frame$PBIO
  frame$BioCP<-frame$CBIO/frame$PBIO
  
  frame$LitCN<-frame$CLIT/frame$NLIT
  frame$LitNP<-frame$NLIT/frame$PLIT
  frame$LitCP<-frame$CLIT/frame$PLIT
  
  frame$SoilCN<-frame$CSOIL/frame$NPORG
  frame$SoilNP<-frame$NPORG/frame$PPORG
  frame$SoilCP<-frame$CSOIL/frame$PPORG
  
  for (y in unique(frame$MOD)) {
    for (k in unique(frame$SIM)) {
      
      ind<-frame$MOD==y&frame$SIM==k
      frame$NMIN_c[ind]<-cumsum(frame$NMIN[ind])
      frame$NUP_c[ind]<-cumsum(frame$NUP[ind])
      frame$NFIX_c[ind]<-cumsum(frame$NFIX[ind])
      frame$NLEACH_c[ind]<-cumsum(frame$NLEACH[ind])
      
      frame$PMIN_c[ind]<-cumsum(frame$PMIN[ind])
      frame$PBMIN_c[ind]<-cumsum(frame$PBMIN[ind])
      frame$PUP_c[ind]<-cumsum(frame$PUP[ind])
      frame$PIN_c[ind]<-cumsum(frame$PWEA[ind])+cumsum(frame$PDEP[ind])
      frame$PLEACH_c[ind]<-cumsum(frame$PLEACH[ind])
    }
  }
  
  return(frame)
}
#dfa<-addvars(dfa)

#---------------------------------------------------------------------------------------------------------------
# CALCULATE effects of CO2

#need 2 years after 15 years of eCO2 to calculate average effect of 15 years (1999-2013)
#15 year effect is taken as mean over 13th to 17th year
maxyear<-2013
ind<-dfa$YEARbio<=maxyear+2 
#separate frames per run
run1 <- subset(dfa[ind,],SIM=='AMB_OBS')
run2 <- subset(dfa[ind,],SIM=='ELE_OBS')

#calculate absolute and relative differences between ambient and eCO2 runs
fr7<-data.frame(run1[,1:4],(run2[,-c(1:5)]-run1[,-c(1:5)]))  #from runs to fr we lose column "SIM"
fr7_rel<-data.frame(run1[,1:4],((run2[,-c(1:5)]/run1[,-c(1:5)])-1)*100)

#calculate relative effects for specific time period with funtion relgr()
#here first year versus mean around 15-years, and complete period
fr7_relgr<-relgr(fr7_rel,1999,maxyear)
fr7_gr<-relgr(fr7,1999,maxyear)

#add grouping variables (C, CN, CNP) in newly created frames
group1<-c("INL","ED2","FAT")
group2<-c("LPJ","OCN","GDA","JUL","PON")
group3<-c("CAB","GDP","ELM","ECA","ORC","POP")
fr7_relgr$MODgr[fr7_relgr$MOD%in%group1]<-"C-only"
fr7_relgr$MODgr[fr7_relgr$MOD%in%group2]<-"CN"
fr7_relgr$MODgr[fr7_relgr$MOD%in%group3]<-"CNP"
fr7_relgr$MODgr<-as.factor(fr7_relgr$MODgr)
fr7_gr$MODgr[fr7_relgr$MOD%in%group1]<-"C-only"
fr7_gr$MODgr[fr7_relgr$MOD%in%group2]<-"CN"
fr7_gr$MODgr[fr7_relgr$MOD%in%group3]<-"CNP"
fr7_gr$MODgr<-as.factor(fr7_relgr$MODgr)


# Figure 1: initial and final effects of eCO2 on 
# a) biomass (p3), b) GPP, NPP and CUE (p2), 
# and c) GPP, NPP and CUE, leaf , root, wood C

ind<-fr7_relgr$PRD=="initial"
vars<-c("GPP","NPP","CUE") 
  
data1<-fr7_relgr[ind,c("MODgr",vars)]
data1<-melt(data1)
colnames(data1)[3]<-c("mean")
levels(data1$MODgr)<-c("C models", "CN models", "CNP models")
  
frame<-aggregate(fr7_relgr[ind,vars],by=list(fr7_relgr$MODgr[ind]),FUN=mean,na.rm=T)  
frame<-melt(frame)
sec<-aggregate(fr7_relgr[ind,vars], by=list(fr7_relgr$MODgr[ind]),FUN=sd,na.rm=T)
frame<-cbind(frame,melt(sec,id.vars=c("Group.1"))[,3]); rm(sec)
colnames(frame)<-c("MODgr","variable","mean","sd")
levels(frame$MODgr)<-c("C models", "CN models", "CNP models")
  
ind<-fr7_relgr$PRD=="final"
vars<-c("GPP","NPP","CUE","CL","CFR","CW")
data2<-fr7_relgr[ind,c("MODgr",vars)]
data2<-melt(data2)
colnames(data2)[3]<-c("mean")
levels(data2$MODgr)<-c("C models", "CN models", "CNP models")
  
frame2<-aggregate(fr7_relgr[ind,vars],by=list(fr7_relgr$MODgr[ind]),FUN=mean,na.rm=T)  
frame2<-melt(frame2)
sec<-aggregate(fr7_relgr[ind,vars], by=list(fr7_relgr$MODgr[ind]),FUN=sd,na.rm=T)
frame2<-cbind(frame2,melt(sec,id.vars=c("Group.1"))[,3]); rm(sec)
colnames(frame2)<-c("MODgr","variable","mean","sd")
levels(frame2$MODgr)<-c("C models", "CN models", "CNP models")
  
#!!!!should use year 2014 here to get 15 years of CO2
#but YEARbio for 2014 is already partly 2014/2015
#so we need to devide by 14.5!!!
ind<-fr7_relgr$PRD=="final"
vars<-c("CBIO")
  
data3<-fr7_gr[ind,c("MODgr",vars)]
data3$variable<-rep("x",nrow(data3)); data3$variable<-factor(data3$variable)
colnames(data3)[2]<-c("mean")
levels(data3$MODgr)<-c("C models", "CN models", "CNP models")

frame3<-aggregate(fr7_gr[ind,vars],by=list(fr7_gr$MODgr[ind]),FUN=mean,na.rm=T)  
frame3<-melt(frame3)
sec<-aggregate(fr7_gr[ind,vars], by=list(fr7_gr$MODgr[ind]),FUN=sd,na.rm=T)
frame3<-cbind(frame3,melt(sec,id.vars=c("Group.1"))[,3]); rm(sec)
colnames(frame3)<-c("MODgr","variable","mean","sd")
levels(frame3$MODgr)<-c("C models", "CN models", "CNP models")

p1<-ggplot(data=frame, aes(x=variable, y=mean, fill=MODgr)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8) +
  geom_point(data=data1, stat="identity",colour="grey20",size=0.75,position=position_dodge(width=0.8), show.legend=F) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.4,lwd=0.4, 
                position=position_dodge(.8)) +
  scale_y_continuous(name="%",limits = c(-15,87),breaks = c(0,20,40,60,80)) +
  scale_x_discrete(name="",labels=c("GPP","NPP","CUE")) +
  scale_fill_manual(values=c("Grey40","Steelblue4","Springgreen4"),name="") +
  labs(title = "Initial response") + theme_bw() +
  theme(strip.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.title=element_text(size=7,face="bold",family="Arial"),
        plot.title = element_text(size=7,face="bold",family="Arial"), #face and family does not work with expression!! must be some trick but cant find it
        legend.position = "right",
        legend.key.width=unit(0.3,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.text=element_text(size=7,family="Arial")) #,face="bold"))
p2<-ggplot(data=frame2, aes(x=variable, y=mean, fill=MODgr)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, show.legend=FALSE) +
  geom_point(data=data2,stat="identity",colour="grey20",size=0.75,position=position_dodge(width=0.8),show.legend=F) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.4,lwd=0.4,      
              position=position_dodge(0.8)) +
  scale_y_continuous(name="%",limits = c(-16, 48),breaks = c(0,10,20,30,40)) +
  scale_x_discrete(name="",labels=c("GPP","NPP","CUE","Leaf C","Fine root C","Wood C")) +
  scale_fill_manual(values=c("Grey40","Steelblue4","Springgreen4")) + 
  labs(title = "Final response") +
  theme_bw() +
  theme(strip.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.title=element_text(size=7,face="bold",family="Arial"),
        plot.title = element_text(size=7,face="bold",family="Arial")) #face and family does not work with expression!! must be some trick but cant find it
p3<-ggplot(data=frame3, aes(x=variable, y=mean/14.5, fill=MODgr)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6, show.legend=FALSE) + 
  geom_point(data=data3, stat="identity",colour="grey20",size=0.75,position=position_dodge(width=0.6),show.legend=F) +
  geom_errorbar(aes(ymin=mean/14.5-sd/14.5, ymax=mean/14.5+sd/14.5),width=.4,lwd=0.4,                  
                position=position_dodge(.6)) +
  scale_y_continuous(name=expression(paste("g C ",m^-2," ",yr^-1))) +
  scale_x_discrete(name="",labels="Biomass C") +
  scale_fill_manual(values=c("Grey40","Steelblue4","Springgreen4"),name="") +
  labs(title = "Final response") +
  theme_bw() +
  theme(strip.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.text=element_text(size=7, face="bold",family="Arial",colour="black"),
        axis.title=element_text(size=7,face="bold",family="Arial"),
        plot.title = element_text(size=7,face="bold",family="Arial")) #face and family does not work with expression!! must be some trick but cant find it


# after re-jigging, I cannot align the top with the bottom
ptop<-ggarrange(p3, p1, ncol=2,widths=c(4,6),labels=c("a","b"),font.label=list(size=7,color="black",face="bold",family="Arial"),common.legend = T,legend="bottom")
print(ggarrange(ptop,p2, nrow=2,labels=c("","c"),font.label=list(size=7,color="black",face="bold",family="Arial")))
ggsave("FIG1_Fleischer_etal_NatGeo.pdf",device=cairo_pdf,width=88,height=90,units=c("mm"))


#-----------------------------------------------------------------------------------------------------------------------------
#FIG3: KEY Process responses of CO2 response in CNP models

#fluxes as cumulative difference
#pools as final difference 
#stoichiometry, allocation and nutrient use as whole time difference
frame<-melt(fr7_gr[,names(fr7_gr)!=("MODgr")],id=c("MOD","PRD"))
frame<-frame[frame$PRD=="wholetime"|frame$PRD=="final",]
frame<-droplevels(frame)

#panel with allocation
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="wholetime",]
frame_sub<-frame_sub[frame_sub$variable%in%c("CALLOC_FR","CALLOC_W"),] #"NUE",
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p1<-ggplot(data=frame_sub, aes(x=variable, y=value*100, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="%",breaks=c(-2,0,2), limits = c(-3, 3)) +
  scale_x_discrete(name="", labels=c("Wood","Fine roots")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "C allocation fraction") +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=7.5,face="bold")) 

#panel with stoichiometry AND PUE
#some things seem odd (BioCN for ECA, how come negative if all else positive?: 
#WoodNP for GDAY is negative but remains positive for BioNP, mingkai explained wood effect but does bio response make sense?
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="wholetime",]
frame_sub<-frame_sub[frame_sub$variable%in%c("LeafCP","BioCP","LitCP","PUE"),]  #"SoilCN","SoilNP","BioCN","BioNP",
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("PUE","LitCP","BioCP","LeafCP"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p2<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name=expression(paste("g C g ",P^-1," ",yr^-1),sep=""), breaks=c(0,1000)) + #, limits = c(-250, 250)) +
  scale_x_discrete(name="", labels=c("PUE","Litter C:P","Biomass C:P","Leaf C:P")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "Plant stoichiometry and P use efficiency") +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=7.5,face="bold")) 

##panels with P pools: Ecosystem retention and shifts in P pools 
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("PECO","PBIO","PLEACH_c","PIN_c","PPORG","PPMIN"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("PECO","PPMIN","PPORG","PBIO","PLEACH_c","PIN_c"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p3<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name=expression(paste("g P ",m^-2," ",yr^-1),sep=""),breaks=c(-1,0,1),limits = c(-1.5,1.5)) +
  scale_x_discrete(name="", labels=c("P ecosystem","P mineral soil","P organic soil","P biomass","P leaching","P input")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "P ecosystem retention") +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=7.5,face="bold")) 

frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("PLAB","PSEC","PMIN_c","PBMIN_c","PUP_c"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("PUP_c","PLAB","PSEC","PBMIN_c","PMIN_c"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p4<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name=expression(paste("g P ",m^-2," ",yr^-1)),breaks=c(-1,0,1),limits = c(-1.7,1.7)) +
  scale_x_discrete(name="", labels=c("P uptake","P labile","P secondary","P bioch. min.","P net min.")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "P plant acquisition") +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=7.5,face="bold")) 

##panels with C pools 
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("CL","CW","CFR"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"),labels=c("CABLE","GDAY","ELM-CTC","CABLE-POP","ELM-ECA","ORCHIDEE"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("CW","CFR","CL"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

#devide by 14.5 to get to annual values(see comments above)
p5<-ggplot(data=frame_sub, aes(x=variable, y=value/14.5, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name=expression(paste("g C ",m^-2," ",yr^-1),sep=""), breaks=c(0,100)) + #, limits = c(-250, 250)) +
  scale_x_discrete(name="", labels=c("C wood","C fine root","C leaf")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "Biomass C gain") +
  theme_bw() +
  theme(plot.margin = unit(c(2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=7.5,face="bold")) 

print(ggarrange(p5, p1, p2, p3, p4, nrow=5,align="v", 
                heights=c(6.5,4.55,6.2,8,8),labels=c("a","b","c","d","e"),font.label=list(size=7,color="black",face="bold",family="Arial")))

ggsave("FIG3_Fleischer_etal_NatGEO.pdf",device=cairo_pdf,width=180,height=120,units=c("mm"))


