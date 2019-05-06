###################################
#
# Data plotting and analysis of preprocessed model output for AmazonFACE modelling project
# run amaface_process_output_all.R before
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

#set working directories
wd<-c('~/SpiderOak Hive/AMAZONAS/FACE-MEI/model_output/all_models')

upload<-T
plot_diagnostic<-F
plot_mainfig<-T
plot_extrafig<-F

maxyear<-2013
model_names<-c("CABLE","CABLE-POP","ELM-CTC","ELM-ECA","GDAY","ORCHIDEE","LPJ-GUESS","CABLE-POP(CN)","O-CN","GDAY(CN)","JULES","ED2","ELM-FATES","InLand")
#uploads data and small post-process steps
if (upload==T) {
  
setwd("~/SpiderOak Hive/AMAZONAS/FACE-MEI/info_for_modellers")
obs<-read.table("observations_to_plot.csv",header=T,sep=",")  
obs$MOD<-"OBS"  
  
setwd(wd)
print("start uploading and preprocessing")
dfa<-readRDS("amaface_allruns_annual_R2_bio.rds")
units<- readRDS('units_frame.rds')

modcod<-c("CAB","POP","PON","LPJ","OCN","GDA","GDP","ECA","ELM","JUL","ORC","INL","ED2","FAT")
dfa$MOD<-factor(dfa$MOD, levels = c("CAB","POP","ELM","ECA","GDP","ORC","LPJ","PON","OCN","GDA","JUL","ED2","FAT","INL"))
cols=data.frame(MOD=c(levels(dfa$MOD)),COL=c(brewer.pal(7,"Greens")[7:2],brewer.pal(6,"Blues")[2:6],brewer.pal(4,"Greys")[c(2:4)]),stringsAsFactors = F)

#add grouping variables (C, CN, CNP, DGVM) and move forward
group1<-c("INL","ED2","FAT")
group2<-c("LPJ","OCN","GDA","JUL","PON")
group3<-c("CAB","GDP","ELM","ECA","ORC","POP")
dfa$MODgr[dfa$MOD%in%group1]<-"C-only"
dfa$MODgr[dfa$MOD%in%group2]<-"CN"
dfa$MODgr[dfa$MOD%in%group3]<-"CNP"
dfa$MODgr<-as.factor(dfa$MODgr)

group1<-c("ED2","LPJ","FAT","POP","PON")
group2<-c("INL","OCN","GDA","JUL","CAB","GDP","ELM","ECA","ORC")
dfa$MODgr2[dfa$MOD%in%group1]<-"DGVM"
dfa$MODgr2[dfa$MOD%in%group2]<-"noDGVM"
dfa$MODgr2<-as.factor(dfa$MODgr2)

dfa<-dfa[,c(1:2,ncol(dfa)-1,ncol(dfa),3:(ncol(dfa)-2))]

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

#function for calculating additional variables
addvars<-function(frame) {
  
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

#function for printing summary statistics (by group)
stat_summ<-function(frame,vars,modgroup) {  
  summ<-aggregate(frame[,vars],by=list(frame[,modgroup]),FUN=mean,na.rm=T)    
  summ<-rbind(summ,aggregate(frame[,vars],by=list(frame$YEARbio>1900),FUN=mean,na.rm=T)) 
  names(summ)<-c("Means",vars)
  summ$Means[nrow(summ)]<-"ALL"
  print(summ)  
  cat("\n")
  
  summ<-aggregate(frame[,vars],by=list(frame[,modgroup]),FUN=sd,na.rm=T)    
  summ<-rbind(summ,aggregate(frame[,vars],by=list(frame$YEARbio>1900),FUN=sd,na.rm=T)) 
  names(summ)<-c("SDs",vars)
  summ$SDs[nrow(summ)]<-"ALL"
  print(summ)  
}

#function for printing, initial, final, percent and actual change for 1 variable (by group)
var_initfin<-function(dat,var) {
  ind1<-dat$SIM=="AMB_OBS"&dfa$YEARbio==1999
  frame<-aggregate(dat[ind1,var],by=list(dat$MODgr[ind1]),FUN=mean,na.rm=T)
  colnames(frame)<-c(var,"init")
  frame$init_sd<-aggregate(dat[ind1,var],by=list(dat$MODgr[ind1]),FUN=sd,na.rm=T)[,2]
  
  ind3<-dat$SIM=="ELE_OBS"&dat$YEARbio==2013
  rel<-(dat[ind3,var]-dat[dat$SIM=="AMB_OBS"&dat$YEARbio==2013,var])
  relp<-((dat[ind3,var]-dat[dat$SIM=="AMB_OBS"&dat$YEARbio==2013,var])/dat[dat$SIM=="AMB_OBS"&dat$YEARbio==2013,var])*100
  
  frame$ac<-aggregate(rel,by=list(dat$MODgr[ind3]),FUN=mean,na.rm=T)[,2]
  frame$ac_sd<-aggregate(rel,by=list(dat$MODgr[ind3]),FUN=sd,na.rm=T)[,2]
  frame$pc<-aggregate(relp,by=list(dat$MODgr[ind3]),FUN=mean,na.rm=T)[,2]
  frame$pc_sd<-aggregate(relp,by=list(dat$MODgr[ind3]),FUN=sd,na.rm=T)[,2]
  
  ind2<-dfa$SIM=="ELE_OBS"&dfa$YEARbio==2013
  frame$fin<-aggregate(dat[ind2,var],by=list(dat$MODgr[ind2]),FUN=mean,na.rm=T)[,2]
  frame$fin_sd<-aggregate(dat[ind2,var],by=list(dat$MODgr[ind2]),FUN=sd,na.rm=T)[,2]
  return(frame)
}

#some fixes for now (need fixing by modelers) (always copy from processing script)
dfa$TAIR<-dfa$TAIR+273.15  #all models provide in celsius now??
dfa$AIRP[dfa$MOD%in%c("OCN","GDP","GDA","ECA","JUL","ELM")]<-dfa$AIRP[dfa$MOD%in%c("OCN","GDP","GDA","ECA","JUL","ELM")]/100 #is in Pa but should be in hPa
#dfa$VPD[dfa$MOD%in%c("CAB","POP")]<-NA #Bernard said this output is irrelevant, not used by CABLE
dfa$PPMIN[dfa$MOD=="ECA"]<-dfa$PLAB[dfa$MOD=="ECA"]+dfa$PSEC[dfa$MOD=="ECA"]+dfa$POCC[dfa$MOD=="ECA"]+dfa$PPAR[dfa$MOD=="ECA"]
dfa$PSOIL[dfa$MOD=="ECA"]<-dfa$PPMIN[dfa$MOD=="ECA"]+dfa$PPORG[dfa$MOD=="ECA"]
dfa$PPMIN[dfa$MOD=="ORC"]<-dfa$PLAB[dfa$MOD=="ORC"]  #+dfa$PSEC[dfa$MOD=="ORC"]+dfa$POCC[dfa$MOD=="ORC"]+dfa$PPAR[dfa$MOD=="ORC"],na.rm=T)
dfa$NFLIT[dfa$MOD=="LPJ"]<-dfa$NFLITA[dfa$MOD=="LPJ"]+dfa$NFLITB[dfa$MOD=="LPJ"] 
#move CR to fine roots in CAB (still check with BP/YP)
dfa$CFR[dfa$MOD=="CAB"]<-dfa$CCR[dfa$MOD=="CAB"] 
dfa$CCR[dfa$MOD=="CAB"]<-NA
dfa$CGFR[dfa$MOD=="CAB"]<-dfa$CGCR[dfa$MOD=="CAB"] 
dfa$CGCR[dfa$MOD=="CAB"]<-NA
#PBMIN flux is difference between PGMIN and PBMIN (not corrected by BP, but check)
dfa$PBMIN[dfa$MOD=="CAB"]<-dfa$PGMIN[dfa$MOD=="CAB"]-dfa$PBMIN[dfa$MOD=="CAB"]  

#write out 20 years
dfa<-dfa[dfa$YEARbio>1998&dfa$YEARbio<2019,]
dfa<-addvars(dfa)

#write out data to be analyes (for external use)
#write.csv(dfa,"amaface_co2runs_annual.csv",row.names = F)
#saveRDS(dfa, "amaface_co2runs_annualbio.rds")

#---------------------------------------------------------------------------------------------------------------
# CALCULATE effects of CO2, dryness, etc

ind<-dfa$YEAR<=maxyear+2 #need min. 2 extra years to calculate average effect of 15 years (1999-2013)

#separate frames per run
run1 <- subset(dfa[ind,],SIM=='AMB_OBS')
run2 <- subset(dfa[ind,],SIM=='ELE_OBS')
# run3 <- subset(dfa[ind,],SIM=='AMB_WET'); # run4 <- subset(dfa[ind,],SIM=='ELE_WET');# run5 <- subset(dfa[ind,],SIM=='AMB_DRY'); # run6 <- subset(dfa[ind,],SIM=='ELE_DRY')

#calculate absolute and relative differences between runs
#we are interested:  7=2-1, 8=3-1, 9=5-1 (effect of CO2, episodic and prolonged drought) 
# and 10=4-3, 11=6-5 (effect of CO2 under no and prolonged water limitation)
#that gives 12=7-10, 13=7-11  (drought CO2 interaction)

fr7<-data.frame(run1[,1:4],(run2[,-c(1:5)]-run1[,-c(1:5)]))  #from runs to fr we lose column "SIM"
fr7_rel<-data.frame(run1[,1:4],((run2[,-c(1:5)]/run1[,-c(1:5)])-1)*100)
# fr8<-data.frame(run1[,1:2],(run3[,-c(1:3)]-run1[,-c(1:3)]))
# fr8_rel<-data.frame(run1[,1:2],((run3[,-c(1:3)]/run1[,-c(1:3)])-1)*100)
# fr9<-data.frame(run1[,1:2],(run5[,-c(1:3)]-run1[,-c(1:3)]))
# fr9_rel<-data.frame(run1[,1:2],((run5[,-c(1:3)]/run1[,-c(1:3)])-1)*100)
# 
# fr10<-data.frame(run3[,1:2],(run4[,-c(1:3)]-run3[,-c(1:3)]))
# fr10_rel<-data.frame(run3[,1:2],((run4[,-c(1:3)]/run3[,-c(1:3)])-1)*100)
# fr11<-data.frame(run5[,1:2],(run6[,-c(1:3)]-run5[,-c(1:3)]))
# fr11_rel<-data.frame(run5[,1:2],((run6[,-c(1:3)]/run5[,-c(1:3)])-1)*100)
# 
# #here only absolute difference of relative effect makes sense, rel would perc of perc
# fr12<-data.frame(fr10[,1:2],(fr7_rel[,-c(1:3)]-fr10_rel[,-c(1:3)]))   #here we lose column "CO2" - this makes sense right?
# fr13<-data.frame(fr11[,1:2],(fr7_rel[,-c(1:3)]-fr11_rel[,-c(1:3)]))

#calculate relative effects for specific time period with funtion relgr()
#here first year versus mean around 15-years, and complete period
fr7_relgr<-relgr(fr7_rel,1999,maxyear)
fr7_gr<-relgr(fr7,1999,maxyear)
#fr13_relgr<-relgr(fr13,1999,c(2011:2015))  #extend time here to test for stabilization of the CO2 effect

saveRDS(fr7, "amaface_co2fr7_annbio.rds")
saveRDS(fr7_rel, "amaface_co2fr7rel_annbio.rds")
}

if (plot_diagnostic==T) {
  
  pdf("Diag_NUP_CO2effect_permodel.pdf",width=10,height=8)
  plotamb<-xyplot(NUP~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylim=c(5,30),ylab=y,xlab=NULL,
                  scales=list(tck=c(-0.5,0),alternating=F),main=paste('ambient NUP'),
                  panel=function(...){
                    panel.xyplot(...)
                    panel.abline(v=xtick,col="grey85")})
  plotres<-xyplot(NUP~YEARbio,data=fr7,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=y,
                  scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('CO2 response NUP'),
                  panel=function(...){
                    panel.xyplot(...)
                    panel.abline(h=0,col='grey40')
                    panel.abline(v=xtick,col="grey85")},
                  key=list(space='top',border=T,columns=5,text=list(levels(run1$MOD)),
                           lines=list(col=cols$COL,lty=1,lwd=2)))
  print(plotamb,split=c(1,1,1,2),more=T)
  print(plotres,split=c(1,2,1,2),more=F)
  dev.off()
  
  
  
  
  # plot CO2 effect of all variables for exploration
  xtick<-seq(2000,maxyear,2)
  ofile <- paste("CO2resp_7_ann_abs",maxyear,"bio.pdf",sep='_')
  pdf(ofile,width=10,height=8)
  
  for (y in names(dfa)[-c(1:5)]) {
    print(y)
    plotamb<-xyplot(get(y)~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=y,xlab=NULL,
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('ambient',y),
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(v=xtick,col="grey85")})
    plotres<-xyplot(get(y)~YEARbio,data=fr7,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=y,
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('CO2 response',y),
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(h=0,col='grey40')
                      panel.abline(v=xtick,col="grey85")},
                    key=list(space='top',border=T,columns=5,text=list(levels(run1$MOD)),
                             lines=list(col=cols$COL,lty=1,lwd=2)))
    print(plotamb,split=c(1,1,1,2),more=T)
    print(plotres,split=c(1,2,1,2),more=F)
  }
  dev.off()
  
  ofile <- paste("CO2resp_7_ann_rel",maxyear,"bio.pdf",sep='_')
  pdf(ofile,width=10,height=8)
  
  for (y in names(dfa)[-c(1:5)]) {
    print(y)
    plotamb<-xyplot(get(y)~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=y,xlab=NULL,
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('ambient',y),
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(v=xtick,col="grey85")})
    plotres<-xyplot(get(y)~YEARbio,data=fr7_rel,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=y,#ylim=c(0,50),
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('CO2 response',y),
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(h=0,col='grey40')
                      panel.abline(v=xtick,col="grey85")},
                    key=list(space='top',border=T,columns=5,text=list(levels(run1$MOD)),
                             lines=list(col=cols$COL,lty=1,lwd=2)))
    print(plotamb,split=c(1,1,1,2),more=T)
    print(plotres,split=c(1,2,1,2),more=F)
  }
  dev.off()
    
  #complete stoichiometry separately again
  frame<-melt(fr7_gr[,names(fr7_gr)!=("MODgr")],id=c("MOD","PRD"))
  frame<-frame[frame$PRD=="final",]
  frame<-droplevels(frame)
  levels(frame$PRD) <- c("15YEARS")
  
  frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP"),]
  frame_sub<-frame_sub[frame_sub$variable%in%c("BioCN","BioNP","CrootCN","CrootNP","FrootCN","FrootNP","LeafCN","LeafNP","WoodCN","WoodNP","SoilCN","SoilNP","LitCN","LitNP"),]
  frame_sub<-droplevels(frame_sub)
  frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")
  
  p3<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    facet_grid(MOD~.) + #scales = "free_y"
    scale_y_continuous(name="CO2 response of stoichiometry") + #, limits = c(-250, 250)) +
    scale_x_discrete(name="") +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    theme(strip.text=element_text(size=12, face="bold"),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=90, vjust=0.25),
          axis.title=element_text(size=14,face="bold"),
          legend.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=10,face="bold"))
  
  pdf("Diag_FIG4_keymechanisms_NP.pdf",width=8,height=6)
  print(p3)
  dev.off()
  
  
  frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP"),]
  frame_sub<-frame_sub[frame_sub$variable%in%c("NBIO","PBIO","NCR","PCR","NFR","PFR","NL","PL","NW","PW"),]
  frame_sub<-droplevels(frame_sub)
  frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")
  
  p3<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    facet_grid(MOD~.) + #scales = "free_y"
    scale_y_continuous(name="CO2 response of N/P content") + #, limits = c(-250, 250)) +
    scale_x_discrete(name="") +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    theme(strip.text=element_text(size=12, face="bold"),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=90, vjust=0.25),
          axis.title=element_text(size=14,face="bold"),
          legend.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=10,face="bold"))
  
  pdf("Diag_FIG4_NPbiomass.pdf",width=8,height=6)
  print(p3)
  dev.off()  
} #end plot_diagnostic()

# main figures for paper
if (plot_mainfig==T) {

# print("calculate evaluation statistics: C, N and P pools")
# use MODgr for C/CN/CNP and MODgr2 for DGVM separation
# vars<-c("CW","CROOT","CL","CSOIL","LAI","CBIO")  
# stat_summ(dfa[dfa$SIM=="AMB_OBS"&dfa$YEARbio==1999,],vars,"MODgr")
# vars<-c("PLAB","PSEC","PPORG","PBIO")  

print("start plotting main figures")

#FIG 2: CO2 effect on key variables and biomass per model and model group
#---------------------------------------------------------------------------------------------------------

#var_initfin() calculates absolut and percent change per group from start to end 
#takes AMB 1999 versus ELE 2013
  
#NOTE that for fluxes the 15 years period (calculated 13 to 17th year is fine), 
#but the final absolute biomass C response need to be after 16 years, year 2014, 
#as pools are updated only the years after  

print(var_initfin(dfa,"CBIO"))
#noteworthy here: C-only start off with lower biomass, 
#so despite much higher percentage change, the final result is the same 
 
#add grouping variables (C, CN, CNP) and move forward (not considered in relgr function..)
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
ggsave("FIG1_CO2effect_keyvars_perGroup_points_NatGeo.pdf",device=cairo_pdf,width=88,height=90,units=c("mm"))

#print out actual numbers for the text, devide by 14.5!
print("annual biomass C change")
(aggregate(fr7_gr[fr7_gr$PRD=="final","CBIO"],by=list(fr7_gr$MODgr[fr7_gr$PRD=="final"]),FUN=mean,na.rm=T)$x)/14.5
(aggregate(fr7_gr[fr7_gr$PRD=="final","CBIO"],by=list(fr7_gr$MODgr[fr7_gr$PRD=="final"]),FUN=sd,na.rm=T)$x)/14.5

#---------------------------------------------------------------------------------------
#Ambient conditions - "evaluation" of initial conditions based on mean of ambient run 

ind<-dfa$SIM=="AMB_OBS"
ann<-aggregate(dfa[ind,c("GPP","LAI","CBIO","CBIO_inc","CUE","NPP","LeafCN","LeafNP","SoilCN","SoilNP","NMIN","PMIN","NPMIN","PLAB","NPORG","PPORG")],by=list(dfa$MOD[ind]),FUN=mean)
names(ann)<-c("Model","GPP*","LAI*","Biomass C*","Bio. C inc.","CUE","NPP","Leaf C:N*","Leaf N:P*","Soil C:N","Soil N:P","N net min.","P net min.","Mineral N","Labile P","Organic N","Organic P")
levels(ann$Model)<-model_names

#dont trust CBIO_inc calculated in dfa frame yet (based on CG_TOT and LITIN...)
#recalculate here on actual CBIO
cbio_inc<-NULL
for (k in levels(dfa$MOD)) {
 cbio_inc<-c(cbio_inc,mean(diff(dfa[ind==T&dfa$MOD==k,c("CBIO")])))
}
ann$`Bio. C inc.`<-cbio_inc

print("range of ambient C sink")
print(cbio_inc)

#convert fluxes and pools to per kg, except bio C increment, leave at g C
ann[,c(2,4,7)]<-ann[,c(2,4,7)]/1000
frame<-melt(ann)

frame$value[frame$variable=="Mineral N"&frame$Model=="ED2"]<-NA
frame$value[frame$variable=="Organic N"&frame$Model=="ED2"]<-NA
frame$value[frame$variable=="N net min."&frame$Model=="ED2"]<-NA
frame$value[frame$variable=="Leaf C:N*"&frame$Model=="ED2"]<-NA
frame$value[frame$variable=="Soil C:N"&frame$Model=="ED2"]<-NA
frame$value[frame$variable=="Leaf N:P*"&frame$Model=="CABLE-POP(CN)"]<-NA
frame$value[frame$variable=="Soil N:P"&frame$Model=="CABLE-POP(CN)"]<-NA
frame$value[frame$variable=="Labile P"&frame$Model=="CABLE-POP(CN)"]<-NA
frame$value[frame$variable=="Organic P"&frame$Model=="CABLE-POP(CN)"]<-NA

frame_obs<-data.frame(variable = levels(frame$variable), Z = c(3.25,5.3,18.8,64,NA,1.29,23,37,17.1,76,NA,NA,17.5,1.6,827.5,NA))

p1<-ggplot(data=frame, aes(x=Model, y=value, fill=Model)) +
  geom_bar(stat="identity",position="identity",width=0.8,colour="black",show.legend=T) +
  facet_wrap(~variable, scales = "free_y") +
  geom_hline(data=frame_obs,aes(yintercept = Z),linetype=2) +
  scale_y_continuous(name="") +
  scale_x_discrete(name="",labels=NULL,breaks=NULL) +
  scale_fill_manual(values=cols$COL) +
  theme_bw() +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

print(p1)
ggsave("Appendix_Fig2A_Ambient.pdf",width=10,height=8)


# Biomass C change, absolute and relative for supplement per model
p1<-ggplot(data=fr7[fr7$YEARbio<=maxyear,],aes(x=YEARbio-1998, y=CBIO/1000, colour=MOD)) + 
  geom_line(size=1.5,show.legend=T) +
  scale_colour_manual(values=cols$COL,guide = guide_legend(override.aes = list(color = "white"))) +
  scale_y_continuous(name="Biomass C change [kg C m-2]") + 
  scale_x_continuous(name="",breaks=c(1,5,10,15)) +
  theme_bw() +
  theme(strip.text=element_text(size=14, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        # legend.title=element_text(size=14,face="bold"),
        # legend.text=element_text(size=12,face="bold"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

p2<-ggplot(data=fr7_rel[fr7_rel$YEARbio<=maxyear,],aes(x=YEARbio-1998, y=CBIO, colour=MOD)) + 
    geom_line(size=1.5,show.legend=T) +
    scale_colour_manual(values=cols$COL, name="Model",labels=model_names) +
    scale_y_continuous(name="Biomass C change [%]") + 
    scale_x_continuous(name="Years of eCO2 [+200ppm]",breaks=c(1,5,10,15)) +
    theme_bw() +
    theme(strip.text=element_text(size=14, face="bold"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14),
          legend.position = "right")

print(ggarrange(p1, p2, nrow=2,align="v",labels=c("a","b")))
ggsave(paste("Appendix_Fig1A_CBiomass_permodel.pdf",sep=""),width=10,height=8)

#print out actual numbers for the text
print("annual biomass C change per model")
a<-fr7[fr7$YEARbio==maxyear,c("MOD","CBIO")]
a$annual<-a[,2]/14.5
a

# FIG2 Suppl: CO2 effect on C allocation (mean and sd), for elevated mean over 3 years and 15 years
#----------------------------------------------------------------------------------------------------

#C allocation:
#CAB only varies with phenology, nothing
#POP initially more leaf, fr not changing at all?
#ELM neglegible changes, little more leaf at expense of fr
#ECA, GDAY, O-CN more FR at expense of leaf/wood
#ECA is changing in AMB run, reducing FR alloc, and ele stays the same...?? 
#ECA changes in wood but cant see line in graph??

frame<-fr7_rel[fr7_rel$YEARbio<=maxyear,c("MOD","YEARbio","CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")]
#frame<-frame[!frame$MOD%in%c("INL"),]
frame<-melt(frame,id.vars = c("MOD","YEARbio"))
frame$variable<-factor(frame$variable,levels=c("CALLOC_L","CALLOC_W","CALLOC_CR","CALLOC_FR"))
levels(frame$variable)<-c("Leaf","Wood","Coarse root","Fine root")
frame$value[frame$MOD%in%c("CAB")][1:45]<-0 #set C allocation change to 0 (L,W,FR, no CR, still old order in dataframe)(its only due to phenology and misleading)
frame$value[frame$MOD%in%c("INL")][1:45]<-0 #set C allocation change to 0 (its only due to phenology and misleading)

levels(frame$MOD)<-model_names

p1<-ggplot(data=frame, aes(x=YEARbio-1998, y=value, col=variable)) +
  geom_line(lwd=1.2) +
  facet_wrap(~MOD,scales = "free_y") +
  scale_y_continuous(name="C allocation change in %") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_colour_manual(values=c("Springgreen4","firebrick4","steelblue4","steelblue2"),name="Biomass \ncompartment") +
  theme_bw() +
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

#combine the line graph above (Co2 effect on alloc over time) with mean alloc in ambient run
ind<-dfa$YEARbio<=maxyear&dfa$SIM=="AMB_OBS"
frame2<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_CR","CALLOC_FR")], by=list(dfa$MOD[ind]),FUN=mean))
frame2<-melt(frame2,id.vars = c("Group.1"))
sec2<-aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_CR","CALLOC_FR")], by=list(dfa$MOD[ind]),FUN=sd)
frame2<-cbind(frame2,melt(sec2,id.vars=c("Group.1"))[,3]); rm(sec2)
colnames(frame2)<-c("MOD","variable","mean","sd")
levels(frame2$variable)<-c("Leaf","Wood","Coarse root","Fine root")
frame2$mean[frame2$MOD%in%c("INL")]<-c(0.3,0.5,NA,0.2) #set fixed allocation in INL manually
frame2$sd[frame2$MOD%in%c("INL")]<-frame2$sd[frame2$MOD%in%c("CAB")] #just to get the SD line
levels(frame2$MOD)<-model_names

p2<-ggplot(data=frame2, aes(x=variable, y=mean*100, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8) + #, show.legend=FALSE) +
  geom_errorbar(aes(ymin=mean*100-sd*100, ymax=mean*100+sd*100),width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_wrap(~MOD) +
  scale_y_continuous(name="Ambient C allocation in %") + 
  scale_x_discrete(name="",labels=NULL) +
  scale_fill_manual(values=c("Springgreen4","firebrick4","steelblue4","steelblue2"),name="Biomass \ncompartment") +
  theme_bw() +
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

print(ggarrange(p2,p1,align="v",nrow=2))
ggsave("Appendix_Fig1C_C_allocation_combined.pdf",width=10,height=8)

#just to check here what the actual allocation fractions are in the ele run, since the relative effect is somewhat misleading
ind<-dfa$YEARbio<=maxyear&dfa$SIM=="AMB_OBS"
frame2<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind]),FUN=mean))
print("range ambient mean FR allocation per model")
max(frame2$CALLOC_FR,na.rm=T); min(frame2$CALLOC_FR,na.rm=T);

ind<-dfa$YEARbio<=maxyear&dfa$SIM=="ELE_OBS"
frame2<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind]),FUN=mean))
print("range elevated mean FR allocation per model")
max(frame2$CALLOC_FR,na.rm=T); min(frame2$CALLOC_FR,na.rm=T);

ind<-dfa$YEARbio<=maxyear&dfa$SIM=="ELE_OBS"
frame2<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind]),FUN=max))
print("range elevated max FR allocation per model")
max(frame2$CALLOC_FR,na.rm=T); min(frame2$CALLOC_FR,na.rm=T);

#Appendix figure: functionality of increased fine root allocation
#is there a saturation of the benefit of ++FR and do the models keep doing FRallocation+ despite this
frame<-fr7[,c("MOD","YEARbio","CALLOC_FR","PUP")]  #fr7$YEARbio<=maxyear
frame$PUP_FRalloc<-frame$PUP/frame$CALLOC_FR
frame<-frame[frame$MOD%in%c("ECA","ORC","GDP"),]

frame<-melt(frame,id.vars = c("MOD","YEARbio"))
#levels(frame$variable)<-c("Fine root allocation","P uptake","P uptake per FR allocation")
levels(frame$MOD)<-model_names

ind<-frame$variable=="PUP"
p1<-ggplot(data=frame[ind,], aes(x=YEARbio-1998, y=value, col=variable)) +
  geom_line(lwd=1.2,show.legend=FALSE) +
  facet_wrap(~MOD,scales = "free_y") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="") +
  scale_colour_manual(values=c("Springgreen4"),name="P uptake") +   #CO2 effect
  theme_bw() +
  ggtitle("(a) CO2 effect on P uptake [g P / m2]") +
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

ind<-frame$variable=="CALLOC_FR"
p2<-ggplot(data=frame[ind,], aes(x=YEARbio-1998, y=value, col=variable)) +
  geom_line(lwd=1.2,show.legend=FALSE) +
  facet_wrap(~MOD,scales = "free_y") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="") +
  scale_colour_manual(values=c("firebrick4"),name="") +
  theme_bw() +
  ggtitle("(b) CO2 effect on fine root allocation [NPP fraction]") +
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

ind<-frame$variable=="PUP_FRalloc"
p3<-ggplot(data=frame[ind,], aes(x=YEARbio-1998, y=value, col=variable)) +
  geom_line(lwd=1.2,show.legend=FALSE) +
  facet_wrap(~MOD,scales = "free_y") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_colour_manual(values=c("steelblue2"),name="") +
  theme_bw() +
  ggtitle("(c) P uptake per fine root allocation under eCO2 (A/B)") +
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

print(ggarrange(p1,p2,p3,align="v",nrow=3))
ggsave("Appendix_PUP_FRC_CO2effect_relative.pdf",width=10,height=8)


# CO2 effect on vegetation C residence times
# check if CBIO_RES is reliable in the models....!
ind<-dfa$YEARbio<=maxyear
frame<-cbind(aggregate(dfa[ind&dfa$SIM=="AMB_OBS",c("CBIO_RES")], by=list(dfa$MOD[ind&dfa$SIM=="AMB_OBS"]),FUN=mean),
             aggregate(dfa[ind&dfa$SIM=="ELE_OBS",c("CBIO_RES")], by=list(dfa$MOD[ind&dfa$SIM=="ELE_OBS"]),FUN=mean)[,2])
names(frame)<-c("MOD","AMB_OBS","ELE_OBS")
frame[,2:3]<-round(frame[,2:3],digits=2)

#testing the effect of turnover changes on the biomass C response 

#Sönke Zaehle: 
#dBIO = NPP - turnover = NPP - tau * BIO
#dBIO/dCO2 = dNPP / dCO2 - dtau/dCO2 * BIO
#So, testing whether dtau/dCO2 = 0, and then testing (assuming  turnover remains tau*BIO)  
#dBIO/dCO2 = dNPP/dCO2 - tau(initial) * BIO 

turn_amb<-dfa[dfa$MOD=="OCN"&dfa$SIM=="AMB_OBS",c("YEARbio","CBIO","CBIO_inc","NPP","CG_TOT","CLITIN_TOT")]
turn_ele<-dfa[dfa$MOD=="OCN"&dfa$SIM=="ELE_OBS",c("YEARbio","CBIO","CBIO_inc","NPP","CG_TOT","CLITIN_TOT")]

#calculate turnover rate (tau)
turn_amb$tau<-turn_amb$CLITIN_TOT/turn_amb$CBIO
turn_ele$tau<-turn_ele$CLITIN_TOT/turn_ele$CBIO
#testing if NPP (CG_TOT) and tau result in dCBIO (CBIO_inc)
turn_amb$dCBIO<-turn_amb$CG_TOT-(turn_amb$tau*turn_amb$CBIO)
turn_ele$dCBIO<-turn_ele$CG_TOT-(turn_ele$tau*turn_ele$CBIO)
#yes, they agree
turn_amb$CBIO_inc==turn_amb$dCBIO
turn_ele$CBIO_inc==turn_ele$dCBIO

#dBIO/dCO2 = dNPP / dCO2 - dtau/dCO2 * BIO 
#adding CO2 effect in one frame (incl. CBIO of ambient run)
turnCO2<-turn_amb[,c("YEARbio","CBIO")]
#calculating dCBIO/dCO2, dNPP/dCO2, dtau/dCO2
turnCO2$dCBIO<-turn_ele$CBIO_inc-turn_amb$CBIO_inc 
turnCO2$dNPP<-turn_ele$NPP-turn_amb$NPP
turnCO2$dtau<-turn_ele$tau-turn_amb$tau
#-- so dtau/dCO2 is not 0 but increasing with eCO2

#now calculating what dBIO/dCO2 would have been IF tau had not changed with eCO2?
#dBIO/dCO2 = dNPP/dCO2 - tau(initial) * BIO 
turnCO2$dCBIO_ft<-turnCO2$dNPP-(turn_amb$tau*turn_ele$CBIO)

#--> Here I must be doing something wrong, dCBIO turns negative, 
#but turnover increases with eCO2 so without this turnover change dCBIo/dCO2 should
write.csv(frame,"biomassC_residence.csv",row.names = F)

#--------------------------------------------------------------------------------------------------
# Suppl on CO2 effect on NPP/GPP/CUE fluxes per model
frame<-melt(fr7_relgr[,names(fr7_relgr)!=("MODgr")],id=c("MOD","PRD"))  #have to do this fix here now sine MODgr was added..
frame<-frame[frame$PRD!="wholetime",]
frame<-droplevels(frame)
levels(frame$PRD) <- c("1st year", "15 years", "CHANGE")

ind<-frame$variable%in%c("GPP","NPP","CUE","CL","CW","CFR")&frame$PRD%in%c("1st year", "15 years")
levels(frame$variable)[which(levels(frame$variable)=="CL")]<-"Leaf C"
levels(frame$variable)[which(levels(frame$variable)=="CW")]<-"Wood C"
levels(frame$variable)[which(levels(frame$variable)=="CFR")]<-"Fine root C"

p1<-ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="Model",labels=model_names) +
  theme_bw() + 
  theme(strip.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

print(p1)
ggsave("Appendix_Fig1B_GPPNPPCUE_permodel.pdf",width=8,height=4)

#--------------------------------------------------------------------------------------------------
#FIG3: KEY Process responses of NUP/NUE and PUP/PUE per model group (YP analysis)
frame<-aggregate(dfa[,c("NPP","GPP","NUP","NUE","PUP","PUE")], by=list(dfa$MODgr,dfa$SIM,dfa$YEARbio),FUN=mean)
names(frame)[1:3]<-c("Group","SIM","Year")
frame<-frame[frame$Year<=maxyear,]
ind<-frame$Group%in%c("CN")
ind2<-frame$Group%in%c("CNP")

p1<-ggplot(data=frame[ind,], aes(x=Year-1998, y=NUP/GPP, colour=SIM)) +
  geom_line(lwd=0.9)+
  geom_point(show.legend=FALSE,size=2) + 
  labs(title = "A   N uptake / GPP (CN models)") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_colour_manual(values=c(cols$COL[11],cols$COL[8]),guide = guide_legend(override.aes = list(color = "white"))) +
  theme_bw() +
  theme(strip.text=element_text(size=12),
        axis.text=element_text(size=12),
        plot.title = element_text(size=14),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

p2<-ggplot(data=frame[ind2,], aes(x=Year-1998, y=NUP/GPP, colour=SIM)) +
  geom_line(lwd=0.9)+
  geom_point(show.legend=FALSE,size=2) + 
  labs(title = "B   N uptake / GPP (CNP models)") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_colour_manual(values=c(cols$COL[11],cols$COL[8]),guide = guide_legend(override.aes = list(color = "white"))) +
  theme_bw() +
  theme(strip.text=element_text(size=12),
        axis.text=element_text(size=12),
        plot.title = element_text(size=14),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

p3<-ggplot(data=frame[ind2,], aes(x=Year-1998, y=PUP/GPP, colour=SIM)) +
  geom_line(lwd=0.9)+
  geom_point(show.legend=FALSE,size=2) + 
  labs(title = "C   P uptake / GPP (CNP models)") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_colour_manual(values=c(cols$COL[11],cols$COL[8]),guide = guide_legend(override.aes = list(color = "white"))) +
  theme_bw() +
  theme(strip.text=element_text(size=12),
        axis.text=element_text(size=12),
        plot.title = element_text(size=14),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

p4<-ggplot(data=frame[ind2,], aes(x=Year-1998, y=PUP/NUP, colour=SIM)) +
  geom_line(lwd=0.9)+
  geom_point(show.legend=FALSE,size=2) + 
  labs(title = "D   P uptake / N uptake (CNP models)") +
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Years of eCO2") +
  scale_shape_manual(values=c(19,1),name="Simulation") +
  scale_colour_manual(values=c(cols$COL[11],cols$COL[8]),name="Simulation",labels=c("Ambient","eCO2")) +
  theme_bw() +
  theme(strip.text=element_text(size=12),
        axis.text=element_text(size=12),
        plot.title = element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

print(ggarrange(p1,p2,p3,p4,heights=c(1,1,1,1),align="v",nrow=4)) 
ggsave("Appendix_Fig3_NUP_PUP_Annual_modelgroup.pdf",width=8,height=8)

# Suppl on CO2 effect on NUE, C:N, PUE, N:P  per mdel
# CO2 effect on key N and P variables per model

frame<-melt(fr7_relgr[,names(fr7_relgr)!=("MODgr")],id=c("MOD","PRD")) #the fix again
frame<-frame[frame$PRD!="wholetime",]
frame<-droplevels(frame)
levels(frame$PRD) <- c("1st YEAR", "15YEARS", "CHANGE")

ind<-frame$variable%in%c("NUP","NUE","LeafCN","WoodCN","BioCN")&frame$PRD%in%c("1st YEAR", "15YEARS")
levels(frame$variable)[which(levels(frame$variable)=="LeafCN")]<-"Leaf C:N"
levels(frame$variable)[which(levels(frame$variable)=="WoodCN")]<-"Wood C:N"
levels(frame$variable)[which(levels(frame$variable)=="BioCN")]<-"Biomass C:N"

p1<-ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=T) + 
  scale_fill_manual(values=cols$COL,guide = guide_legend(override.aes = list(fill = "white", colour = "white"))) +
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
    theme_bw() +
  theme(strip.text=element_text(size=14, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

ind2<-frame$variable%in%c("PUP","PUE","LeafNP","WoodNP","BioNP")&frame$PRD%in%c("1st YEAR", "15YEARS")
levels(frame$variable)[which(levels(frame$variable)=="LeafNP")]<-"Leaf N:P"
levels(frame$variable)[which(levels(frame$variable)=="WoodNP")]<-"Wood N:P"
levels(frame$variable)[which(levels(frame$variable)=="BioNP")]<-"Biomass N:P"
frame$value[frame$variable=="Leaf N:P"&frame$MOD=="PON"]<-NA
frame$value[frame$variable=="Wood N:P"&frame$MOD=="PON"]<-NA
frame$value[frame$variable=="Biomass N:P"&frame$MOD=="PON"]<-NA

p2<-ggplot(data=frame[ind2,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=T) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="Model",labels=model_names) +
  theme_bw() +
  theme(strip.text=element_text(size=14, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14),
        legend.position = "right")

print(ggarrange(p1,p2,heights=c(1,1),align="v",nrow=2)) 
ggsave("Appendix_Fig3A_CO2_N_Pvars_permodel.pdf",width=10,height=8)

#NOW show the actual changes in key N and P values in amb/ele after 15 years
vars<-c("NUE","LeafCN","WoodCN","BioCN","PUE","LeafNP","WoodNP","BioNP")
ind<-dfa$YEARbio>2010&dfa$YEARbio<2016&(dfa$MODgr=="CN"|dfa$MODgr=="CNP")  #same frame as in fr7_gr/relgr... 
fres<-aggregate(dfa[ind,vars], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=mean,na.rm=T)
names(fres)[1:2]<-c("Group","SIM")
fres<-melt(fres)
sec<-aggregate(dfa[ind,vars], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=sd,na.rm=T)
fres<-cbind(fres,melt(sec,id.vars=c("Group.1","Group.2"))[,-c(1:3)]); rm(sec)
colnames(fres)<-c("MOD","SIM","variable","mean","sd")
levels(fres$SIM)<-c("Ambient","eCO2")
levels(fres$variable)<-c("NUE","Leaf C:N","Wood C:N","Bio. C:N","PUE","Leaf N:P","Wood N:P","Bio. N:P")

ind<-fres$variable%in%c("NUE","Leaf C:N","Wood C:N","Bio. C:N")

p1<-ggplot(data=fres[ind==T,],aes(x=MOD,y=mean,fill=SIM))+ #position="dodge",
  geom_bar(stat ="identity", position="dodge",show.legend = F) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_grid(variable~.,scales="free") +
  scale_y_continuous(name="") +
  scale_x_discrete(name="",labels=model_names[1:11]) +
  scale_fill_manual(values=c(cols$COL[11],cols$COL[8])) +
  theme_bw() +
  theme(strip.text=element_text(size=14),
        axis.text=element_text(size=14),
        plot.title = element_text(size=16),
        axis.text.x=element_text(angle=90, vjust=0.25,hjust=1))

ind2<-fres$variable%in%c("PUE","Leaf N:P","Wood N:P","Bio. N:P")&fres$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC")

p2<-ggplot(data=fres[ind2==T,],aes(x=MOD,y=mean,fill=SIM))+ #position="dodge",
  geom_bar(stat ="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_grid(variable~.,scales="free") +
  scale_y_continuous(name="") +
  scale_x_discrete(name="",labels=model_names[1:6]) +
  scale_fill_manual(values=c(cols$COL[11],cols$COL[8]),name="Simulation") +
  theme_bw() +
  theme(strip.text=element_text(size=14),
        axis.text=element_text(size=14),
        plot.title = element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.position = "right",
        axis.text.x=element_text(angle=90, vjust=0.25,hjust=1))

print(ggarrange(p1,p2,nrow=2)) 
ggsave("Appendix_Fig3B_absolute_NUE_CN_PUE_NP_permodel.pdf",width=8,height=12)

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
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="",breaks=c(-2,0,2), limits = c(-3, 3)) +
  scale_x_discrete(name="", labels=c("Wood","Fine roots")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = "C allocation fraction [%]") +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(family="Arial",size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(family="Arial",size=7,face="bold"),
        plot.title=element_text(family="Arial",size=7,face="bold")) 

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
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="", breaks=c(0,1000)) + #, limits = c(-250, 250)) +
  scale_x_discrete(name="", labels=c("PUE","Litter C:P","Biomass C:P","Leaf C:P")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = expression(paste("Plant stoichiometry and P use efficiency [g C g ",P^-1," (",yr^-1,")]"),sep="")) +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(family="Arial",size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(family="Arial",size=7,face="bold"),
        plot.title=element_text(family="Arial",size=7,face="bold")) 

##panels with P pools: Ecosystem retention and shifts in P pools 
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("PECO","PBIO","PLEACH_c","PIN_c","PPORG","PPMIN"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("PECO","PPMIN","PPORG","PBIO","PLEACH_c","PIN_c"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p3<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="",breaks=c(-1,0,1),limits = c(-1.5,1.5)) +
  scale_x_discrete(name="", labels=c("P ecosystem","P mineral soil","P organic soil","P biomass","P leaching","P input")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = expression(paste("P ecosystem retention [g P ",m^-2," (",yr^-1,")]"),sep="")) +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(family="Arial",size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(family="Arial",size=7,face="bold"),
        plot.title=element_text(family="Arial",size=7,face="bold")) 

frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("PLAB","PSEC","PMIN_c","PBMIN_c","PUP_c"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("PUP_c","PLAB","PSEC","PBMIN_c","PMIN_c"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

p4<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="",breaks=c(-1,0,1),limits = c(-1.7,1.7)) +
  scale_x_discrete(name="", labels=c("P uptake","P labile","P secondary","P bioch. min.","P net min.")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  labs(title = expression(paste("P plant acquisition [g P ",m^-2," (",yr^-1,")]"),sep="")) +
  theme_bw() +
  theme(plot.margin = unit(c(-2,2,2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(family="Arial",size=7,face="bold"),
        axis.text.x=element_text(angle=0),
        axis.title=element_text(family="Arial",size=7,face="bold"),
        plot.title=element_text(family="Arial",size=7,face="bold")) 

##panels with C pools 
frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
frame_sub<-frame_sub[frame_sub$variable%in%c("CL","CW","CFR"),] 
frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"),labels=c("CABLE","GDAY","ELM-CTC","CABLE-POP","ELM-ECA","ORCHIDEE"))
frame_sub<-droplevels(frame_sub)
frame_sub$variable<-factor(frame_sub$variable,levels=c("CW","CFR","CL"))
frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")

#devide by 14.5 to get to annual values(see comments above)
p5<-ggplot(data=frame_sub, aes(x=variable, y=value/14.5, fill=sign)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
  coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
  scale_y_continuous(name="", breaks=c(0,100)) + #, limits = c(-250, 250)) +
  #scale_x_discrete(name="", labels=c("C wood","C fine root","C leaf")) +
  scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
  #labs(title = expression(paste("C gain [g C ",m^-2," ",yr^-1,"]"),sep="")) +
  theme_bw() +
  theme(plot.margin = unit(c(2,2,-2,2), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(family="Arial",size=7,face="bold"),
        axis.text.x=element_text(angle=0,family="Arial"),
        axis.title=element_text(family="Arial",size=7,face="bold"),
        plot.title=element_text(family="Arial",size=7,face="bold")) 

print(ggarrange(p5, p1, p2, p3, p4, nrow=5,align="v", 
                heights=c(6.5,4.5,5.5,8,8),labels=c("a","b","c","d","e"),font.label=list(size=7,color="black",face="bold",family="Arial")))
ggsave("FIG3_keymechanisms_NP_NatGEO.pdf",device=cairo_pdf,width=180,height=120,units=c("mm"))

} #end if plotmainfig

extrafont::loadfonts(device = "pdf")
loadfonts(device = "pdf", quiet = FALSE)

#--------------------------------------------------------------------------------------------------
# all extra figures
if (plot_extrafig==T) {

  #FIG3: KEY Process responses of CO2 response in CNP models : short version for presentations
  #fluxes as cumulative difference
  #pools as final difference 
  #stoichiometry, allocation and nutrient use as whole time difference
  frame<-melt(fr7_gr[,names(fr7_gr)!=("MODgr")],id=c("MOD","PRD"))
  frame<-frame[frame$PRD=="wholetime"|frame$PRD=="final",]
  frame<-droplevels(frame)
  
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
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
    scale_y_continuous(name="", breaks=c(0,1000)) + #, limits = c(-250, 250)) +
    scale_x_discrete(name="", labels=c("PUE","Litter C:P","Biomass C:P","Leaf C:P")) +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    labs(title = "Plant stoichiometry and P use efficiency") +
    theme_bw() +
    theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=0),
          axis.title=element_text(size=14,face="bold"))
  
  ##panels with P pools: Ecosystem retention and shifts in P pools 
  frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
  frame_sub<-frame_sub[frame_sub$variable%in%c("PBIO","PPORG","PPMIN"),] 
  frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
  frame_sub<-droplevels(frame_sub)
  frame_sub$variable<-factor(frame_sub$variable,levels=c("PPMIN","PPORG","PBIO"))
  frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")
  
  p3<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
    scale_y_continuous(name="",breaks=c(-1,0,1),limits = c(-1.5,1.5)) +
    scale_x_discrete(name="", labels=c("P mineral soil","P organic soil","P biomass")) +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    labs(title = "P ecosystem retention") +
    theme_bw() +
    theme(plot.margin = unit(c(-2,2,-2,2), "pt"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=0),
          axis.title=element_text(size=14,face="bold"))
  
  frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
  frame_sub<-frame_sub[frame_sub$variable%in%c("PLAB","PSEC","PMIN_c","PBMIN_c","PUP_c"),] 
  frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"))
  frame_sub<-droplevels(frame_sub)
  frame_sub$variable<-factor(frame_sub$variable,levels=c("PUP_c","PLAB","PSEC","PBMIN_c","PMIN_c"))
  frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")
  
  p4<-ggplot(data=frame_sub, aes(x=variable, y=value, fill=sign)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
    scale_y_continuous(name="",breaks=c(-1,0,1),limits = c(-1.7,1.7)) +
    scale_x_discrete(name="", labels=c("P uptake","P labile","P secondary","P bioch. min.","P net min.")) +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    labs(title = "P plant acquisition") +
    theme_bw() +
    theme(plot.margin = unit(c(-2,2,2,2), "pt"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=0),
          axis.title=element_text(size=14,face="bold"))
  
  ##panels with C pools 
  frame_sub<-frame[frame$MOD%in%c("CAB","GDP","ELM","ECA","ORC","POP")&frame$PRD=="final",]
  frame_sub<-frame_sub[frame_sub$variable%in%c("CBIO"),] 
  frame_sub$MOD<-factor(frame_sub$MOD,levels=c("CAB","GDP","ELM","POP","ECA","ORC"),labels=c("CABLE","GDAY","ELM-CTC","CABLE-POP","ELM-ECA","ORCHIDEE"))
  frame_sub<-droplevels(frame_sub)
  frame_sub$variable<-factor(frame_sub$variable,levels=c("CBIO"))
  frame_sub[["sign"]] = ifelse(frame_sub[["value"]] >= 0, "positive", "negative")
  
  #devide by 14.5 to get to annual values(see comments above)
  p5<-ggplot(data=frame_sub, aes(x=variable, y=value/14.5, fill=sign)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) + 
    coord_flip() + facet_grid(.~MOD) + #scales = "free_y"
    scale_y_continuous(name="", breaks=c(0,100)) + #, limits = c(-250, 250)) +
    scale_x_discrete(name="", labels=c("C biomass")) +
    scale_fill_manual(values = c("positive" = "steelblue4", "negative" = "firebrick4")) +
    labs(title = "C gain") +
    theme_bw() +
    theme(plot.margin = unit(c(2,2,-2,2), "pt"),
          strip.text = element_text(size=10, face="bold"),
          axis.text=element_text(size=12, face="bold"),
          axis.text.x=element_text(angle=0),
          axis.title=element_text(size=14,face="bold"))
  
  print(ggarrange(p5, p2, p3, p4, nrow=4, align="v", 
                  heights=c(3.8,6,5,8), labels=c("a","b","c","d")))
  setwd("~/SpiderOak Hive/AMAZONAS/MEETINGS/AGU2018_Washington")
  ggsave("FIG3_keymechanisms_NP_for_ppt.png",width=7,height=8)
  
  
  
     
  
 #plot initial conditions of C,N and P pools in one FIG against OBS in separate panel 
frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
ind<-frame$variable%in%c("CW","CSOIL","CROOT") #"CL",
frame<-frame[ind,];frame<-droplevels(frame);
frame$variable<-factor(frame$variable,levels=c("CW","CROOT","CSOIL")) #"CL",

levels(frame$variable)[which(levels(frame$variable)=="CW")]<-"C wood"
#levels(frame$variable)[which(levels(frame$variable)=="CL")]<-"C Leaves"
levels(frame$variable)[which(levels(frame$variable)=="CSOIL")]<-"C soil"
levels(frame$variable)[which(levels(frame$variable)=="CROOT")]<-"C root"

p1<-ggplot(data=frame, aes(x=MOD, y=value/1000, fill=variable)) +
    geom_bar(stat="identity",size=0.3,position=position_dodge(), width=0.8, colour="black", show.legend=T) +
    scale_y_continuous(name="C pools [kgC/m2]", limits=c(0,30)) + 
    scale_x_discrete(name="") +
    scale_fill_manual(values=brewer.pal(4,"Greys")[2:4],name="") +
    theme_bw() +
    theme(strip.text=element_text(size=18, face="bold"),
           axis.text=element_text(size=18, face="bold"),
           axis.title=element_text(size=18,face="bold"),
           legend.title=element_text(size=20,face="bold"),
           legend.text=element_text(size=18,face="bold"),
           legend.position=c(.9,.79))

frame2<-melt(obs)
ind<-frame2$variable%in%c("CW","CSOIL","CROOT") #"CL",
frame2<-frame2[ind,];frame2<-droplevels(frame2)
frame2$variable<-factor(frame2$variable,levels=c("CW","CROOT","CSOIL")) #"CL",

levels(frame2$variable)[which(levels(frame2$variable)=="CW")]<-"C wood"
#levels(frame2$variable)[which(levels(frame2$variable)=="CL")]<-"C leaves"
levels(frame2$variable)[which(levels(frame2$variable)=="CSOIL")]<-"C soil"
levels(frame2$variable)[which(levels(frame2$variable)=="CROOT")]<-"C root"

p1_obs<-ggplot(data=frame2,aes(x=MOD,y=value/1000,fill=variable)) +
        geom_bar(stat="identity",size=0.3, position=position_dodge(), width=0.8, colour="black", show.legend=F) +
        scale_y_continuous(name="", limits=c(0,30)) + 
        scale_x_discrete(name="") +
        scale_fill_manual(values=brewer.pal(4,"Greys")[2:4]) +
        theme_bw() +
        theme(strip.text=element_text(size=18, face="bold"),
               axis.text=element_text(size=18, face="bold"),
               axis.title=element_text(size=20,face="bold"))

widC<-c(6.55, 1)
# pdf("FIG1_initial_Cpools.pdf",width=10,height=4)
# print(ggarrange(p1, p1_obs, widths=widC, align="h", ncol=2)) #, nrow = 1) #labels = c("A", "B"),
# dev.off()

# N pools
frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
ind<-frame$variable%in%c("NPMIN","NPORG")&frame$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC","OCN","LPJ","GDA","JUL")
frame<-frame[ind,];frame<-droplevels(frame)
frame$variable<-factor(frame$variable,levels=c("NPORG","NPMIN"))

levels(frame$variable)[which(levels(frame$variable)=="NPMIN")]<-"N min."
levels(frame$variable)[which(levels(frame$variable)=="NPORG")]<-"N org."

p2<-ggplot(data=frame, aes(x=MOD, y=value, fill=variable)) +
    geom_bar(stat="identity", size=0.3, position=position_dodge(), width=0.8, colour="black") + 
    scale_y_continuous(name="N pools [gN/m2]",limits=c(0,2850)) + #, limits = c(-250, 250)) +
    scale_x_discrete(name="") +
    scale_fill_manual(values=brewer.pal(3,"Greys")[1:2],name="") +
    theme_bw() +
    theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"),
        legend.position=c(.82,.85))

frame2<-melt(obs)
ind<-frame2$variable%in%c("NPMIN","NPORG")
frame2<-frame2[ind,];frame2<-droplevels(frame2)
frame2$variable<-factor(frame2$variable,levels=c("NPORG","NPMIN"))
levels(frame2$variable)[which(levels(frame2$variable)=="NPORG")]<-"N org."
levels(frame2$variable)[which(levels(frame2$variable)=="NPMIN")]<-"N min."

p2_obs<-ggplot(data=frame2,aes(x=MOD,y=value,fill=variable)) +
  geom_bar(stat="identity", size=0.3, position=position_dodge(), width=0.8, colour="black", show.legend=FALSE) +
  scale_y_continuous(name="",limits=c(0,2850)) + 
  scale_x_discrete(name="") +
  scale_fill_manual(values=brewer.pal(3,"Greys")[c(1:2)],name="") +
  theme_bw() +
  theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"))

widN<-c(5.5, 1)
# pdf("FIG1_initial_Npools.pdf",width=10,height=4)
# print(ggarrange(p2, p2_obs, widths=widN, align="h", ncol=2)) #, nrow = 1) #labels = c("A", "B"),
# dev.off()

# P pools
frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
ind<-frame$variable%in%c("PSEC","PPORG","PLAB")&frame$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC") #"POCC"
frame<-frame[ind,];frame<-droplevels(frame)
frame$variable<-factor(frame$variable,levels=c("PSEC","PPORG","PLAB"))

levels(frame$variable)[which(levels(frame$variable)=="PLAB")]<-"P lab."
levels(frame$variable)[which(levels(frame$variable)=="PPORG")]<-"P org."
levels(frame$variable)[which(levels(frame$variable)=="PSEC")]<-"P sec."

p3<-ggplot(data=frame, aes(x=MOD, y=value, fill=variable)) +
  geom_bar(stat="identity", size=0.3, position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  scale_y_continuous(name="P pools [gP/m2]",breaks=seq(0,100,10),limits=c(0,100)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=brewer.pal(3,"Greys"),name="") +
  theme_bw() +
  theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"),
        legend.position=c(.8,.79))

frame3<-melt(obs)
ind<-frame3$variable%in%c("PSOIL_extract","PSOIL_extract_organic","PSOIL_resin") #"PSOIL_extract_inorganic","PSOIL_total",
frame3<-frame3[ind,];frame3<-droplevels(frame3)
frame3$variable<-factor(frame3$variable,levels=c("PSOIL_extract","PSOIL_extract_organic","PSOIL_resin"))

levels(frame3$variable)[which(levels(frame3$variable)=="PSOIL_extract")]<-"P non-occ."
levels(frame3$variable)[which(levels(frame3$variable)=="PSOIL_extract_organic")]<-"P org."
levels(frame3$variable)[which(levels(frame3$variable)=="PSOIL_resin")]<-"P resin"

p3_obs<-ggplot(data=frame3,aes(x=MOD,y=value,fill=variable)) +
  geom_bar(stat="identity", size=0.3,position=position_dodge(), width=0.8, colour="black") +
  scale_y_continuous(name="",breaks = c(seq(0,10,10)),limits=c(0,90)) + 
  scale_x_discrete(name="") +
  scale_fill_manual(values=brewer.pal(3,"Greys"),name="") +
  theme_bw() +
  theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"))
  #      legend.position=c(.8,.79))

widP<-c(10.7, 6)
# pdf("FIG1_initial_Ppools.pdf",width=10,height=4)
# print(ggarrange(p3, p3_obs, widths=widP, align="h", ncol=2)) #, nrow = 1) #labels = c("A", "B"),
# dev.off()

# initial C, N and P pools
pdf("FIGX_initial_CNPpools.pdf",width=10,height=12)
print(ggarrange(ggarrange(p1, p1_obs, widths=widC, align="h", ncol=2),
ggarrange(p2, p2_obs, widths=widN, align="h", ncol=2),
ggarrange(p3, p3_obs, widths=widP, align="h", ncol=2), nrow=3))
dev.off()



#plot for Rich
  
  ofile <- paste("CO2resp_GPP_NPP_NEP",maxyear,"bio.pdf",sep='_')
  pdf(ofile,width=10,height=8)

    plotgpp<-xyplot(GPP/1000~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=expression(paste("GPP [kg C ",m^-2," ",yr^-1,"]"),sep=""),xlab=NULL,
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main='', 
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(v=xtick,col="grey85")})
    plotnpp<-xyplot(NPP/1000~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=expression(paste("NPP [kg C ",m^-2," ",yr^-1,"]"),sep=""),xlab=NULL,
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main='',
                    panel=function(...){
                    panel.xyplot(...)
                    panel.abline(v=xtick,col="grey85")},
                    key=list(space='top',border=T,columns=5,text=list(levels(run1$MOD)),
                             lines=list(col=cols$COL,lty=1,lwd=2)))
    plotnep<-xyplot(NEP/1000~YEARbio,data=run1,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab=expression(paste("NEP [kg C ",m^-2," ",yr^-1,"]"),sep=""),
                    scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main='',
                    panel=function(...){
                      panel.xyplot(...)
                      panel.abline(h=0,col='grey40')
                      panel.abline(v=xtick,col="grey85")})
    grid.arrange(plotgpp,plotnpp,plotnep,ncol=1)
    
  dev.off()
  

ofile <- paste("fr7rel_init_final15_bio.pdf",sep='')

frame<-melt(fr7_relgr,id=c("MOD","PRD"))
levels(frame$PRD) <- c("1st YEAR", "15YEARS", "CHANGE")

pdf(ofile,width=10,height=8)

#combine variables in one plot
ind<-frame$variable%in%c("GPP","NPP","NUP","PUP")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CBIO","GPP","NPP","CUE")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

#combine variables in one plot
ind<-frame$variable%in%c("CBIO","PUP","PUE","LeafNP","WoodNP")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))


#combine variables in one plot
ind<-frame$variable%in%c("CBIO","NUP","NUE","LeafCN","WoodNP")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

#combine variables in one plot
ind<-frame$variable%in%c("GPP","NPP","CSOIL","CBIO")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

#combine variables in one plot
ind<-frame$variable%in%c("GPP","NPP","CUE","GPP_CBIO","RAU_CBIO")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CALLOC_L","CALLOC_W","CALLOC_FR") #,"CALLOC_CR")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("NUP","NUE","NPP_N","NCON","LeafCN","WoodCN","SoilCN")
frame$value[ind==T&(frame$MOD=="ECA")]<-0  #ECA is mucking up the graph (leaching)

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("PUP","PUE","NPP_P","LeafNP","SoilNP") #,"PLAB"
#frame$value[ind==T&(frame$MOD=="ORC")]<-0  #ORC is mucking up the graph

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("NLEACH","PLEACH","NECO","PECO")  #"PPAR" is nothing"RO","DRAIN",

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("ET","EC","ES","T","SW","GCd","WUE")  #"PPAR" is nothing"RO","DRAIN",

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

dev.off()

#now plot only relative effect at the end of the experiment:

frame2<-fr7[fr7$YEARbio==2013,c("YEARbio","MOD","CBIO")]
ind<-frame2$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC")
frame2<-frame2[ind,]
frame2$MOD<-factor(frame2$MOD,c("POP","CAB","GDP","ELM","ECA","ORC"))

cols2<-cols[c(2,1,5,3,4,6),]

ofile <- paste("Pmodels_CBIO.pdf",sep='')
pdf(ofile,width=12,height=6)

ggplot(data=frame2, aes(x=MOD, y=CBIO/1000, fill=MOD)) +
  geom_bar(stat="identity", width=0.8, colour="black") + #, show.legend=FALSE) + position=position_dodge()
 # facet_grid(PRD~.) +
  scale_y_continuous(name="Biomass C in kg/m2 [ELE-AMB]") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols2$COL,name="MODEL") +
  theme(strip.text=element_text(size=22, face="bold"),
        axis.text=element_text(size=22, face="bold"),
        axis.title=element_text(size=24,face="bold"),
        legend.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20,face="bold"))

dev.off()


ofile <- paste("fr7rel_final15_bio.pdf",sep='')
pdf(ofile,width=10,height=4.5)

ind<-frame$variable%in%c("GPP","NPP","CBIO","PUP","PUE")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("NPP","LAI","LMA","CL")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CALLOC_L","CALLOC_W","CALLOC_FR") #,"CALLOC_CR")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CL_RES","CW_RES","CFR_RES","CCR_RES")
ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CSOIL_RES","RHET_CSOIL")
ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

dev.off()

#changes in C, N, P pools (only long-term and absolute effect)
frame<-melt(fr7_gr,id=c("MOD","PRD"))
frame<-frame[frame$PRD=="final",]
frame<-droplevels(frame)
levels(frame$PRD) <- c("15YEARS")

ofile <- paste("fr7abs_final15_bio.pdf",sep='')

pdf(ofile,width=10,height=8)

ind<-frame$variable%in%c("CECO","CBIO","CSOIL","CROOT")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in gC/m2") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("CW","CL","CROOT","CSOIL")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in gC/m2") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("NPMIN","NBIO","NPORG")  #"PPAR" is nothing

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in gN/m2") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("PECO","PBIO","PPORG","PSEC","PLAB")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in gP/m2") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

ind<-frame$variable%in%c("PECO","PBIO","PPORG","PPMIN","PSEC","POCC","PLAB")

ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response in gP/m2") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"))

dev.off()



#ambient versus CO2 response of leaf plasticity
ind<-dfa$SIM=="AMB_OBS"
ann<-aggregate(dfa[ind,c("LeafCN","LeafNP")],by=list(dfa$MOD[ind]),FUN=mean)
names(ann)[1]<-"MOD"
frame<-melt(ann)

frame$value[frame$variable=="LeafNP"&frame$MOD=="GDA"]<-NA
frame$value[frame$variable=="LeafCN"&frame$MOD=="ED2"]<-NA
frame_obs<-data.frame(variable = levels(frame$variable), Z = c(23,37))

p1<-ggplot(data=frame, aes(x=MOD, y=value, fill=MOD)) +
  geom_bar(stat="identity",position="identity",width=0.8,colour="black",show.legend=F) +
  facet_wrap(~variable, scales = "free_y", nrow=2) +
  geom_hline(data=frame_obs,aes(yintercept = Z),linetype=2) +
  scale_y_continuous(name="Ambient conditions") +
  scale_x_discrete(name="",labels=NULL,breaks=NULL) +
  scale_fill_manual(values=cols$COL) +
  theme_bw() +
  theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"))

frame<-melt(fr7_gr,id=c("MOD","PRD"))
frame<-frame[frame$PRD=="final",]
frame<-droplevels(frame)
levels(frame$PRD) <- c("15YEARS")

ind<-frame$variable%in%c("LeafCN","LeafNP")

p2<-ggplot(data=frame[ind,], aes(x=MOD, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black", show.legend=T) + 
  facet_wrap(~variable,scale="free_y",nrow=2) +
  scale_y_continuous(name="CO2 response") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="",labels=NULL,breaks=NULL) +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme_bw() +
  theme(strip.text=element_text(size=18, face="bold"),
        axis.text=element_text(size=18, face="bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18,face="bold"))

pdf("FIG_leafplasticity_ambient_CO2.pdf",width=8,height=6)
print(ggarrange(ggarrange(p1, p2, widths=c(1,1), align="h", ncol=2)))
dev.off()










frame<-melt(fr7_relgr,id=c("MOD","PRD"))
levels(frame$PRD) <- c("1st YEAR", "10YEARS", "CHANGE")

ind<-frame$PRD=="10YEARS"
frame<-frame[ind==T,]
frame<-droplevels(frame)

levels(frame$variable)[which(levels(frame$variable)=="GPP")]<-"Photosynthesis"
levels(frame$variable)[which(levels(frame$variable)=="GCd")]<-"Conductance"
levels(frame$variable)[which(levels(frame$variable)=="CECO")]<-"Ecosystem Carbon"
levels(frame$variable)[which(levels(frame$variable)=="CALLOC_BG")]<-"Root allocation"
levels(frame$variable)[which(levels(frame$variable)=="RHET_CSOIL")]<-"Soil Carbon respiration"
levels(frame$variable)[which(levels(frame$variable)=="CBIO")]<-"Biomass Carbon"


pdf("bars_keyvars.pdf",width=12,height=8)

ind<-frame$variable%in%c("Photosynthesis","Root allocation","Biomass Carbon","Conductance")
ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  # facet_grid(PRD~.) +
  scale_y_continuous(name="CO2 response after 10 years in %") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  theme(strip.text=element_text(size=14, face="bold"),
        axis.text=element_text(size=14, face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12,face="bold"))

dev.off()

pdf("bars_keyvars2.pdf",width=16,height=6)

ind<-frame$variable%in%c("Photosynthesis","Root allocation","Biomass Carbon","Soil Carbon respiration")
ggplot(data=frame[ind,], aes(x=variable, y=value, fill=MOD)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
  # facet_grid(PRD~.) +
  scale_y_continuous(name="") + #, limits = c(-250, 250)) +
  scale_x_discrete(name="") +
  scale_fill_manual(values=cols$COL,name="MODEL") +
  ggtitle("CO2 response in % after 10 years") + 
  theme(strip.text=element_text(size=22, face="bold"),
        plot.title = element_text(size=28, face="bold"),
        axis.text=element_text(size=22, face="bold"),
        axis.title=element_text(size=22,face="bold"),
        legend.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=18,face="bold"))

dev.off()       





































# # CO2 effect on CECO and CBIOMASS (1)
# #----------------------------------------------------------------------------------------------------
# pdf(paste("CBIO_CECO_abs_",maxyear,".pdf",sep=""),width=10,height=8)
# 
# p1<-xyplot(CBIO/1000~YEARbio,data=fr7,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab="",xlab=NULL,
#            scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('CO2 induced biomass C'),
#            panel=function(...){
#              panel.xyplot(...)
#              panel.abline(v=xtick,col="grey85")})
# p2<-xyplot(CECO/1000~YEARbio,data=fr7,groups=MOD,lty=1,pch=19,col=cols$COL,type=c("b","g"),ylab="",
#            scales=list(tck=c(-0.5,0),alternating=F,relation='free'),main=paste('CO2 induced ecosystem C'),
#            panel=function(...){
#              panel.xyplot(...)
#              panel.abline(h=0,col='grey40')
#              panel.abline(v=xtick,col="grey85")},
#            key=list(space='top',border=T,columns=5,text=list(levels(fr7_rel$MOD)),
#                     lines=list(col=cols$COL,lty=1,lwd=2)))
# print(p1,split=c(1,1,1,2),more=T)
# print(p2,split=c(1,2,1,2),more=F)
# dev.off()


# # assess T and CO2 interactions: HOW
# frame<-data.frame(MOD=fr7$MOD,YEAR=fr7$YEARbio,TAIR=dfa$TAIR[dfa$SIM=="AMB_OBS"],GPP=fr7$GPP)
# 
# ggplot(data=frame, aes(x=TAIR, y=GPP, col=MOD)) +
#   geom_line()
# 
# 
# # assess T and CO2 interactions: daily GPP vs TAIR: gams for ambient and elevated CO2
# setwd('~/SpiderOak Hive/AMAZONAS/FACE-MEI/model_output/all_models')
# 
# dfa<-readRDS("amaface_allruns_daily_R2.rds")
# units<- readRDS('units_frame.rds')
# date <- Sys.Date()
# 
# modcod<-c("CAB","POP","LPJ","OCN","GDA","GDP","ECA","ELM","JUL","ORC","INL","ED2")
# dfa$MOD<-factor(dfa$MOD, levels = c("CAB","POP","ELM","ECA","GDP","ORC","LPJ","OCN","GDA","JUL","ED2","INL"))
# #define colours to be sued throughout, use paired which has 12 colours, keep each model fixed with one colour like this
# cols=data.frame(MOD=c(levels(dfa$MOD)),COL=brewer.pal(12,"Paired")[1:length(modcod)],stringsAsFactors = F)
# 


# # plot C allocation over time, only long term effect of % change
# 
# frame<-melt(fr7_rel,id=c("MOD","YEARbio"))
# frame$MODgr[frame$MOD%in%c("CAB","JUL","INL","ED2","FAT")]<-"fixed allocation"
# frame$MODgr[frame$MOD%in%c("ECA","GDP","GDA","POP")]<-"dynamic allocation1"
# frame$MODgr[frame$MOD%in%c("ORC","LPJ","OCN","ELM")]<-"dynamic allocation2"
# frame$MODgr<-as.factor(frame$MODgr)
# 
# ofile <- paste("fr7rel_Calloc.pdf",sep='')
# pdf(ofile,width=10,height=8)
# 
# for (i in levels(frame$MODgr)) {
#   
# ind<-frame$variable%in%c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")&frame$MODgr==i
# p1<-ggplot(data=frame[ind,], aes(x=YEARbio, y=value, group=variable)) +
#   geom_line(aes(col=variable)) + #, show.legend=FALSE) + 
#   facet_grid(MOD~.) +
#   scale_y_continuous(name="CO2 response in %") + #, limits = c(-250, 250)) +
#   scale_x_continuous(name="") +
#   ggtitle(i) +
#   scale_fill_manual(values=c("red","blue"),name="allocation") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# print(p1)
# 
# }
# dev.off()
# 
# 

# 
# 
# 
# 
# 
# # CNP stocks and stoiciometry evaluation plots, should sit somewhere else than here
# #------------------------------------------------------------------------------------------------------------
# 
# ofile <- "evaluation_CNPstocks_2000.pdf"
# 
# pdf(ofile,width=10,height=8)
# 
# ind<-dfa$YEARbio==2000&dfa$SIM=="AMB_OBS"
# #somehow add the observations     ,c("OBS",1999,20,0,0,0,0,5.3))
# 
# frame<-melt(dfa[ind,c("MOD","YEARbio","CBIO","CFR","CCR","CL","CSOIL","LAI","LMA")],id=c("MOD","YEARbio"))
# #frame$value[frame$MOD=="JUL"]<-0  #JUL is mucking up the graph
# 
# frame$value[frame$variable=="CBIO"]<-frame$value[frame$variable=="CBIO"]/1000
# frame$value[frame$variable=="CCR"]<-frame$value[frame$variable=="CCR"]/1000
# frame$value[frame$variable=="CSOIL"]<-frame$value[frame$variable=="CSOIL"]/1000
# frame$value[frame$variable=="CFR"]<-frame$value[frame$variable=="CFR"]/1000
# frame$value[frame$variable=="CL"]<-frame$value[frame$variable=="CL"]/1000
# 
# 
# ggplot(data=frame, aes(x=YEARbio, y=value, fill=MOD)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
#   facet_wrap(~variable, scales = "free_y", nrow=3) +
#   scale_y_continuous(name="") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c(cols$COL,"black")) +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# 
# dev.off()
# 
# 
# #plot initial conditios
# 
# pdf("initial_Cpools.pdf")
# frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
# 
# ind<-frame$variable%in%c("CW","CL","CSOIL","CROOT")
# par(mfrow=(c(2,1)))
# 
# ggplot(data=frame[ind,], aes(x=MOD, y=value/1000, fill=variable)) +
#   geom_bar(stat="identity", position=position_stack(), width=0.8, colour="black") + #, show.legend=FALSE) + 
#   #  facet_grid(PRD~.) +
#   scale_y_continuous(name="Initial C pools in kgC/m2", limits=c(0,50)) + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c("darkgreen","brown","orange","yellowgreen"),name="MODEL") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# dev.off()
# 
# 
# #plot initial P pools 
# pdf("initial_Ppools_b.pdf")
# frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
# 
# ind<-frame$variable%in%c("PLAB","PPORG")#&frame$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC")
# par(mfrow=(c(2,1)))
# 
# ggplot(data=frame[ind,], aes(x=MOD, y=value, fill=variable)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
#   #  facet_grid(PRD~.) +
#   scale_y_continuous(name="Initial P pools in gP/m2") +#, limits=c(0,100)) + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c("turquoise","orange","brown","chocolate","darkgreen"),name="MODEL") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# dev.off()
# 
# 
# pdf("initial_Npools.pdf")
# frame<-melt(dfa[dfa$SIM%in%c("AMB_OBS")&dfa$YEARbio==1999,])
# 
# ind<-frame$variable%in%c("NPMIN","NBIO","NPORG")&frame$MOD%in%c("CAB","POP","ELM","ECA","GDP","ORC","OCN","LPJ","GDA","JUL","ED2")
# par(mfrow=(c(2,1)))
# 
# ggplot(data=frame[ind,], aes(x=MOD, y=value, fill=variable)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8, colour="black") + #, show.legend=FALSE) + 
#   #  facet_grid(PRD~.) +
#   scale_y_continuous(name="Initial N pools in gN/m2") +#, limits=c(0,100)) + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c("turquoise","orange","darkgreen"),name="MODEL") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# dev.off()
# 

# #FIG2: alternative, only barchart
# #now plot only relative effect at the end of the experiment:
# frame2<-fr7[fr7$YEARbio==2013,c("YEARbio","MOD","CBIO")]
# 
# p1<-ggplot(data=frame2, aes(x=MOD, y=CBIO/1000, fill=MOD)) +
#   geom_bar(stat="identity", width=0.8, colour="black") + #, show.legend=FALSE) + position=position_dodge()
#   scale_y_continuous(name="Biomass C in kg/m2 [ELE-AMB]") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=cols$COL,name="MODEL") +
#   theme(strip.text=element_text(size=22, face="bold"),
#         axis.text=element_text(size=22, face="bold"),
#         axis.title=element_text(size=24,face="bold"),
#         legend.title=element_text(size=22,face="bold"),
#         legend.text=element_text(size=20,face="bold"))
# 
# pdf(paste("FIG2B_CO2_CBIO_abs_bars.pdf",sep=""),width=10,height=5)
# print(p1)
# dev.off()


} #end if plotallfig=T


# 
# #xyplot of biomass C response
# xtick<-seq(1999,maxyear,1)
# #have to move x-title and legend away from each other still
# plotabs<-xyplot(CBIO/1000~YEARbio-1998,data=fr7[fr7$YEARbio<=maxyear,],groups=MOD,lty=1,lwd=2.5,pch=19,col=cols$COL,type=c("b","g"),
#                 ylab=list("Biomass C change [kgC/m2]",cex=1.5),xlab=list('',cex=1.5),
#                 scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#                 panel=function(...){
#                   panel.xyplot(...)
#                   #    panel.abline(h=0,col='grey40')
#                   panel.abline(v=xtick-1998,col="grey85")
#                   panel.abline(h=seq(0,3,0.5),col="grey85")},
#                 key=list(space='bottom',border=T,columns=5,text=list(levels(fr7$MOD)),
#                          lines=list(col=cols$COL,lty=1,lwd=2.5)))
# 
# plotrel<-xyplot(CBIO~YEARbio-1998,data=fr7_rel[fr7_rel$YEARbio<=maxyear,],groups=MOD,lty=1,lwd=2.5,pch=19,col=cols$COL,type=c("b","g"),
#                 ylab=list("Biomass C change [%]",cex=1.5),xlab=list('Years of CO2 fumigation',cex=1.5),
#                 scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#                 panel=function(...){
#                   panel.xyplot(...)
#                   #    panel.abline(h=0,col='grey40')
#                   panel.abline(v=xtick-1998,col="grey85")
#                   panel.abline(h=seq(0,3,0.5),col="grey85")})
# #   key=list(space='bottom',border=T,columns=5,text=list(levels(fr7_rel$MOD)),
# #          lines=list(col=cols$COL,lty=1,lwd=2.5)))

# print(plotabs,split=c(1,1,1,2),more=T)
# print(plotrel,split=c(1,2,1,2),more=F)
# dev.off()
# 
# 
# 
# 
# 
# #over time per group
# pdf("Diagnostics_Fig1B_GPP_LAI_annual_Group.pdf",width=6,height=8)
# frame<-cbind(aggregate(fr7_rel[,c("GPP","LAI","GPP_LAI")], by=list(fr7_rel$YEARbio,fr7_rel$MODgr),FUN=mean))
# names(frame)[1:2]<-c("YEARbio","MODgr")
# xtick<-seq(1999,maxyear,1)
# plot1<-xyplot(GPP~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr,lty=1,lwd=2.5,pch=19,col=c("Grey40","Steelblue4","Springgreen4"),type=c("b","g"),
#               ylab=list("GPP % ",cex=1.5),xlab=list('',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")})
# 
# plot2<-xyplot(LAI~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr,lty=1,lwd=2.5,pch=19,col=c("Grey40","Steelblue4","Springgreen4"),type=c("b","g"),
#               ylab=list("LAI % ",cex=1.5),xlab=list('',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")})
# 
# 
# plot3<-xyplot(GPP_LAI~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr,lty=1,lwd=2.5,pch=19,col=c("Grey40","Steelblue4","Springgreen4"),type=c("b","g"),
#               ylab=list("GPP/LAI % ",cex=1.5),xlab=list('Years of CO2 fumigation',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 #    panel.abline(h=0,col='grey40')
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")},
#               key=list(space='top',border=T,columns=3,text=list(levels(frame$MODgr)),
#                        lines=list(col=c("Grey40","Steelblue4","Springgreen4"),lty=1,lwd=3)))
# 
# grid.arrange(plot1,plot2,plot3, nrow=3)
# dev.off()


# # for DGVMs - leave that out for now
# frame<-cbind(aggregate(fr7_rel[,c("GPP","LAI","GPP_LAI")], by=list(fr7_rel$MODgr2),FUN=mean))
# frame<-melt(frame,id.vars = c("Group.1"))
# sec<-aggregate(fr7_rel[,c("GPP","LAI","GPP_LAI")], by=list(fr7_rel$MODgr2),FUN=sd)
# frame<-cbind(frame,melt(sec,id.vars=c("Group.1"))[,3]); rm(sec)
# colnames(frame)<-c("MODgr2","variable","mean","sd")
# 
# pdf("Suppl_FIG2_GPP_LAI_Group2.pdf",width=4,height=4)
# 
# ggplot(data=frame, aes(x=MODgr2, y=mean, fill=variable)) +
#   #facet_grid(.~MODgr) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8) + #, show.legend=FALSE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   scale_y_continuous(name="CO2 response in % (mean + sd)") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="MODEL") +
#   scale_fill_manual(values=c("Steelblue4","Springgreen4","Grey40"),name="Variables") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# dev.off()

# pdf("Suppl_FIG2_GPP_LAI_annual_Group2.pdf",width=6,height=8)
# frame<-cbind(aggregate(fr7_rel[,c("GPP","LAI","GPP_LAI")], by=list(fr7_rel$YEARbio,fr7_rel$MODgr2),FUN=mean))
# names(frame)[1:2]<-c("YEARbio","MODgr2")
# xtick<-seq(1999,maxyear,1)
# plot1<-xyplot(GPP~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr2,lty=1,lwd=2.5,pch=19,col=c("Springgreen4","Steelblue4"),type=c("b","g"),
#               ylab=list("GPP % ",cex=1.5),xlab=list('',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")})
# 
# plot2<-xyplot(LAI~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr2,lty=1,lwd=2.5,pch=19,col=c("Springgreen4","Steelblue4"),type=c("b","g"),
#               ylab=list("LAI % ",cex=1.5),xlab=list('',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")})
# 
# 
# plot3<-xyplot(GPP_LAI~YEARbio-1998,data=frame[frame$YEARbio<=maxyear,],groups=MODgr2,lty=1,lwd=2.5,pch=19,col=c("Springgreen4","Steelblue4"),type=c("b","g"),
#               ylab=list("GPP/LAI % ",cex=1.5),xlab=list('Years of CO2 fumigation',cex=1.5),
#               scales=list(tck=c(-0.5,0),alternating=F,relation='free',x=list(cex=1.5),y=list(cex=1.5)),main=paste(''),
#               panel=function(...){
#                 panel.xyplot(...)
#                 #    panel.abline(h=0,col='grey40')
#                 panel.abline(v=xtick-1998,col="grey85")
#                 panel.abline(h=seq(0,30,5),col="grey85")},
#               key=list(space='top',border=T,columns=2,text=list(levels(frame$MODgr2)),
#                        lines=list(col=c("Springgreen4","Steelblue4"),lty=1,lwd=3)))
# grid.arrange(plot1,plot2,plot3, nrow=3)
# dev.off()

# calculate ambient, after 3 and 15 years (leave that out, not very visually informative)
# #calculate mean and sd for each model for one run or frame, for partiuclar set of variables
# ind<-dfa$YEARbio<=maxyear
# frame2<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=mean))
# frame2<-melt(frame2,id.vars = c("Group.1","Group.2"))
# sec2<-aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=sd)
# frame2<-cbind(frame2,melt(sec2,id.vars=c("Group.1","Group.2"))[,4]); rm(sec2)
# colnames(frame2)<-c("MOD","SIM","variable","mean","sd")
# levels(frame2$SIM)<-c("AMB_OBS","ELE_OBS")
# 
# ind<-dfa$YEARbio<=maxyear
# frame<-cbind(aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=mean))
# frame<-melt(frame,id.vars = c("Group.1","Group.2"))
# sec<-aggregate(dfa[ind,c("CALLOC_L","CALLOC_W","CALLOC_FR","CALLOC_CR")], by=list(dfa$MOD[ind],dfa$SIM[ind]),FUN=sd)
# frame<-cbind(frame,melt(sec,id.vars=c("Group.1","Group.2"))[,4]); rm(sec)
# colnames(frame)<-c("MOD","SIM","variable","mean","sd")
# levels(frame$SIM)<-c("AMB_OBS","ELE_OBS_15")
# 
# frame2<-rbind(frame2,frame[frame$SIM=="ELE_OBS_15",])
# 
# p1<-ggplot(data=frame2, aes(x=MOD, y=mean*100, fill=SIM)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8) + #, show.legend=FALSE) +
#   geom_errorbar(aes(ymin=mean*100-sd*100, ymax=mean*100+sd*100),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   facet_grid(variable~.) +
#   scale_y_continuous(name="C allocation in % (mean + sd)") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="MODEL") +
#   scale_fill_manual(values=c("firebrick4","steelblue2","steelblue4"),name="Simulation") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# pdf("Suppl_FIG2_CO2_Callocation_bars.pdf",width=10,height=6)
# print(p1)
# dev.off()

#fraction per tissue, this is another way of looking at it, how much C is really there,
#and also are the alloc variables correct....
#do not seem to be, e.g. ELM

# #this is not right yet...... fr7 on a fraction is confusing......!!!!
# frame<-fr7[fr7$YEARbio<=maxyear,c("MOD","YEARbio","CL","CW","CFR","CCR","CBIO")]
# frame$CL_CBIO<-frame$CL/frame$CBIO
# frame$CW_CBIO<-frame$CW/frame$CBIO
# frame$CFR_CBIO<-frame$CFR/frame$CBIO
# frame$CCR_CBIO<-frame$CCR/frame$CBIO
# 
# #frame<-frame[!frame$MOD%in%c("INL"),]
# frame<-frame[,-7]
# frame<-melt(frame,id.vars = c("MOD","YEARbio"))
# ind<-frame$variable%in%c("CL_CBIO","CW_CBIO","CFR_CBIO","CCR_CBIO")
# 
# p1<-ggplot(data=frame[ind,], aes(x=YEARbio-1998, y=value, col=variable)) +
#   geom_line() +
#   facet_wrap(~MOD,scales = "free_y") + #free_y
#   scale_y_continuous(name="C allocation in % ", limits = c(-0.15, 1.05)) +
#   scale_x_continuous(name="Years of eCO2") +
#   scale_colour_manual(values=c("Springgreen4","firebrick4","steelblue2","steelblue4"),name="Allocation") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# pdf("Suppl_Tissue_CBIO_CO2effectannual.pdf",width=10,height=6)
# print(p1)
# dev.off()

# ind<-dfa$YEARbio<=maxyear&dfa$SIM=="AMB_OBS"
# frame2<-cbind(aggregate(dfa[ind,c("CL","CW","CFR","CCR","CBIO")], by=list(dfa$MOD[ind]),FUN=mean))
# frame2$CL_CBIO<-frame2$CL/frame2$CBIO
# frame2$CW_CBIO<-frame2$CW/frame2$CBIO
# frame2$CFR_CBIO<-frame2$CFR/frame2$CBIO
# frame2$CCR_CBIO<-frame2$CCR/frame2$CBIO
# frame2<-melt(frame2)
# 
# # frame2<-melt(frame2,id.vars = c("Group.1"))
# # sec2<-aggregate(dfa[ind,c("CL","CW","CFR","CCR")], by=list(dfa$MOD[ind]),FUN=sd)
# # frame2<-cbind(frame2,melt(sec2,id.vars=c("Group.1"))[,3]); rm(sec2)
# colnames(frame2)<-c("MOD","variable","value") #,"sd")
# #levels(frame2$variable)<-c("Leaf","Wood","Fine root","Coarse root")
# frame2<-frame2[!frame2$MOD%in%c("INL"),]
# 
# ind<-frame2$variable%in%c("CL_CBIO","CW_CBIO","CFR_CBIO","CCR_CBIO")
# 
# p2<-ggplot(data=frame2[ind,], aes(x=variable, y=value*100, fill=variable)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8) + #, show.legend=FALSE) +
#   # geom_errorbar(aes(ymin=mean*100-sd*100, ymax=value*100+sd*100),width=.2,                    # Width of the error bars
#   #               position=position_dodge(.9)) +
#   facet_wrap(~MOD) +
#   scale_y_continuous(name="Ambient C per tissue in % (mean + sd)") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="",labels=NULL) +
#   scale_fill_manual(values=c("Springgreen4","firebrick4","steelblue2","steelblue4"),name="Allocation") +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=12, face="bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title=element_text(size=12,face="bold"),
#         legend.text=element_text(size=10,face="bold"))
# 
# pdf("Suppl_FIG2_ambient_Callocation.pdf",width=10,height=6)
# print(p2)
# dev.off()

# #try to combine the NUP/PUP story with the CN model response in C:N and NUE (since not shown in paper paper)
# #want to show absolute values in of NUE/C:N here amb/ele after 15 years
# 
# dfa$CFR_CBIO<-(dfa$CFR/dfa$CBIO)*100
# vars<-c("NUE","LeafCN","WoodCN","PUE","LeafNP","WoodNP","CFR_CBIO","CBIO","CFR")
# ind<-dfa$YEARbio>2010&dfa$YEARbio<2016  #same frame as in fr7_gr/relgr... 
# frame2<-aggregate(dfa[ind,vars], by=list(dfa$MODgr[ind],dfa$SIM[ind]),FUN=mean,na.rm=T)
# names(frame2)[1:2]<-c("Group","SIM")
# frame2<-melt(frame2)
# sec<-aggregate(dfa[ind,vars], by=list(dfa$MODgr[ind],dfa$SIM[ind]),FUN=sd,na.rm=T)
# frame2<-cbind(frame2,melt(sec,id.vars=c("Group.1","Group.2"))[,-c(1:3)]); rm(sec)
# colnames(frame2)<-c("MODgr","SIM","variable","mean","sd")
# levels(frame2$MODgr)<-c("C-only models", "CN models", "CNP models")
# levels(frame2$SIM)<-c("Ambient","eCO2")
# 
# 
# ind<-frame2$MODgr%in%c("CN models")
# ind2<-frame2$variable%in%c("NUE","LeafCN","WoodCN")
# 
# p5<-ggplot(data=frame2[ind==T&ind2==T,],aes(x=SIM,y=mean,position="dodge",fill=SIM))+
#   geom_bar(stat ="identity") +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   facet_wrap(~variable,scales="free") +
#   scale_y_continuous(name="") +
#   scale_x_discrete(name="",labels=NULL) +
#   scale_fill_manual(values=c("Grey80","Black"),name="Simulation") +
#   labs(title = "Responses to CO2 \nCN models") +
#   theme_bw() +
#   theme(strip.text=element_text(size=16, face="bold"),
#         axis.text=element_text(size=16, face="bold"),
#         axis.text.x=element_text(angle=45, vjust=0.5),
#         axis.title=element_text(size=16,face="bold"),
#         plot.title = element_text(size=16,face="bold"))
# 
# pdf("Diagnostics_FIG3_actual_NUE_CN_CNmodels.pdf",width=6,height=4)
# print(p5) 
# dev.off()

# combining all P models is not that informative since only few really responde
# for CN models it makes sense, since they behave similarly

# ind<-frame2$MODgr%in%c("CNP models")
# ind2<-frame2$variable%in%c("PUE","LeafNP","WoodNP")
# 
# p6<-ggplot(data=frame2[ind==T&ind2==T,],aes(x=SIM,y=mean,position="dodge",fill=SIM))+
#   geom_bar(stat ="identity") +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   facet_wrap(~variable,scales="free") +
#   scale_y_continuous(name="") +
#   scale_x_discrete(name="",labels=NULL) +
#   scale_fill_manual(values=c("Grey80","Black")) +
#   labs(title = "Responses to CO2 \nCN models") +
#   theme_bw() +
#   theme(strip.text=element_text(size=16, face="bold"),
#         axis.text=element_text(size=16, face="bold"),
#         axis.text.x=element_text(angle=45, vjust=0.5),
#         axis.title=element_text(size=16,face="bold"),
#         plot.title = element_text(size=16,face="bold"))

#it doesnt work combined...
# pdf("FIG3b_combined_NUP_PUP_GPP_NUE_PUE_modelgroup.pdf",width=8,height=6)
# grid.arrange(arrangeGrob(p1, ncol=1),
#              arrangeGrob(p5, ncol=1),
#              arrangeGrob(p2, ncol=1),
#              arrangeGrob(p3, ncol=1), heights=c(1,1,1,1)) #,widths=c(1,1)) 
# dev.off()


# #FIG3: same thing but now all the points per group against GPP
# frame<-aggregate(dfa[,c("NPP","GPP","NUP","NUE","PUP","PUE")], by=list(dfa$MODgr,dfa$SIM,dfa$YEARbio),FUN=mean)
# names(frame)[1:3]<-c("Group","SIM","Year")

# CUE leave out: 
# ggplot(data=frame, aes(x=Year, y=NPP/GPP, shape=SIM)) +
#   geom_line(aes(linetype=SIM))+
#   geom_point(show.legend=FALSE) + 
#   facet_grid(Group~.) +
#   scale_y_continuous(name="CUE = NPP/GPP") + #, limits = c(-250, 250)) +
#   scale_x_discrete(name="Year") +
#   scale_shape_manual(values=c(19,1),name="Simulation") +
#   theme(strip.text=element_text(size=18, face="bold"),
#         axis.text=element_text(size=18, face="bold"),
#         axis.title=element_text(size=20,face="bold"),
#         legend.title=element_text(size=20,face="bold"),
#         legend.text=element_text(size=18,face="bold"))

# #now show GPP versus NUP/GPP, the same as in the line graphs but against GPP
# #(not NPP like YP, which makes less sense since NPP is downregulated in CAB e.g. and not GPP)
# ind<-frame$Group%in%c("CN")
# levels(frame$SIM)<-c("Ambient","eCO2")
# p1<-ggplot(data=frame[ind,], aes(x=NUP/GPP, y=GPP/1000, shape=SIM)) +
#   geom_point(show.legend=F,size=3) + 
#   labs(title = "CN models: Annual NUP/GPP vs. GPP") +
#   scale_y_continuous(name="",limits = c(2.3, 3.5)) +
#   scale_x_continuous(name="") +
#   scale_shape_manual(values=c(19,1),name="Simulation") +
#   theme_bw() +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=10, face="bold"),
#         plot.title = element_text(size=14,face="bold"))
# 
# ind<-frame$Group%in%c("CNP")
# p2<-ggplot(data=frame[ind,], aes(x=NUP/GPP, y=GPP/1000, shape=SIM)) +
#   geom_point(show.legend=F,size=3) + 
#   labs(title = "CNP models: Annual NUP/GPP vs. GPP") +
#   scale_y_continuous(name="",limits = c(2.3, 3.5)) +
#   scale_x_continuous(name="") +
#   scale_shape_manual(values=c(19,1),name="Simulation") +
#   theme_bw() +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=10, face="bold"),
#         plot.title = element_text(size=14,face="bold"))
# 
# p3<-ggplot(data=frame[ind,], aes(x=PUP/GPP, y=GPP/1000, shape=SIM)) +
#   geom_point(show.legend=F,size=3) + 
#   labs(title = "CNP models: Annual PUP/GPP vs. GPP") +
#   scale_y_continuous(name="",limits = c(2.3, 3.5)) +
#   scale_x_continuous(name="") +
#   scale_shape_manual(values=c(19,1),name="Simulation") +
#   theme_bw() +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=10, face="bold"),
#         plot.title = element_text(size=14,face="bold"))
# 
# p4<-ggplot(data=frame[ind,], aes(x=PUP/NUP, y=GPP/1000, shape=SIM)) +
#   geom_point(size=3) + 
#   labs(title = "CNP models: Annual PUP/NUP vs. GPP") +
#   scale_y_continuous(name="",limits = c(2.3, 3.5)) +
#   scale_x_continuous(name="") +
#   scale_shape_manual(values=c(19,1),name="Simulation") +
#   theme_bw() +
#   theme(strip.text=element_text(size=12, face="bold"),
#         axis.text=element_text(size=10, face="bold"),
#         legend.position = "bottom",
#         legend.title=element_text(size=14,face="bold"),
#         legend.text=element_text(size=12,face="bold"),
#         plot.title = element_text(size=14,face="bold"))
# 
# pdf("Suppl_Fig3B_NUP_PUP_GPP_modelgroup.pdf",width=8,height=6)
# grid.arrange(arrangeGrob(p1, ncol=1),
#              arrangeGrob(p2, ncol=1),
#              arrangeGrob(p3, ncol=1),
#              arrangeGrob(p4, ncol=1), heights=c(1,1))
# dev.off()
# #for CFR the story is a bit complicated.. in relative terms CFR is increasing by more than CBio, thus 
# #more goes there, but if we do CFR/CBIO, the difference is marginal
# 
# frame<-melt(fr7_relgr[,names(fr7_relgr)!=("MOD")],id=c("MODgr","PRD")) #the fix again
# levels(frame$PRD) <- c("1st YEAR", "15YEARS", "CHANGE")
# 
# ind<-frame$MODgr%in%c("CN models")
# ind2<-frame$variable%in%c("CFR_CBIO")
# 
# p2<-ggplot(data=frame[ind==T&ind2==T,], aes(x=variable, y=mean, fill=SIM)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8, show.legend=FALSE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   scale_y_continuous(name="") +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c("Grey80","Black")) +
#   labs(title = "Responses to CO2 \nCN models") +
#   theme_bw() +
#   theme(strip.text=element_text(size=16, face="bold"),
#         axis.text=element_text(size=16, face="bold"),
#         axis.text.x=element_text(angle=45, vjust=0.5),
#         axis.title=element_text(size=16,face="bold"),
#         plot.title = element_text(size=16,face="bold"))
# 
# ind<-frame$MODgr%in%c("CNP models")
# ind2<-frame$variable%in%c("CFR_CBIO")
# 
# p4<-ggplot(data=frame[ind==T&ind2==T,], aes(x=variable, y=mean, fill=SIM)) +
#   geom_bar(stat="identity", position=position_dodge(), width=0.8, show.legend=FALSE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   scale_y_continuous(name="") +
#   scale_x_discrete(name="") +
#   scale_fill_manual(values=c("Grey80","Black")) +
#   labs(title = "Responses to CO2 \nCN models") +
#   theme_bw() +
#   theme(strip.text=element_text(size=16, face="bold"),
#         axis.text=element_text(size=16, face="bold"),
#         axis.text.x=element_text(angle=45, vjust=0.5),
#         axis.title=element_text(size=16,face="bold"),
#         plot.title = element_text(size=16,face="bold"))
# 
# grid.arrange(arrangeGrob(p1, ncol=1),
#              arrangeGrob(p2, ncol=1),
#              arrangeGrob(p3, ncol=1),
#              arrangeGrob(p4, ncol=1), heights=c(1,1)) 
# 



