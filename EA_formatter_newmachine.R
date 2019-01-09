#setwd('C:\\Users\\Devin\\Documents\\Soil data')
#source('sample_ID_splitter.R')
source('XRF_formatting_functions.R')

prepnewEAdat=function(dfr){
  dfr=dfr[,c(1,3,5,6,7,8,9)]
  names(dfr)=c('run.date','ID','Wt.mg','N.mg','C.mg','N','C')
  dfr=dfr[-1,]
  dfr1=dfr[,(1:2)]
  dfr2=as.data.frame(lapply(dfr[,-(1:2)],makenumeric))
  newdfr=cbind(dfr1,dfr2)
  head(newdfr)
  newdfr$ID=gsub(' ','.',newdfr$ID)
  newdfr$ID=typofixer(newdfr$ID)
  newdfr=strfun(newdfr)
  laststd=min(which(newdfr$site!='std' & newdfr$site!='bypass'))-1
  newdfr$line=c(rep(0,laststd),seq(1,length(newdfr$ID)-laststd))
  newdfr
}
prepnewEAdat2=function(dfr){
  dfr=dfr[,c(1,3,5,8,9,10,11)] 
  # not sure why Area columns exported only in September
  names(dfr)=c('run.date','ID','Wt.mg','N.mg','C.mg','N','C')
  dfr=dfr[-1,]
  dfr1=dfr[,(1:2)]
  dfr2=as.data.frame(lapply(dfr[,-(1:2)],makenumeric))
  newdfr=cbind(dfr1,dfr2)
  head(newdfr)
  newdfr$ID=gsub(' ','.',newdfr$ID)
  newdfr$ID=typofixer(newdfr$ID)
  newdfr=strfun(newdfr)
  laststd=min(which(newdfr$site!='std' & newdfr$site!='bypass'))-1
  newdfr$line=c(rep(0,laststd),seq(1,length(newdfr$ID)-laststd))
  newdfr
}

plotEAdat=function(dfr){
  intc=mean(dfr$N[dfr$line==0]) # ostensible intercept of regression
  plot(N~line,data=dfr[dfr$site=='std',])
  dfrlm=lm(N~line,data=dfr[dfr$site=='std',])
  abline(summary(dfrlm)$coef[1],summary(dfrlm)$coef[2])
  # good enough
  # should the N correction factor be coef2/coef1 or coef2/intc?
  # I think coef2/intc, but they're very close here
}
fixEAdat=function(dfr){
  intc=mean(dfr$N[dfr$line==0]) # ostensible intercept of regression
  dfrlm=lm(N~line,data=dfr[dfr$site=='std',])
  ncf=summary(dfrlm)$coef[2]/intc
  dfr$Nadj=dfr$N/(1+ncf*dfr$line)
  dfr
} 

plotnewEAdat=function(dfr){
  print(dfr$Nadj[dfr$site=='std']) # better, anyway
  plot(Nadj~N,data=dfr[dfr$site!='std',],col=stand)
  abline(0,1) 
}

ea1 = read.table('tests_7-12-18.txt',header=T,sep='\t',
                 as.is=T)
#head(ea1)
#sapply(ea1,class)
ea1=prepnewEAdat(ea1)
head(ea1)
ea1=strfun(ea1)
unique(ea1$site)
plotEAdat(ea1)
ea1=fixEAdat(ea1)
plotnewEAdat(ea1[ea1$site!='bypass',]) # very small correction

# for subsequent analyses, add standards to the curve throughout the run
#   rather than correcting for N drift
# averages the error over the whole run
# if it looks bad, can also remove those standards and correct
#   (decide before exporting results)

ea2=read.table('BOPItEJPN_7-20-18.txt',header=T,sep='\t',
               as.is=T)

ea2=prepnewEAdat(ea2)
ea2$Nadj=ea2$N

ea3=read.table('BOE_JPE2plus_7-25-18_withlaststd.txt',header=T,sep='\t',
               as.is=T)

ea3=prepnewEAdat(ea3)
ea3$Nadj=ea3$N
plotEAdat(ea3)

ea4=read.table('EA_7-27-18_fixedtypos.txt',header=T,sep='\t',
               as.is=T)
ea4=prepnewEAdat(ea4)
ea4$Nadj=ea4$N

ea5=read.table('8-7-18_EAresults.txt',header=T,sep='\t',
               as.is=T)
ea5=prepnewEAdat(ea5)
ea5$Nadj=ea5$N

ea6=read.table('EA_9-7-18_noadam.txt',header=T,sep='\t',
               as.is=T)
ea6=prepnewEAdat2(ea6)
ea6$Nadj=ea6$N

ea7=read.table('EA_9-11-18.txt',header=T,sep='\t',
               as.is=T)
ea7=prepnewEAdat2(ea7)
ea7$Nadj=ea7$N

ea8=read.table('EA_9-13-18_ok.txt',header=T,sep='\t',
               as.is=T)
ea8=prepnewEAdat2(ea8)
ea8$Nadj=ea8$N

ea9=read.table('EA_9-10-18.txt',header=T,sep='\t',
               as.is=T)
ea9=prepnewEAdat2(ea9)
ea9$Nadj=ea9$N

ea10=read.table('EA_10-2-18.txt',header=T,sep='\t',
               as.is=T)
ea10=prepnewEAdat2(ea10)
ea10$Nadj=ea10$N

ea12=ea1[,which(names(ea1)%in% names(ea2))]
newea=rbind(ea12,ea2,ea3,ea4,ea5,ea6,ea7,ea8,ea9,ea10)
unique(newea$site)
 
newea$ID=typofixer(newea$ID)
newea$ID=gsub('_2','',newea$ID)
newea$ID=gsub('_3','',newea$ID)
newea$ID=gsub('_1','',newea$ID)
newea$ID=gsub('_4','',newea$ID)
# Add position to stands that were only sampled in one position in 2004
newea$ID[newea$stand=='It.E2'&
         newea$year=='04']=gsub('It.E2','It.E2.L',
                              newea$ID[newea$stand=='It.E2'&newea$year=='04'])
newea$ID[newea$stand=='It.E1'&
         newea$year=='04']=gsub('It.E1','It.E1.T',
                              newea$ID[newea$stand=='It.E1'&newea$year=='04'])
newea$ID[newea$stand=='JP.E1'&
         newea$year=='04']=gsub('JP.E1','JP.E1.E',
                              newea$ID[newea$stand=='JP.E1'&newea$year=='04'])
newea$ID[newea$stand=='Bp.E2'&
         newea$year=='04']=gsub('Bp.E2','Bp.E2.E',
                              newea$ID[newea$stand=='Bp.E2'&newea$year=='04'])

# looked at difference in composites with and without L (apparent OM contam)
# same N (0.008%), C is .124% with L vs .107% without
# range for other 3 reps is .87 to .11 % C
# just leave the one with all four positions for consistency with other depths
newea=newea[-which(newea$ID=='Eu.E2.1.40-60.16_noL'),]


newea=strfun(newea)
newea=droplevels(newea[newea$site %in% c('Bp','BO','JP','It','Vg','Eu',
                                         'Cr','TM'),])
head(newea$run.date)
EAstrfun2=function(dat){
  strs=strsplit(as.character(dat$run.date),'_')
  mo=as.character(lapply(strs,function(L){L[2]}))
  day=as.character(lapply(strs,function(L){L[3]}))
  dat$CNevalday=paste(mo,day,sep='-')
  dat
}
newea=EAstrfun2(newea)
#table(newea$CNevalday)
# fix a few samples where a repeat messed up the file name
newea$CNevalday[newea$CNevalday=='2-9'&newea$stand=='Cr.A1']='9-11'
newea$CNevalday[newea$CNevalday=='2-9']='9-7'
#write.csv(newea,'new_CN_7-12and20.csv')
#write.csv(newea,'new_CN_7-12to25-18.csv')
#write.csv(newea,'CN18_to7-27.csv')
#saveRDS(newea,'CN18_to7-27.Rds')
#write.csv(newea,'CN18_to8-17.csv')
#write.csv(newea,'CN18_to9-13.csv')
#write.csv(newea,'CN18_to_10-2.csv')
# 1-8-19
#write.csv(newea,'CN18_to_10-2_remake.csv')
