# Data cleaner for old EA
# Correct N drift for each EA run individually
# Basic idea: at start of run (position 0), N_measured = N_actual
# As the run progresses, N_meas increases or decreases
#   by some proportion of N_act
# e.g. at position 30, N_meas=N_act*1.08
1.08/30
# N_meas = N_act*(1+0.036*line number)
# N_act = N_meas/(1+ncf*line) where ncf=slope as proportion of intercept
# Use standards to measure ncf (N correction factor)
# Within each run, need to make a judgment as to start & stop of
#   linear relationship between position and N_meas
# Don't correct if change during run is < 5%?
setwd('C:\\Users\\Devin\\Documents\\soil_analysis')
source('sample_ID_splitter.R')

setwd('C:\\Users\\Devin\\Documents\\soil_analysis\\EA data')

prepEAdat=function(dfr){
  if(is.element('X',names(dfr))){dfr=dfr[,-which(names(dfr)=='X')]}
  names(dfr)=c('run.date','Wt.mg','ID','N','N.mg','N.area',
             'C','C.mg','C.area')
  dfr$ID=gsub(' ','.',dfr$ID)
  dfr=strfun(dfr)
  unique(dfr$site)
  laststd=min(which(dfr$site!='std'))-1
  dfr$line=c(rep(0,laststd),seq(1,length(dfr$ID)-laststd))
  dfr
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

aug4=read.csv('Aug_04_17.csv')
#head(aug4)
aug4=prepEAdat(aug4)
plotEAdat(aug4)
aug4=fixEAdat(aug4)
plotnewEAdat(aug4) # minimal change

jul14=read.csv('Jul_14_17.csv')
jul14=prepEAdat(jul14)
plotEAdat(jul14)
jul14=jul14[jul14$line<55,] #machine conked out after that
#jul14=fixEAdat(jul14)
#plotnewEAdat(jul14)
# Super linear relationship between std N and position no. 
#   but only after first standard (line = 13)
jul14lm=lm(N~line,data=jul14[jul14$site=='std' & jul14$line>12,])
ncf14=summary(jul14lm)$coef[2]/jul14$N[jul14$line==13]
jul14$Nadj=jul14$N
jul14$Nadj[jul14$line>12]=jul14$N[jul14$line>12]/
  (1+ncf14*(jul14[jul14$line>12,]$line-13))
plotnewEAdat(jul14)
# maybe that original linearity was just by chance, but go with this.

nov01=read.csv('Nov_01_17.csv') 
nov01=prepEAdat(nov01)
plotEAdat(nov01) # Nice
nov01=fixEAdat(nov01)
plotnewEAdat(nov01)

jul11=read.csv('Jul_11_17.csv')
jul11=prepEAdat(jul11)
plotEAdat(jul11)
jul11$site[jul11$site=='std'&jul11$N<3]=NA # not sure what happened there
# I think the mass was just wrong
jul11=jul11[jul11$line<55,] # Guess 54 is how many samples you can run
plotEAdat(jul11) # good enough
jul11=fixEAdat(jul11)
plotnewEAdat(jul11) # still kind of high, but ok I guess
# pretty minimal correction

jun14=read.csv('Jun_14_17_fixed.csv')
jun14=prepEAdat(jun14)
plotEAdat(jun14) # 3 points = super linear
jun14=fixEAdat(jun14)
plotnewEAdat(jun14)
# hm
jun14$ID[jun14$site=='TM']
jun14$N[jun14$site=='TM']
jun14$Nadj[jun14$site=='TM']
sd(c(.1108,.1213,.1501))/mean(c(.1108,.1213,.1501))
# 0.1597153
sd(c(.0887,.0815,.1001))/mean(c(.0887,.0815,.1001))
# 0.1040921 SE of replicates got smaller, anyway 
# but they got worse for Cr, ugh
sd(c(.2136,.2286,.2239))/mean(c(.2136,.2286,.2239))
sd(c(0.17906404, 0.15602440, 0.15164233))/
  mean(c(0.17906404, 0.15602440, 0.15164233))
# still smallish though (went from 3 to 9%)

#jul07=read.csv('Jul_07_17.csv')
#jul07=prepEAdat(jul07)
#jul07$site[jul07$N>9]=NA
#jul07$site[jul07$N<4 & jul07$site=='std']=NA
#plotEAdat(jul07)
# Nah there is no relationship with position, 
#   just some crazy standard values

oct23=read.csv('Oct_23_17_dubious.csv')
oct23=prepEAdat(oct23)
plotEAdat(oct23)
oct23=oct23[oct23$line<20,]
oct23=fixEAdat(oct23)
plotnewEAdat(oct23) 
# go with it I guess

jun29=read.csv('Jun_29_17.csv')
jun29=prepEAdat(jun29)
plotEAdat(jun29) # N goes way off at end of run
jun29lm=lm(N~line,data=jun29[jun29$site=='std' & jun29$line>26,])
ncf29=summary(jun29lm)$coef[2]/jun29$N[jun29$line==27]
jun29$Nadj=jun29$N
jun29$Nadj[jun29$line>26]=jun29$N[jun29$line>26]/
  (1+ncf29*(jun29[jun29$line>26,]$line-27))
plotnewEAdat(jun29) 
# That's probably not quite right, but the duplicates are now closer 

aug2=read.csv('Aug_2_17_nofixneeded.csv')
aug2=prepEAdat(aug2)
jul28=read.csv('Jul_28_17_nofixneeded.csv')
jul28=prepEAdat(jul28)
jul07=read.csv('Jul_07_17_nofixneeded_errsremoved.csv')
jul07=prepEAdat(jul07)
jul20=read.csv('Jul_20_17_nofix.csv')
jul20=prepEAdat(jul20)
jul26=read.csv('Jul_26_17_nofix.csv')
jul26=prepEAdat(jul26)


alladj=rbind(aug4,jul14,jul11,nov01,jun14,jun29,oct23)
alloth=rbind(aug2,jul28,jul07,jul20,jul26)
alloth$Nadj=alloth$N
alladj=alladj[,which(names(alladj)%in% names(alloth))]
allEA=rbind(alladj,alloth)
unique(allEA$site)

setwd('C:\\Users\\Devin\\Documents\\soil_analysis')
source('XRF_formatting_functions.R')

allEA$ID=typofixer(allEA$ID)
allEA$ID=gsub('_2','',allEA$ID)
allEA$ID=gsub('_3','',allEA$ID)
allEA$ID=gsub('_1','',allEA$ID)
allEA$ID=gsub('hr','',allEA$ID)
allEA$ID=gsub('_0','',allEA$ID)
allEA$ID=gsub('_4','',allEA$ID)
allEA=allEA[!is.na(allEA$site),]
allEA$ID[allEA$site=='TM'|allEA$site=='Cr']=paste(
   allEA$ID[allEA$site=='TM'|allEA$site=='Cr'],'16',sep='.'
 )
unique(allEA$year)
# remake year column with substitutions
allEA=strfun(allEA)
# Add position to stands that were only sampled in one position in 2004
allEA$ID[allEA$stand=='It.E2'&
         allEA$year=='04']=gsub('It.E2','It.E2.L',
                              allEA$ID[allEA$stand=='It.E2'&allEA$year=='04'])
allEA$ID[allEA$stand=='It.E1'&
         allEA$year=='04']=gsub('It.E1','It.E1.T',
                              allEA$ID[allEA$stand=='It.E1'&allEA$year=='04'])
allEA$ID[allEA$stand=='JP.E1'&
         allEA$year=='04']=gsub('JP.E1','JP.E1.E',
                              allEA$ID[allEA$stand=='JP.E1'&allEA$year=='04'])
allEA$ID[allEA$stand=='JP.E2'&
           allEA$year=='04']=gsub('JP.E2','JP.E2.E',
                                  allEA$ID[allEA$stand=='JP.E2'&allEA$year=='04'])
allEA$ID[allEA$stand=='Bp.E2'&
         allEA$year=='04']=gsub('Bp.E2','Bp.E2.E',
                              allEA$ID[allEA$stand=='Bp.E2'&allEA$year=='04'])

allEA=strfun(allEA)
allEA=droplevels(allEA[allEA$site %in% c('Bp','BO','JP','It','Vg','Eu',
                              'Cr','TM'),])

# get month and day--all are in 2017
EAstrfun=function(dat){
  strs=strsplit(as.character(dat$run.date),'-')
  firsts=as.character(lapply(strs,function(L){L[1]}))
  seconds=as.character(lapply(strs,function(L){L[2]}))
  dat$CNevalday=paste(seconds,firsts,sep='-')
  dat
}

allEA=EAstrfun(allEA) 
head(allEA)
unique(allEA$year)

#write.csv(allEA,'2017_EA_data.csv')
