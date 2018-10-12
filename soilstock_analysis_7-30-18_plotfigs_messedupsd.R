# set your working directory (or open a project containing all the files)
setwd('C:\\Users\\Devin\\Documents\\Soil data')

# useful packages:
library(dplyr)
library(lattice)
library(nlme)
library(reshape)
library(RColorBrewer)
source('xrf_formatting_functions.R')

#widedats=readRDS('all_data_7-30-18.Rds')
#widedats=readRDS('all_data_8-17-18.Rds')
widedats=readRDS('all_data_9-13-18.Rds')
widedats$CN=widedats$C/widedats$N # this ratio can be useful, or just noisy
#xyplot(depth~CN|stand,data=widedats,ylim=c(90,0),xlim=c(0,40))
# C:N decreases with depth? why? burning off C to get N which stays in OM?
widedats1=widedats[,-which(names(widedats) %in% 
                             c('stand','depth','year','site',
                               'LU','position','rep','elt',
                               'inc_from','inc_to','newest','oldmeas',
                               'P_ppm','Fe_pct','Ca_ppm',
                               'K_ppm','Mg_ppm'))]

# replace NAs with detection limit
widedats1=mutate(widedats1,Ca2=ifelse(is.na(Ca) & !is.na(Ca_dl),
                                      Ca_dl,Ca),
                 P2=ifelse(is.na(P) & !is.na(P_dl),P_dl,P),
                 Mg2=ifelse(is.na(Mg) & !is.na(Mg_dl),Mg_dl,Mg))

# seperate ppm columns are unneccessary
# "Melt" the data frame into a long format
dats=melt(widedats1,id.vars=c('ID','run.date','CNevalday','CNavgd',
                              'measured_in','Method',
                              'Eval.Date','evaldate',
                              'xrfavgd','recent','avgBD','BDsd'),
              na.rm=F,variable_name = 'element')
# Restore the columns made from the ID
dats=strfun(dats)
dats=eltfun(dats)
dats$stdonly=as.character(lapply(strsplit(as.character(dats$stand),'[.]'),
                              function(x){x[2]}))
dats$unit=as.character(lapply(strsplit(as.character(dats$element),'_'),
                                 function(x){x[2]}))
dats$unit[is.na(dats$unit)]='ppm'
dats$unit[dats$element %in% 
               c('Al','Si','C','N','N.old','C.WB','Fe')]='pct'
dats$value[dats$element %in% 
               c('Mg','P','Mg2','P2','Ca2','S','Cl','K','Ca')]=
  dats$value[dats$element %in% 
               c('Mg','P','Mg2','P2','Ca2','S','Cl','K','Ca')]*10000
dats$value[dats$element=='Fe']= dats$value[dats$element =='Fe']/10000
#dats=readRDS('all_data_long_7-22-18.Rds')
#unique(dats$element)

# Take out crazy outliers
#plot(value~as.numeric(stand),data=dats[dats$element=='Ca',],col=year)
#dats$ID[dats$element=='Ca'& dats$value>5000] #JP.E2.T.2.0-10.16
# Also two of Eu.E1
dats$value[dats$ID=='JP.E2.T.2.0-10.16'&
             dats$element %in% c('Ca','Mg')]=NA 
# piece of dolomite got in here, maybe

dats$value[dats$element=='P' & dats$value>600 & 
             dats$stand!='Eu.E1']=NA
# a little more subjective-- three points >> others in their stand 


JPAs=dats[dats$stand=='JP.A',]
plot(value~depth,data=JPAs[JPAs$element=='C',],col=position)
E2s=droplevels(dats[dats$stand=='Eu.E2',])
plot(value~as.numeric(elt),col=rep,
     data=E2s[E2s$element=='K'&E2s$depth==5,])
# Two crazy high potassium values in toco
# lots of P bdl; all < 50 ppm
# one high Ca also in toco, but some Ees high, too
# can't measure the Mg here
plot(value~as.numeric(elt),data=dats[dats$element=='C' & dats$depth==5&
                             dats$stand=='It.E1',],
     pch=as.numeric(year),col=rep)




# Average the row positions to get 4 reps per site (5 if there's a soil pit measurement)
# This way, we aren't weighting the individual cores and the composited cores the same
dats4=group_by(dats,stand,year,depth,rep,element) %>%
  mutate(repval=mean(value,na.rm=T),repsd=sd(value,na.rm=T),
         repn=n(),repwt=ifelse(rep==5,2,1))
dats4=distinct(dats4,stand,year,depth,rep,element,.keep_all = T)
# In future stats, maybe weight pit samples more than core samples when both taken
# Pit samples were taken when cores dubious (i.e. contaminated with surface OM)

# An analogous version with the "wide" data, mostly useful if you want to compare
#   different elements within the same sample
widedats4=group_by(widedats,stand,year,depth,rep,site,LU) %>%
  summarise_if(is.numeric, mean,na.rm=T)

plot(C~Fe_pct,data=widedats4,col=year)
plot(C~N,data=widedats4,col=year) # generally pretty consistent
#   This suggests (but does not prove, of course) that most of the N in the sites
#     is in the form of organic matter, not inorganic N from fertilizer
#   We would expect this to be the case because this is a pretty low-N system
#     so plants and microbes will rapidly take up avaiable N into their biomass
#   The relationship is more tenuous at low values of N
#   You can come up with a statistical analysis of this relationship in the 
#     different depths, years, and vegetation types if you want

# Calculate pseudo-tau relative to Zr (or other element) in 2004
dats4=group_by(dats4,stand,rep,depth) %>%
  mutate(Zrrat=ifelse(sum(!is.na(repval[element=='Zr'& year=='04']))>0 &
                        sum(!is.na(repval[element=='Zr'& year=='16']))>0,
           repval[element=='Zr'&year=='16']/
           repval[element=='Zr'&year=='04'],NA))
summary(dats4$Zrrat)
unique(dats4$stand[dats4$Zrrat>2])
dats4=group_by(dats4,stand,rep,depth,element)%>%
  mutate(tau=ifelse(!is.na(Zrrat),
                    ((repval[year=='16']/repval[year=='04'])/Zrrat)-1,
                    NA))


##### Further analysis of different depths--skim this before focusing on a single
##### depth across all sites or some other composite metric
################
# Now that we have equivalent, composited replicates in each site, run all the t-tests
# Exclude the stands with incomplete data (Eu, JP.A)
#   and the extra XRF data (intensity, etc.))
ttests=group_by(droplevels(dats4[!is.element(dats4$unit,c('dl','err','int')) &
                        #!is.element(dats4$stand,c('Eu.N','Eu.E1','Eu.E2','JP.A'))&
                          !is.element(dats4$site,c('TM','Cr')) &
                          dats4$LU!='A' & !is.element(dats4$element,c('Mg2','P2')),]),
                stand,depth,element,unit,site,LU) %>%
  summarise(mn16=mean(repval[year=='16'],na.rm=T),
            mn04=mean(repval[year=='04'],na.rm=T),
            n16=sum(!is.na(repval[year=='16'])),
            n04=sum(!is.na(repval[year=='04'])),
            sd16=sd(repval[year=='16'],na.rm=T),
            sd04=sd(repval[year=='04'],na.rm=T),
            tstat=ifelse((n16>2 & n04 > 2),
                      t.test(repval[year=='16'],
                             repval[year=='04'])$statistic,NA),
            pval=ifelse((n16>2 & n04 > 2),
                        t.test(repval[year=='16'],
                               repval[year=='04'])$p.value,NA))

# Where did anything change?
table(ttests$stand[ttests$pval<.05 & !is.element(ttests$site,c('TM','Cr'))],
      ttests$depth[ttests$pval<.05& !is.element(ttests$site,c('TM','Cr'))])
# Something significant just about everywhere, but meaningful?

shorttests=droplevels(ttests[ttests$element %in% 
                              c('C','N','P','K','Ca','Mg'),])
# How many elements changed at a different stand and depth?
# Remember that p < 0.05 is arbitrary and may have no biological significance
#   but this arbitrary threshold helps us think about where measurable changes occurred
#   and how those changes might be correlated
#   and to focus on depths and stands of interest
table(shorttests$stand[shorttests$pval<.05], 
      shorttests$depth[shorttests$pval<.05])

# Add up the rows or columns of the table with apply()
apply(table(shorttests$stand[shorttests$pval<.05], 
            shorttests$depth[shorttests$pval<.05]),1,sum)
# Most changes: Bp.E1, then JP.N and E1, BO.P, all Its
apply(table(shorttests$stand[shorttests$pval<.05], 
            shorttests$depth[shorttests$pval<.05]),2,sum)
#  most changes in 0-10 and 20-40

# Where did N change?
table(ttests$stand[ttests$pval<.05 & ttests$element=='N'],
      ttests$depth[ttests$pval<.05& ttests$element=='N'])

# visualize the direction of changes by multiplying and adding matrices
# matrix of (which ones increased significantly)-matrix of (decreased signif.ly)
ttable=function(dfr){
table(dfr$depth[dfr$pval<.05],dfr$stand[dfr$pval<.05])*
  table(dfr$depth[dfr$tstat>0],dfr$stand[dfr$tstat>0])+
  table(dfr$depth[dfr$pval<.05],dfr$stand[dfr$pval<.05])*
  table(dfr$depth[dfr$tstat<0],dfr$stand[dfr$tstat<0])*-1
}
# 0: no significant change; 1: significant increase between 2004 and 2016, 
#   -1: significant decrease in concentration at this depth
# try a different alpha threshold: "almost significant", p<0.1
laxttable=function(dfr){
  table(dfr$depth[dfr$pval<.1],dfr$stand[dfr$pval<.1])*
    table(dfr$depth[dfr$tstat>0],dfr$stand[dfr$tstat>0])+
    table(dfr$depth[dfr$pval<.1],dfr$stand[dfr$pval<.1])*
    table(dfr$depth[dfr$tstat<0],dfr$stand[dfr$tstat<0])*-1
}
ttable(ttests[ttests$element=='N',])
laxttable(ttests[ttests$element=='N',])

gen_ttable=function(elmtcol,groupcol,pvalcol,tstatcol,a){
  table(elmtcol[pvalcol<a],groupcol[pvalcol<a])*
    table(elmtcol[tstatcol>0],groupcol[tstatcol>0])+
    table(elmtcol[pvalcol<a],groupcol[pvalcol<a])*
    table(elmtcol[tstatcol<0],groupcol[tstatcol<0])*-1
}

gen_ttable(ttests[ttests$element=='P',]$depth,
           ttests[ttests$element=='P',]$stand,
           ttests[ttests$element=='P',]$pval,
           ttests[ttests$element=='P',]$tstat,0.05)

gen_ttable(shorttests[shorttests$site=='It',]$depth,
           shorttests[shorttests$site=='It',]$element,
           shorttests[shorttests$site=='It',]$pval,
           shorttests[shorttests$site=='It',]$tstat,0.05) 
# this doesn't quite work as intended, because tstat table isn't all 0s and 1s

gen_ttable(ttests[ttests$depth==30,]$stand,
           ttests[ttests$depth==30,]$element,
           ttests[ttests$depth==30,]$pval,
           ttests[ttests$depth==30,]$tstat,0.05)


Centests=shorttests[shorttests$site %in% c('BO','Vg'),]
gen_ttable(Centests$depth,Centests$element,Centests$pval,
           Centests$tstat,0.05)
# N increased at depth, C decreased, K increased a lot, Ca at surface,
# P increased in some depths

# what increased at depth? 
#   N in several (BO.E and P, Bp.E1, It.E1, JP.N),
#   Also some P (BO.E, It.E1 and N) and K (BO.P, JP.N, both Vgs)
table(shorttests$stand[shorttests$depth>30 & shorttests$tstat>0 
                           &shorttests$pval<0.05],
      shorttests$element[shorttests$depth>30 & shorttests$tstat>0 
                         &shorttests$pval<0.05])
# what decreased (tstat <0)?
# less: K in Bp.E1 of course, also JP.E1; Mg in Bp.E1 and It.E2 (not legit?)
#   C in BO.P, Bp.E1, JP.P


# N: increase at most depths in It.E1 and JP.N, 
#   distribution becomes deeper in BO.E? 
# No depletion below 40 cm, just increase in some sites
# In contrast, C depleted at depth in both pastures
# Few sites have opposing signif changes at different depths; no strong evidence
#   for vertical redistribution nutrients


# For ease of visualization, get a single average value per stand, depth, 
#   year, and element
datsmean=group_by(dats4[!is.element(dats4$unit,c('dl','err','int')),],
                  stand,depth,year,element,unit,stdonly,site,LU,avgBD,BDsd) %>%
  summarise(mn=mean(repval,na.rm=T),sd=sd(repval,na.rm=T),
            taumn=mean(tau,na.rm=T),sdtau=sd(tau,na.rm=T))
datsmnok=droplevels(datsmean[datsmean$site!='TM' &
                               datsmean$site!='Cr'& datsmean$LU!='A',])

# Use the lattice functions to visualize changes at different depths
# Compare to the t-test tables to see what differences between years 
#   are statistically significant
datsmnok=datsmnok[order(datsmnok$depth),]
xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='C' ,],
       ylim=c(90,0),xlab='C (g / 100 g)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# Try different nutrients and subsets of stands; change units as needed
xyplot(depth~mn*10|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='N',],#&
                       #datsmnok$stand!='Bp.E2',],
       ylim=c(90,0),
       xlab='N (g / kg)',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# Bulk density vs depth
xyplot(depth~avgBD|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='N',],#&datsmnok$stand!='JP.P',],
                       ylim=c(90,0),
       xlab='Bulk density (g / cm3)',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
xyplot(depth~mn/avgBD|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='C',],ylim=c(90,0),
       xlab='C/Bulk density',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
plot(avgBD~mn,data=datsmnok[datsmnok$element=='C',],col=site,
     pch=as.numeric(year),xlab='Mean C content, % (depth and site)',
     ylab='Mean bulk density, g cm-3',las=1)
legend('topright',bty='n',col=c(as.numeric(unique(datsmnok$site)),1,1),
       pch=c(rep(16,6),1,2),legend=c(unique(as.character(datsmnok$site)),
                                     '04','16'))

# Don't talk about Mg, it's a mess
# Ca from lower depths is dropped?
# K increases at all depths in BO.P, decreases in Bp.E1
#   big K changes in Bp.E2 are not significant

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca2'&datsmnok$stand!='Eu.E1',],
       ylim=c(90,0),
       xlab='Ca (mg/kg)',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# Lots of missing data
# Ca in JP goes deeper in the pasture than in the eucalyptus 
#  (need to check on timing/amount of application; could indicate different rates 
#  of leaching vs. uptake by the plants to keep it at the surface)
xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='K'&datsmnok$stand!='Bp.E1',],ylim=c(90,0),
       xlab='K (mg/kg)',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='P2',],ylim=c(90,0),
       xlab='P (mg/kg)',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# P is maybe the most shallowly distributed
# Added in fertilizer and doesn't move, or actively cycled?

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Fe' ,],
       ylim=c(90,0),xlab='Fe (g / 100 g)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))


# Looking at pseudo-tau with depth
xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Zr' ,],
       ylim=c(90,0),xlab='Zr (mg / kg)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~taumn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='K'&datsmnok$stand!='Bp.E1',],
       ylim=c(90,0), xlab='K pseudo-tau',as.table=T,#layout=c(4,3),
       par.settings = list(superpose.line = list(col = 'gray20',lwd = 2)))
xyplot(depth~taumn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='C',],
       ylim=c(90,0), xlab='C pseudo-tau',as.table=T,
       par.settings = list(superpose.line = list(col = 'gray20',lwd = 2)))
xyplot(depth~taumn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='P2',],
       ylim=c(90,0), xlab='P pseudo-tau',as.table=T,
       par.settings = list(superpose.line = list(col = 'gray20',lwd = 2)))


# Increases in It.N: feasible?
#Its=datsmnok[datsmnok$site=='It',]
Its=datsmnok[datsmnok$stand %in% c('It.E1','It.N'),]
#Its$stand=factor(Its$stand,levels=c('It.E2','It.E1','It.N'))
Its$stand=factor(Its$stand,levels=c('It.E1','It.N'))
bgColors=c('blue3','blue3','springgreen')
myStripStyle <- function(which.panel, factor.levels, ...) {
  panel.rect(0, 0, 1, 1,
             col = bgColors[which.panel],
             border = 1)
  panel.text(x = 0.5, y = 0.5,
             lab = c('Low-OM Euc','High-OM Euc',
                     'High-OM Cerrado')[which.panel],
             col = 'gray90')
}  
trellis.par.set(strip.background = list(col = 'grey80'),
                par.strip.text=list(cex=.8))
xyplot(depth~mn*10|stand,groups=year,type='l',ylab='Depth (cm)',
       data=Its[Its$element=='N',],ylim=c(90,0),
       xlab='N (g / kg)',as.table=T,layout=c(2,1),#layout=c(3,1),
       par.settings = list(superpose.line = list(
         col = c('orange','darkslateblue'),lwd = 2),
                           fontsize=list(text=20)),
       #scales=list(cex=1.5), #that's relative to fontsize, keep as is
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE,
                     cex=.8),
       #strip.background = list(col = c('blue3','blue3','springgreen')),
      #strip=strip.custom(factor.levels=c('Euc A','Euc B',
      #                                   'Cerrado')))
      strip=strip.custom(factor.levels=c('Eucalyptus','Cerrado')))
xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       #data=Its[Its$element=='C',],xlab='C (g / 100 g)',
       #data=Its[Its$element=='K',],xlab='K (mg / kg)',
       data=Its[Its$element=='P',],xlab='P (mg / kg)',
       ylim=c(90,0),as.table=T,layout=c(2,1),#layout=c(3,1),
       par.settings = list(superpose.line = list(
         col = c('orange','darkslateblue'),lwd = 2),
         fontsize=list(text=20)),
       strip=strip.custom(factor.levels=c('Eucalyptus','Cerrado')))
#       strip=strip.custom(factor.levels=c('Euc A','Euc B',
 #                                         'Cerrado')))


xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=Its[Its$element=='Ca',],ylim=c(90,0),
       xlab='Ca (mg / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=Its[Its$element=='K',],ylim=c(90,0),
       xlab='K (mg / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
xyplot(depth~avgBD|stand,groups=year,type='l',ylab='Depth (cm)',
       data=Its[Its$element=='K',],ylim=c(90,0),
       xlab='Bulk density (g/cm3)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# add points at 4,depth where change is significant?

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca'&
                       datsmnok$stand%in% c('JP.E2','JP.N','It.E1'),],
       ylim=c(90,0),
       xlab='Ca (mg / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(
         col = c('orange','darkslateblue'),lwd = 2),
         fontsize=list(text=20)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE,
                     cex=.8),
       strip=strip.custom(factor.levels=c('Euc B','Euc C',
                                          'Cerrado')))

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca'&
                       datsmnok$stand%in% c('JP.E2','JP.N','It.E1'),],
                       #datsmnok$stand%in% c('JP.E2','It.N','It.E1'),],
ylim=c(90,0),xlab='Ca (mg/kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))


abdatsmn=droplevels(datsmean[is.element(datsmean$site, c('TM','Cr','JP'))&
                               !is.element(datsmean$stand,c('JP.P'))&
                               datsmean$year==16,])
xyplot(depth~mn|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='C' ,],
       ylim=c(90,0),xlab='C (g / 100 g)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))


# End of depth distribution plots
################



# Bulk density
#	Convert concentrations to stocks (Mg ha-1) for a given depth increment
dats = group_by(dats,stand,depth) %>% 
  mutate(oldBD=mean(avgBD[year=='04']),newBDsd=mean(BDsd[year=='16']),
         BD16=mean(avgBD[year=='16']))

# For JP.P, assume no BD change below 20 cm to get stocks despite
#   inexplicably missing bulk density data
dats$avgBD[dats$stand=='JP.P'&dats$depth>20&dats$year=='04']=
  dats$BD16[dats$stand=='JP.P'&dats$depth>20&dats$year=='04']

datsstk=group_by(dats[!is.element(dats$unit,c('dl','err','int')),],
                 stand,inc_to,year,rep,element,
                 site,LU,avgBD,oldBD,newBDsd) %>%
  mutate(inc=inc_to-inc_from,  repval=mean(value,na.rm=T),repvar=var(value,na.rm=T),
            stock=ifelse(unit=='pct',repval*avgBD*inc,repval*avgBD*inc*.0001),
         # stock/area = concentration of element * density of soil * depth of soil
         #   stockunit=ifelse(unit=='pct','Mg.ha','kg.ha'),
         stockunit='Mg/ha',
         # variation in bulk density: an attempt at uncertainty quantification
         # only the bulk density values from 2016 have standard deviation
            stocklo=ifelse(unit=='pct',
                        repval*(avgBD-newBDsd)*inc,repval*(avgBD-newBDsd)*inc*.0001),
            stockhi=ifelse(unit=='pct',
                        repval*(avgBD+newBDsd)*inc,repval*(avgBD+newBDsd)*inc*.0001),
         # another way: suppose bulk density didn't change since 2004
            oldBDstock=ifelse(unit=='pct',repval*oldBD*inc,repval*oldBD*inc*.0001))
datsstk=distinct(datsstk,stand,inc_to,year,rep,element,.keep_all=T) # one per rep
datsstk=group_by(datsstk,stand,inc_to,year,element) %>%
  mutate(depmnsq=mean(repval,na.rm=T)^2,depvar=var(repval,na.rm=T),
         stockvar=depmnsq*(depvar/depmnsq+newBDsd^2/avgBD^2))

# For the paper, create more representative confidence intervals/distributions of 
#   possible stocks given the uncertainty in both bulk density and concentration measurements
#   by bootstrapping (resampling many times with replacement) different combinations
#   of the bulk density and concentration values measured in a given stand, year, and depth
#   prior to adding up the stocks. Not necessary for now.

# Summarize at two depths: 20 cm (approximately the A horizon, maximum root influence)
#   and 100 cm
# Add up the layers of soil to get the stock
# variance of stock = stock^2 * (var(conc)/conc^2 + var(BD)/BD2 +
#  2*cov(conc-BD)/conc*BD)
examp1=datsstk[datsstk$stand=='It.N'&datsstk$inc_to==20&
                 datsstk$element=='N'&datsstk$year=='16',]
cov(examp1$avgBD,examp1$repval) #0 because just one avgBD
# just skip that part? 1st-order expansion? Do this later?

dats2deps=group_by(datsstk[!is.element(datsstk$site,c('TM','Cr'))&
                             datsstk$stand!='JP.A',],stand,year,rep,
                   element,site,LU,unit,stockunit) %>% 
  summarise(ndeps20=length(unique(inc_to[inc_to<=20])),
            stock20=ifelse(ndeps20==2,sum(stock[inc_to<=20]),NA),
            stocklo20=ifelse(ndeps20==2,sum(stocklo[inc_to<=20]),NA),
            stockhi20=ifelse(ndeps20==2,sum(stockhi[inc_to<=20]),NA),
            oldBDstock20=ifelse(ndeps20==2,sum(oldBDstock[inc_to<=20]),NA),
            conc20=ifelse(ndeps20==2,
                          sum(repval[inc_to<=20]*inc[inc_to<=20]*avgBD[inc_to<=20])/
                            sum(inc[inc_to<=20]*avgBD[inc_to<=20]),NA),
            # weight concentrations by mass of soil in each layer
            # to get average density of the whole 1-m block of soil
            BD20=ifelse(ndeps20==2,
                       sum(avgBD[inc_to<=20]*inc[inc_to<=20])/
                             sum(inc[inc_to<=20]),NA),
            #BDsd20=ifelse(ndeps20==2,
            #           sum(BDsd[inc_to<=20]*inc[inc_to<=20])/
            #             sum(inc[inc_to<=20]),NA),
            BDsd20=ifelse(ndeps20==2,sqrt(sum(BDsd[inc_to<=20])^2),NA),
            stocksd20=ifelse(ndeps20==2,sqrt(sum(stockvar[inc_to<=20])),NA),
            # Variance of sum = sum of variances; check how best to do this
            ndeps100=length(unique(inc_to)),
            stock100=ifelse(ndeps100==5,sum(stock),NA),
            stocklo100=ifelse(ndeps100==5,sum(stocklo),NA),
            stockhi100=ifelse(ndeps100==5,sum(stockhi),NA),
            oldBDstock100=ifelse(ndeps100==5,sum(oldBDstock),NA),
            conc100=ifelse(ndeps100==5,
                          sum(repval*inc*avgBD)/sum(inc*avgBD),NA),
            conc60to100=mean(repval[inc_to==100],na.rm=T),
            BD100=ifelse(ndeps100==5,sum(avgBD*inc)/sum(inc),NA),
            #BDsd100=ifelse(ndeps100==5,sqrt(sum(BDsd^2)),NA),
           # weight by inc
           suminc=sum(inc),
           BDsd100=ifelse(ndeps100==5,sqrt(sum((BDsd^2)*inc)/sum(inc)),NA),
           stocksd100=ifelse(ndeps100==5,sqrt(sum(stockvar*inc)/sum(inc)),NA),
           stockratio=stock20/stock100,concratio=conc20/conc100,
           concrat2=conc20/conc60to100)

hist(dats2deps$stockratio) # > 20% = surface enrichment, less = surface depletion
# (but how much less is significant will depend on the variability of the measurements)
plot(concratio~concrat2,data=dats2deps[dats2deps$element=='P',])
# more accurate depiction of superficiality, 
#   or more prone to outlier weirdness?

# A pretty straightforward analysis: are there directional trends in concentrations
#   of different nutrients (could modify it to look at stocks and 100 cm instead of 20)
#   between eucalyptus and other vegetation, and between years?
##########
# Are the inter-year trends different in different vegetation types (interaction term)?
simple20=droplevels(dats2deps[dats2deps$stand %in% 
                                c('BO.E','BO.P','Vg.E','Vg.N',
                                  'JP.E1','JP.N','It.E1','It.N'),])
simple20=mutate(simple20,LU2=ifelse(LU=='E','E','O'))
# for this analysis, O = "other"
table(simple20$LU2[simple20$element=='C'],simple20$year[simple20$element=='C']) 
# Not balanced--rep 5
Ccsimp.lme=lme(conc20~year*LU2,random=~1|site,
               data=simple20[simple20$element=='C',], na.action=na.omit)
summary(Ccsimp.lme) # nothing is significant 
Ccsimp.lme=lme(conc20~year+LU2,random=~1|site,
               data=simple20[simple20$element=='C',], na.action=na.omit)
summary(Ccsimp.lme) # now both have positive slopes and p=.06

# Are the model assumptions (normality of residuals, homogeneity of variance) met?
qqnorm(resid(Ccsimp.lme)) #ooh not at all
qqline(resid(Ccsimp.lme)) 
Ccsimp.lme=lme(log(conc20)~year+LU2,random=~1|site,
               data=simple20[simple20$element=='C',], na.action=na.omit)
summary(Ccsimp.lme) # residuals pretty normal, but nothing is significant 

Ncsimp.lme=lme(conc20~year+LU2,random=~1|site,
               data=simple20[simple20$element=='N',], na.action=na.omit)
summary(Ncsimp.lme) # definitely different between land uses (higher for O)
# no interaction so leave that off
# residuals ok, don't transform conc20
plot(resid(Ncsimp.lme)~simple20$year[simple20$element=='N' &
                                       !is.na(simple20$conc20)]) 
# variance not quite the same in different groups
plot(resid(Ncsimp.lme)~as.factor(simple20$LU2[simple20$element=='N' &
                                       !is.na(simple20$conc20)])) 
# different medians, same 25-75% range, ok.


Pcsimp.lme=lme(conc20~year+LU2,random=~1|site,
               data=simple20[simple20$element=='P',], na.action=na.omit)
summary(Pcsimp.lme) 
Kcsimp.lme=lme(log(conc20)~year*LU2,random=~1|site,
               data=simple20[simple20$element=='K',], na.action=na.omit)
summary(Kcsimp.lme) # try these with and without the interaction term
Cacsimp.lme=lme(log(conc20)~year*LU2,random=~1|site,
               data=simple20[simple20$element=='Ca',], na.action=na.omit)
summary(Cacsimp.lme) 
Cacsimp2.lme=lme(log(conc20)~year+LU2,random=~1|site,
                data=simple20[simple20$element=='Ca' & 
                                simple20$site!='JP',], na.action=na.omit)
summary(Cacsimp2.lme) # are the differences all due to JP?

Nbcsimp.lme=lme(conc20~year+LU2,random=~1|site,
               data=simple20[simple20$element=='Nb',], na.action=na.omit)
summary(Nbcsimp.lme) # signif higher in euc and in 2004
Nbcsimp.lme=lme(conc20~year*LU2,random=~1|site,
                data=simple20[simple20$element=='Nb',], na.action=na.omit)
summary(Nbcsimp.lme) # no interaction
xyplot(conc20~year|site,data=simple20[simple20$element=='C',],groups=LU2)
xyplot(conc20~year|site,data=simple20[simple20$element=='C',])
# larger differences between land uses than years, but not directional


Crat.lme=lme(stockratio~year*LU2,random=~1|site,
             data=simple20[simple20$element=='C',], na.action=na.omit)
summary(Crat.lme) 
bwplot(stockratio~LU|site+year,data=dats2deps) # look all kind of the same to me
##########

# Across all sites, does C accumulate over time?
allC20.lme=lme(stock20~year,random=~1|site/stand,
               data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20.lme) # yes, increases (marginally, adding Eu)
plot(resid(allC20.lme)~dats2deps$stock20[dats2deps$element=='C' &
                                           !is.na(dats2deps$stock20)])
# residuals still highest at highest C
plot(resid(allC20.lme)~dats2deps$year[dats2deps$element=='C' &
                                           !is.na(dats2deps$stock20)])
# more spread in 04?
qqnorm(resid(allC20.lme))
qqline(resid(allC20.lme)) # upper tail off but probably ok


allC20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
               data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20LU.lme) # yes, increases, but not in pasture
qqnorm(resid(allC20LU.lme))
qqline(resid(allC20LU.lme)) # probably ok
plot(resid(allC20LU.lme)~dats2deps$year[dats2deps$element=='C' &
                                        !is.na(dats2deps$stock20)]) #ok
plot(resid(allC20LU.lme)~dats2deps$LU[dats2deps$element=='C' &
                                          !is.na(dats2deps$stock20)]) 
# less ok, tiny spread for pasture
# which way is correct?
allC20LUsite.lme=lme(stock20~year*LU,random=~1|site/stand,#site,
                 data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20LUsite.lme) # p for increase =.07, rest insignif; .05 w/o intrxn
qqnorm(resid(allC20LUsite.lme))
qqline(resid(allC20LUsite.lme)) # upper tail off
plot(resid(allC20LUsite.lme)~dats2deps$year[dats2deps$element=='C' &
                                          !is.na(dats2deps$stock20)]) 
# now these are maybe more different
plot(resid(allC20LUsite.lme)~dats2deps$LU[dats2deps$element=='C' &
                                        !is.na(dats2deps$stock20)]) 

allP20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='P',],na.action = na.omit)
summary(allP20LU.lme) # decreases, only in native veg
allP20cratLU.lme=lme(concratio~year*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='P',],na.action = na.omit)
summary(allP20cratLU.lme) 
#no change between years, or land uses, even when each in a one-way lm

allC100LUsite.lme=lme(stock100~year*LU,random=~1|site,
                     data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC100LUsite.lme) # p for increase =.07, rest insignif; .05 w/o intrxn
qqnorm(resid(allC100LUsite.lme))
qqline(resid(allC100LUsite.lme)) # probably ok
plot(resid(allC100LUsite.lme)~dats2deps$year[dats2deps$element=='C' &
                                              !is.na(dats2deps$stock100)]) 
# prob ok
plot(resid(allC100LUsite.lme)~dats2deps$LU[dats2deps$element=='C' &
                                            !is.na(dats2deps$stock100)]) 

allC100LUstd.lme=lme(stock100~year*LU,random=~1|stand,
                      data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC100LUstd.lme) # that seems to have less power


allN100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allN100LU.lme) 
# 20 cm: increases only in natveg (signif intrxn);
#   within site not stand, both natveg and interaction signif
# 100 cm within site, only intrxn signif again
# with Eu, now increase is significant and larger in N
# should these be proportional increases?
qqnorm(resid(allN20LU.lme))
qqline(resid(allN20LU.lme)) # nice
qqnorm(resid(allN100LU.lme))
qqline(resid(allN100LU.lme)) # tails way off, one crazy outlier
plot(resid(allN20LU.lme)~dats2deps$year[dats2deps$element=='N' &
                                          !is.na(dats2deps$stock20)]) #ok
plot(resid(allN20LU.lme)~dats2deps$LU[dats2deps$element=='N' &
                                        !is.na(dats2deps$stock20)]) #eh

allP100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='P2',],na.action = na.omit)
summary(allP100LU.lme) # less in pasture, decreases in native
# p < .06
# now p=0.1 for increase between years, still decreases in native though
#   that one crazy value for JP.N

allK100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='K' &
                                   dats2deps$site!='Bp',],na.action = na.omit)
summary(allK100LU.lme) # with Bp: decreases in 2016, but increases 
# without: only signif things are in native veg (start with less, gain)
# now with Eu, increases overall, too
# residuals not as bad as some others, but still a big outlier
# lower AIC with stand within site vs just stand


allNcratLU.lme=lme(concratio~year*LU,random=~1|site/stand,
                     data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allNcratLU.lme) # concratio decreases with time
# residuals pretty normal
allNcratLUstd.lme=lme(concratio~year*LU,random=~1|stand,
                      +                      data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
AIC(allNcratLU.lme,allNcratLUstd.lme) # lower for just stand

allCacratLU.lme=lme(concratio~year*LU,random=~1|stand,
                   data=dats2deps[dats2deps$element=='Ca',],na.action = na.omit)
summary(allCacratLU.lme) # increases, nice; by more in euc, much less in nat
# more in pasture to start with 
# K: no change over time for euc, but signif year*natveg intrxn 
#   (distribution gets deeper in nat, because of JP.N weirdness)

allNcratLU2.lme=lme(concrat2~year*LU,random=~1|stand,
                   data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allNcratLU2.lme) # barely decreases, p=.054


# More t-tests, but on stocks rather than each depth
# creating a lot of columns again
tstock=group_by(dats2deps,stand,element,unit,stockunit,site,LU) %>%
  summarise(stock100_16=mean(stock100[year=='16'],na.rm=T),
            stock100_04=mean(stock100[year=='04'],na.rm=T),
            oldBDstock100_16=mean(oldBDstock100[year=='16'],na.rm=T),
            oldBDstock100_04=mean(oldBDstock100[year=='04'],na.rm=T),
            BD100_16=mean(BD100[year=='16'],na.rm=T),
            BD100_04=mean(BD100[year=='04'],na.rm=T),
            BD20_16=mean(BD20[year=='16'],na.rm=T),
            BD20_04=mean(BD20[year=='04'],na.rm=T),
            BDsd20=mean(BDsd20[year=='16'],na.rm=T),
            BDsd100=mean(BDsd100[year=='16'],na.rm=T),
            stock20_16=mean(stock20[year=='16'],na.rm=T),
            stock20_04=mean(stock20[year=='04'],na.rm=T),
            conc100_16=mean(conc100[year=='16'],na.rm=T),
            conc100_04=mean(conc100[year=='04'],na.rm=T),
            conc20_16=mean(conc20[year=='16'],na.rm=T),
            conc20_04=mean(conc20[year=='04'],na.rm=T),
            n100_16=sum(!is.na(stock100[year=='16'])),
            n100_04=sum(!is.na(stock100[year=='04'])),
            sd100_16=sd(stock100[year=='16'],na.rm=T),
            sd100_04=sd(stock100[year=='04'],na.rm=T),
            #sd100_16=mean(stocksd100[year=='16'],na.rm=T),
            #sd100_04=mean(stocksd100[year=='04'],na.rm=T),
            n20_16=sum(!is.na(stock20[year=='16'])),
            n20_04=sum(!is.na(stock20[year=='04'])),
            #sd20_16=sd(stock20[year=='16'],na.rm=T),
            #sd20_04=sd(stock20[year=='04'],na.rm=T),
            sd20_16=mean(stocksd20[year=='16'],na.rm=T),
            sd20_04=mean(stocksd20[year=='04'],na.rm=T),
            sdc100_16=sd(conc100[year=='16'],na.rm=T),#fix these later; weight by mass of soil
            sdc100_04=sd(conc100[year=='04'],na.rm=T),
            sdc20_16=sd(conc20[year=='16'],na.rm=T),
            sdc20_04=sd(conc20[year=='04'],na.rm=T),
            tstat100=ifelse((n100_16>2 & n100_04 > 2),
                         t.test(stock100[year=='16'],
                                stock100[year=='04'])$statistic,NA),
            pval100=ifelse((n100_16>2 & n100_04 > 2),
                        t.test(stock100[year=='16'],
                               stock100[year=='04'])$p.value,NA),
            tstat20=ifelse((n20_16>2 & n20_04 > 2),
                            t.test(stock20[year=='16'],
                                   stock20[year=='04'])$statistic,NA),
            pval20=ifelse((n20_16>2 & n20_04 > 2),
                           t.test(stock20[year=='16'],
                                  stock20[year=='04'])$p.value,NA),
            tstatc20=ifelse((n20_16>2 & n20_04 > 2),
                           t.test(conc20[year=='16'],
                                  conc20[year=='04'])$statistic,NA),
            pvalc20=ifelse((n20_16>2 & n20_04 > 2),
                          t.test(conc20[year=='16'],
                                 conc20[year=='04'])$p.value,NA),
            tstatc100=ifelse((n100_16>2 & n100_04 > 2),
                            t.test(conc100[year=='16'],
                                   conc100[year=='04'])$statistic,NA),
            pvalc100=ifelse((n100_16>2 & n100_04 > 2),
                           t.test(conc100[year=='16'],
                                  conc100[year=='04'])$p.value,NA),
            tstatrat=ifelse((n100_16>2 & n100_04 > 2),
                             t.test(stockratio[year=='16'],
                                    stockratio[year=='04'])$statistic,NA),
            pvalrat=ifelse((n100_16>2 & n100_04 > 2),
                            t.test(stockratio[year=='16'],
                                   stockratio[year=='04'])$p.value,NA),
            tstatcrat=ifelse((n100_16>2 & n100_04 > 2),
                            t.test(concratio[year=='16'],
                                   concratio[year=='04'])$statistic,NA),
            pvalcrat=ifelse((n100_16>2 & n100_04 > 2),
                           t.test(concratio[year=='16'],
                                  concratio[year=='04'])$p.value,NA),
            rat_16=mean(stockratio[year=='16'],na.rm=T),
            rat_04=mean(stockratio[year=='04'],na.rm=T),
            sdrat_16=sd(stockratio[year=='16'],na.rm=T),
            sdrat_04=sd(stockratio[year=='04'],na.rm=T),
            ratc_16=mean(concratio[year=='16'],na.rm=T),
            ratc_04=mean(concratio[year=='04'],na.rm=T),
            sdratc_16=sd(concratio[year=='16'],na.rm=T),
            sdratc_04=sd(concratio[year=='04'],na.rm=T))

# Test the stock variance
summary(dats2deps$stocksd100[dats2deps$element=='C'])
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
} # From Macro on StackExchange
testIt=dats2deps[dats2deps$stand=='It.E1'&dats2deps$element=='C',] 
summary(testIt$stocksd100)
t.test2(mean(testIt$stock100[testIt$year=='04'],na.rm=T),
        mean(testIt$stock100[testIt$year=='16'],na.rm=T),
        mean(testIt$stocksd100[testIt$year=='04'],na.rm=T),
        mean(testIt$stocksd100[testIt$year=='16'],na.rm=T),
        4,4,0,FALSE)

# Where did the ratio of the top 20 to the full 100 cm change?
# Don't spend too much time on this; I think it's ambiguous
#   because it mixes bulk density changes and concentration changes
#   when we mostly care about concentration. 
tstktable=function(dfr){
  table(dfr$element[dfr$pvalrat<.05],dfr$stand[dfr$pvalrat<.05])*
    table(dfr$element[dfr$tstatrat>0],dfr$stand[dfr$tstatrat>0])+
    table(dfr$element[dfr$pvalrat<.05],dfr$stand[dfr$pvalrat<.05])*
    table(dfr$element[dfr$tstatrat<0],dfr$stand[dfr$tstatrat<0])*-1
}

shorttstk=droplevels(tstock[tstock$element %in% 
                              c('C','N','P2','K','Ca2','Mg2','Fe','S','Al',
                                'Cl','Nb','Zr'),])# &
                            #  !is.element(tstock$stand,
                             #             c('Eu.N','Eu.E1','Eu.E2','JP.A')),])

#tstktable(shorttstk)
gen_ttable(shorttstk$element,shorttstk$stand,shorttstk$pvalrat,shorttstk$tstatrat,0.05)
# mostly decreases: C and N in BO.E, Bp.E1; such N in It.E1,JP.N
# P only changed in It.N (decrease = deeper)
shorttstk$rat_16[shorttstk$element=='P'&shorttstk$stand=='It.N'] # just 24%

gen_ttable(shorttstk$element,shorttstk$stand,shorttstk$pvalc100,shorttstk$tstatc100,0.05)
# average concentration: most changes in JP.N, 
#   It.N (Al, P, S, Cl, Fe--one changed and affected others?) 
#   and It.E1 (C,N, Al, P, S, Cal, Zr, Nb)
# Nb increases signif. in Bp.E1, JP.N; decreases It.E1, Vg.N
# Zr changes the same direction but not signif in Vg.N
gen_ttable(shorttstk$element,shorttstk$stand,
           shorttstk$pvalcrat,shorttstk$tstatcrat,0.05)
# most significant changes are decreases (deeper distribution in 2016)
# Ca does get shallower in JP.E1 and 2; C but not N in Vg.E
# N deeper in BO.E, Bp.E1, It.E1
# Still driven by BD changes?
gen_ttable(shorttstk$element,shorttstk$stand,
           shorttstk$pval100,shorttstk$tstat100,0.05)
# P2 does increase in It.E1 and N, Vg.N; decr JP.N and Bp.E1
# Ca2 only in It.E2, Vg.E; not in the JPs

plot(stock100_16~stock100_04,col=LU,data=shorttstk)
abline(0,1)
#stkIts=tstock[tstock$site=='It' & 
#                tstock$element %in% c('C','P','K','N','Zr','K','Ca','Al'),]
stkIts=tstock[tstock$stand %in% c('It.E1','It.N') & 
                tstock$element %in% c('C','P','K','N','Zr','K','Ca','Al'),]
stkIts$stand=factor(stkIts$stand,levels=c('It.E1','It.N'))
stkIts$yrratc20=stkIts$conc20_16/stkIts$conc20_04
stkIts$yrrats20=stkIts$stock20_16/stkIts$stock20_04
stkIts$yrchgc20=(stkIts$conc20_16-stkIts$conc20_04)*100/stkIts$conc20_04
stkIts$yrchgc100=(stkIts$conc100_16-stkIts$conc100_04)*100/stkIts$conc100_04
#plot(yrratc20~as.numeric(stand),data=stkIts, type='n',
#     xaxt='n',xlim=c(0.6,3.4),las=1,xlab='',
#     ylab='Ratio of [element] at 20 cm (2016/2004)')
par(mar=c(4.5,5,2,2))
plot(yrchgc100~as.numeric(stand),data=stkIts, type='n',
     xaxt='n',xlim=c(0.7,2.3),las=1,xlab='',cex.axis=1.3,ylab='')
     #ylab='12-year change in concentration to 20 cm, %')
text(as.numeric(stkIts$stand),stkIts$yrchgc100,
     labels=stkIts$element,cex=1.5)
axis(side=1,at=c(1,2),line=1.5,
     labels=c('Eucalyptus','Paired\nCerrado'),tick=F,cex.axis=1.3)
title(ylab= '% change in concentration, 0-100 cm',cex.lab=1.3, line=3.5)
abline(h=0, lty=3)

plot(yrrats20~as.numeric(stand),data=stkIts, type='n',
     xaxt='n',xlim=c(0.6,3.4),las=1,xlab='',
     ylab='Ratio of stock to 20 cm (2016/2004)')
text(as.numeric(stkIts$stand),stkIts$yrrats20,
     labels=stkIts$element)
axis(side=1,at=c(1,2,3),
     labels=c('Low organic\nEucalyptus',
              'High organic\nEucalyptus,\nrecently harvested',
              'High organic\nCerrado reserve'),line=2,tick=F)

yrdiffstockplot_site=function(sub_ttests){
  #palette(brewer.pal(6,'Dark2'))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  #sub_ttests=droplevels(sub_ttests)
  plot(stock20_16~stock20_04,data=sub_ttests,type='n',las=1,
       #xlab=paste(unique(sub_ttests$element),' (',
       #            unique(sub_ttests$stockunit),') ','to 20 cm, 2004',sep=''),
       #ylab=paste(unique(sub_ttests$element),' (',
       #            unique(sub_ttests$stockunit),') ','to 20 cm, 2016',sep=''),
       xlab='2004',ylab='2016',cex.lab=1.6,cex.axis=1.5,
       #xlim=c(min(sub_ttests$stock20_04-sub_ttests$sd20_04,na.rm=T)*.97,
       #      max(sub_ttests$stock20_04+sub_ttests$sd20_04,na.rm=T)*1.03),
       #ylim=c(min(sub_ttests$stock20_16-sub_ttests$sd20_16,na.rm=T)*.97,
       #      max(sub_ttests$stock20_16+sub_ttests$sd20_16,na.rm=T)*1.03)
       xlim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-20 cm',
                      sep=' '))
  segments(sub_ttests$stock20_04,sub_ttests$stock20_16-sub_ttests$sd20_16,
           sub_ttests$stock20_04,sub_ttests$stock20_16+sub_ttests$sd20_16,
           col='gray60',lwd=2)
  segments(sub_ttests$stock20_04-sub_ttests$sd20_04,sub_ttests$stock20_16,
           sub_ttests$stock20_04+sub_ttests$sd20_04,sub_ttests$stock20_16,
           col='gray60',lwd=2)
  points(stock20_16~stock20_04,data=sub_ttests,col=LU,
         #pch=as.numeric(site)*2+13,cex=2)
         pch=as.numeric(site)+14,cex=3)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
}
yrdiffstockplot100_site=function(sub_ttests){
  palette(rainbow(6))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  plot(stock100_16~stock100_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',
       xlim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
            max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
            max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05))
  abline(0,1)
  segments(sub_ttests$stock100_04,sub_ttests$stock100_16-sub_ttests$sd100_16,
           sub_ttests$stock100_04,sub_ttests$stock100_16+sub_ttests$sd100_16,col='gray60')
  segments(sub_ttests$stock100_04-sub_ttests$sd100_04,sub_ttests$stock100_16,
           sub_ttests$stock100_04+sub_ttests$sd100_04,sub_ttests$stock100_16,col='gray60')
  points(stock100_16~stock100_04,data=sub_ttests,col=site,
         pch=as.numeric(LU)+14,cex=2)
  legend('topleft',bty='n',
         legend=paste(unique(sub_ttests$element),' stock to 100 cm (',
                      unique(sub_ttests$stockunit),') ',sep=''))
  palette('default')
}

yrdiffstockplot100_LU=function(sub_ttests){
  #palette(brewer.pal(6,'Dark2'))
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  #sub_ttests=droplevels(sub_ttests)
  plot(stock100_16~stock100_04,data=sub_ttests,type='n',las=1,
       #xlab=paste(unique(sub_ttests$element),' (',
       #            unique(sub_ttests$stockunit),') ','to 20 cm, 2004',sep=''),
       #ylab=paste(unique(sub_ttests$element),' (',
       #            unique(sub_ttests$stockunit),') ','to 20 cm, 2016',sep=''),
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       #xlim=c(min(sub_ttests$stock20_04-sub_ttests$sd20_04,na.rm=T)*.97,
       #      max(sub_ttests$stock20_04+sub_ttests$sd20_04,na.rm=T)*1.03),
       #ylim=c(min(sub_ttests$stock20_16-sub_ttests$sd20_16,na.rm=T)*.97,
       #      max(sub_ttests$stock20_16+sub_ttests$sd20_16,na.rm=T)*1.03)
       xlim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-100 cm',
                      sep=' '))
#  segments(sub_ttests$stock100_04,sub_ttests$stock100_16-2*sub_ttests$sd100_16,
#           sub_ttests$stock100_04,sub_ttests$stock100_16+2*sub_ttests$sd100_16,
#           col='gray60',lwd=2)
#  segments(sub_ttests$stock100_04-2*sub_ttests$sd100_04,sub_ttests$stock100_16,
#           sub_ttests$stock100_04+2*sub_ttests$sd100_04,sub_ttests$stock100_16,
#           col='gray60',lwd=2)
# Replace std devs by std errors
    segments(sub_ttests$stock100_04,sub_ttests$stock100_16-2*sub_ttests$sd100_16/sqrt(3),
           sub_ttests$stock100_04,sub_ttests$stock100_16+2*sub_ttests$sd100_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$stock100_04-2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           sub_ttests$stock100_04+2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           col='gray60',lwd=2)
  points(stock100_16~stock100_04,data=sub_ttests,col=LU,
         #pch=as.numeric(site)*2+13,cex=2)
         #pch=as.numeric(site)+14,cex=2)#3)
         pch=as.numeric(site)+19,cex=2,bg=LU)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
  par(mar=c(4,4,2,2))
}
yrdiffstockplot_site(tstock[tstock$element=='N',])
#legend('bottomright',pch=c(16,17,15),bty='n',
#       legend=c('Eucalyptus','Native vegetation','Pasture'))
yrdiffstockplot100_LU(tstock[tstock$element=='N',])
legend('bottomright',pch=16,bty='n',#cex=1.6,
       col=c('blue3','springgreen','darkgoldenrod1'),
       legend=c('Eucalyptus','Native vegetation','Pasture'))
legend('right',pch=seq(20,25), bty='n',#cex=1.6,
       legend=c('BO','Bp','It','JP','Eu','Vg'),pt.bg=1)
yrdiffstockplot100_LU(tstock[tstock$element=='C',])
yrdiffstockplot100_LU(tstock[tstock$element=='P2',])
legend('topleft',bty='n',cex=1.8,
       legend='P (Mg / ha)\n0-100 cm')

yrdiffstockplot100_LU(tstock[tstock$element=='K'&tstock$site!='Bp',])
#yrdiffstockplot100_LU(tstock[tstock$element=='Ca2'&tstock$stock100_16<1,])

yrdiffoldBDstockplot100_LU=function(sub_ttests){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  plot(oldBDstock100_16~oldBDstock100_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16),na.rm=T)*.9,
              max(c(sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16),na.rm=T)*.9,
              max(c(sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-100 cm',
                      sep=' '))
  # Replace std devs by std errors
  segments(sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16-2*sub_ttests$sd100_16/sqrt(3),
           sub_ttests$oldBDstock100_04,sub_ttests$oldBDstock100_16+2*sub_ttests$sd100_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$oldBDstock100_04-2*sub_ttests$sd100_04/sqrt(3),sub_ttests$oldBDstock100_16,
           sub_ttests$oldBDstock100_04+2*sub_ttests$sd100_04/sqrt(3),sub_ttests$oldBDstock100_16,
           col='gray60',lwd=2)
  points(oldBDstock100_16~oldBDstock100_04,data=sub_ttests,col=LU,
         pch=as.numeric(site)+19,cex=2,bg=LU)
  palette('default')
  par(mar=c(4,4,2,2))
}

yrdiffoldBDstockplot100_LU(tstock[tstock$element=='C',])

yrdiffstockplot_site(tstock[tstock$element=='C',])
yrdiffstockplot_site(tstock[tstock$element=='K'&tstock$site!='Bp',])
yrdiffstockplot_site(tstock[tstock$element=='Ca',])
yrdiffstockplot_site(tstock[tstock$element=='Ca'&tstock$site!='JP',])
# Increased a little in all sites, but not significantly, maybe

yrdiffstockplot100_site(tstock[tstock$element=='N',])
yrdiffstockplot_site(tstock[tstock$element=='K' & 
                              tstock$stand %in% c('BO.E','BO.P','Vg.E','Vg.N',
                                                  'JP.E1','JP.N','It.E1','It.N'),])
yrdiffstockplot_site(tstock[tstock$element=='P',])

stockC=tstock[tstock$element=='C',]
# get a single subset to look at BD
palette(rainbow(6))
  stockC$LU=factor(stockC$LU,levels=c('P','E','N'))
  #stockC=droplevels(stockC)
  plot(BD20_16~BD20_04,data=stockC,type='n',las=1,
       xlab='2004',ylab='2016',
       xlim=c(min(stockC$BD20_16-stockC$BDsd20,na.rm=T)*.97,
              max(stockC$BD20_16+stockC$BDsd20,na.rm=T)*1.03),
       ylim=c(min(stockC$BD20_16-stockC$BDsd20,na.rm=T)*.97,
              max(stockC$BD20_16+stockC$BDsd20,na.rm=T)*1.03))
  abline(0,1)
  segments(stockC$BD20_04,stockC$BD20_16-stockC$BDsd20,
           stockC$BD20_04,stockC$BD20_16+stockC$BDsd20,col='gray60')
  points(BD20_16~BD20_04,data=stockC,col=site,
         pch=as.numeric(LU)+14,cex=2)
  legend('bottomright',pch=c(16,17,15),bty='n',
         legend=c('Eucalyptus','Native vegetation','Pasture'))
  legend('topleft',bty='n',
         legend=expression(paste('Bulk density, 0-20 cm (g ',cm^-3,')')))
  palette('default')

# This was a function, but no reason for that
  stockC$LU=factor(stockC$LU,levels=c('P','E','N'))
  stockC$site=factor(stockC$site,
                     levels=c('BO','Bp','It','JP','Eu','Vg'))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  par(mar=c(5,5,2,2))
  plot(BD100_16~BD100_04,data=stockC,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(stockC$BD100_16,stockC$BD100_04),na.rm=T)*.9,
              max(c(stockC$BD100_16,stockC$BD100_04),na.rm=T)*1.05),
       ylim=c(min(c(stockC$BD100_16,stockC$BD100_04),na.rm=T)*.9,
              max(c(stockC$BD100_16,stockC$BD100_04),na.rm=T)*1.05))
  abline(0,1)
  segments(stockC$BD100_04,stockC$BD100_16-2*stockC$BDsd100,
           stockC$BD100_04,stockC$BD100_16+2*stockC$BDsd100,
           col='gray60',lwd=2)
  #points(BD100_16~BD100_04,data=stockC,col=LU,
  #       pch=as.numeric(site)+14,cex=3)
  points(BD100_16~BD100_04,data=stockC,col=LU,
         pch=as.numeric(site)+19,cex=2,bg=LU)
  #legend('topleft',bty='n',cex=1.8,
  #       legend=expression(paste('Mean bulk density, 0-100 cm (g ',cm^-3,')')))
  palette('default')
  par(mar=c(4,4,2,2))
  

# plots for individual reports?
  yrdiffstockplot100_LU(tstock[tstock$element=='C' &tstock$site=='It',])
Itstock=tstock[tstock$site=='It',]
Itsum=Itstock[Itstock$element %in% c('C','N','K','P','Ca'),
              c('stand','element','stock20_04','stock20_16','stock100_04',
                'stock100_16','BD20_04','BD100_16','sd100_04','sd100_16',
                'sd20_04','sd20_16','conc100_04','conc100_16',
                'conc20_04','conc20_16','sdc20_04','sdc20_16')]

Itsum100=Itstock[Itstock$element %in% c('C','N','K','P','Ca2'),
                 c('stand','element','stock100_04','stock100_16',
                   'sd100_04','sd100_16','tstat100','pval100')]
Itsum100
# Is that right? stocks increase significantly for almost everything?
# Error in the sd calculation?
Itsum100=Itstock[Itstock$element %in% c('C','N','K','P','Ca2'),
                 c('stand','element','stock100_04','stock100_16',
                   'sd100_04','sd100_16','BD100_04','BD100_16','BDsd100')]

# Test pit vs not
# temporary fix: 
widedats$ID=as.character(widedats$ID)
widedats$ID[widedats$ID=='It.N.T.A.60-100.16'|
              widedats$ID=='It.N.TA.60-100.16']='It.N.A.5.60-100.16'
widedats$ID[widedats$ID=='It.N.T.B.60-100.16'|
              widedats$ID=='It.N.TB.60-100.16']='It.N.B.5.60-100.16'
widedats$ID[widedats$ID=='Vg.N.1.B.60-100.16']='Vg.N.B.1.60-100.16'
widedats$ID[widedats$ID=='Vg.N.2.B.60-100.16']='Vg.N.B.2.60-100.16'
widedats$ID[widedats$ID=='Vg.N.1.B.0-10.16']='Vg.N.B.1.0-10.16'
widedats$ID[widedats$ID=='Vg.N.2.B.0-10.16']='Vg.N.B.2.0-10.16'
widedats=strfun(widedats)

widedats=group_by(widedats,stand,year,depth) %>%
  mutate(maxrep=max(rep,na.rm=T))
pits=widedats[widedats$maxrep==5,]
plot(C~rep,data=pits,col=stand,pch=as.numeric(as.factor(depth)))
plot(C~rep,data=pits[pits$depth==50,],
     col=stand,pch=as.numeric(year))
# a bunch of values are missing, ugh
# for 40-60, pit samples have much less C in It.E1, also less in It.N
plot(P~rep,data=pits[pits$depth==50,],
     col=stand,pch=as.numeric(year))

# The plots to check which samples are completed
CNs=widedats[!is.na(widedats$C),]
CNs=group_by(CNs,stand,depth,year)%>% 
  mutate(nreps=n())
CNs$nreps=as.character(CNs$nreps)
CNs=ungroup(CNs)
plot(as.numeric(stand)~depth,data=CNs,
     yaxt='n',ylab='',type='n',xlim=c(5,90),main='CN analyzed')
text(CNs$depth[CNs$year=='04'],
     as.numeric(CNs$stand[CNs$year=='04']),CNs$nreps[CNs$year=='04'])
axis(side=2, at=seq(length(levels(CNs$stand))),
     labels=levels(CNs$stand),las=1)
text(CNs$depth[CNs$year=='16']+3,
     as.numeric(CNs$stand[CNs$year=='16']),CNs$nreps[CNs$year=='16'],col=2)
# still a couple of typos

XRF=widedats[!is.na(widedats$K),]
XRF=group_by(XRF,stand,depth,year)%>% 
  mutate(nreps=n())
XRF$nreps=as.character(XRF$nreps)
XRF=ungroup(XRF)
plot(as.numeric(stand)~depth,data=XRF,
     yaxt='n',ylab='',type='n',xlim=c(5,90),main='XRF analyzed')
text(XRF$depth[XRF$year=='04'],
     as.numeric(XRF$stand[XRF$year=='04']),XRF$nreps[XRF$year=='04'])
axis(side=2, at=seq(length(levels(XRF$stand))),
     labels=levels(XRF$stand),las=1)
text(XRF$depth[XRF$year=='16']+3,
     as.numeric(XRF$stand[XRF$year=='16']),XRF$nreps[XRF$year=='16'],col=2)

boths=widedats[!is.na(widedats$K) & !is.na(widedats$C),]
boths=group_by(boths,stand,depth,year)%>% 
  mutate(nreps=n())
boths$nreps=as.character(boths$nreps)
boths=ungroup(boths)
plot(as.numeric(stand)~depth,data=boths,
     yaxt='n',ylab='',type='n',xlim=c(5,90),main='Both analyzed')
text(boths$depth[boths$year=='04'],
     as.numeric(boths$stand[boths$year=='04']),boths$nreps[boths$year=='04'])
axis(side=2, at=seq(length(levels(boths$stand))),
     labels=levels(boths$stand),las=1)
text(boths$depth[boths$year=='16']+3,
     as.numeric(boths$stand[boths$year=='16']),boths$nreps[boths$year=='16'],col=2)

