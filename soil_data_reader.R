# set your working directory (or open a project containing all the files)
#setwd('C:\\Users\\Devin\\Documents\\Soil data')

# useful packages:
library(dplyr)
library(lattice)
library(nlme)
library(reshape)
library(RColorBrewer)
library(MASS)
library(ggplot2)
source('xrf_formatting_functions.R')

gen_ttable=function(elmtcol,groupcol,pvalcol,tstatcol,a){
  table(elmtcol[pvalcol<a],groupcol[pvalcol<a])*
    table(elmtcol[tstatcol>0],groupcol[tstatcol>0])+
    table(elmtcol[pvalcol<a],groupcol[pvalcol<a])*
    table(elmtcol[tstatcol<0],groupcol[tstatcol<0])*-1
}

qqr=function(lmeobject){
  qqnorm(resid(lmeobject),main=lmeobject$call,cex.main=.8)
  qqline(resid(lmeobject))
}


#widedats=readRDS('all_data_9-13-18.Rds')
#widedats=readRDS('all_data_9-29-18.Rds')
widedats=readRDS('all_data_10-4-18.Rds')
#widedats=readRDS('all_data_10-10-18.Rds') # separate N correction for JP.E2
# Makes a really big N stock change
# Go with consistent conversion for all stands (all_data_10-4)
# temporary fix: 
#widedats$ID=as.character(widedats$ID)
#widedats$ID[widedats$ID=='It.N.T.A.60-100.16'|
#              widedats$ID=='It.N.TA.60-100.16']='It.N.A.5.60-100.16'
#widedats$ID[widedats$ID=='It.N.T.B.60-100.16'|
#              widedats$ID=='It.N.TB.60-100.16']='It.N.B.5.60-100.16'
#widedats$ID[widedats$ID=='Vg.N.1.B.60-100.16']='Vg.N.B.1.60-100.16'
#widedats$ID[widedats$ID=='Vg.N.2.B.60-100.16']='Vg.N.B.2.60-100.16'
#widedats$ID[widedats$ID=='Vg.N.1.B.0-10.16']='Vg.N.B.1.0-10.16'
#widedats$ID[widedats$ID=='Vg.N.2.B.0-10.16']='Vg.N.B.2.0-10.16'
#widedats=strfun(widedats)

widedats$CN=widedats$C/widedats$N # this ratio can be useful, or just noisy
widedats1=widedats[,-which(names(widedats) %in% 
                             c('stand','depth','year','site',
                               'LU','position','rep','elt',
                               'inc_from','inc_to','newest','oldmeas'))]#,
                               #'P_ppm','Fe_pct','Ca_ppm',
                               #'K_ppm','Mg_ppm'))]

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
# Add a biome column (Atlantic Forest or Cerrado)
dats=mutate(dats,biome=ifelse(site %in% c('BO','Vg','Eu'),'AF','Cer'))
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

# Take out crazy outliers
dats$value[dats$ID=='JP.E2.T.2.0-10.16'& 
             dats$element %in% c('Ca','Mg','Ca2','Mg2')]=NA 
# piece of dolomite got in here, maybe
dats$value[dats$ID=='Eu.E1.L.2.0-10.16'& 
             dats$element %in% c('Ca','Ca2')]=NA 
# also one point that really changes the whole analysis
#   probably real but not representative (chunk of fertilizer, e.g.)

dats$value[dats$element=='P' & dats$value>600 & 
             dats$stand!='Eu.E1']=NA
dats$value[dats$element=='P2' & dats$value>600 & 
             dats$stand!='Eu.E1']=NA
# a little more subjective-- three points >> others in their stand 

# Pit samples from Vg.E.04 were called rep 5 for XRF, rep 1 for C&N
#  similar to other samples, probably makes no difference?
#   Except that C is higher for 10-20 and 20-40 depths
# Rep 1 was missing so pit used instead to get 4

# Remove duplicates if just a few simple samples were analyzed 
#   in addition to composite samples
# group by rep; if there exist positions and no-positions for same rep and year
#   just take the no-positions
datspre=group_by(dats,stand,year,depth,element,rep) %>%
  mutate(napos=sum(is.na(position[!is.na(value)])),
         rmsimp=ifelse(napos>0 & !is.na(position),1,0))
datspre=ungroup(datspre)
# Average the row positions to get 4 reps per site (5 if there's a soil pit measurement)
# This way, we aren't weighting the individual cores and the composited cores the same
dats4=group_by(datspre[datspre$rmsimp==0,],stand,year,depth,rep,element) %>%
#dats4=group_by(dats,stand,year,depth,rep,element) %>%
  mutate(repval=mean(value,na.rm=T),repsd=sd(value,na.rm=T),
         repn=sum(!is.na(value)),repwt=ifelse(rep==5,2,1))
dats4=distinct(dats4,stand,year,depth,rep,element,.keep_all = T)
# In future stats, maybe weight pit samples more than core samples when both taken
# Pit samples were taken when cores dubious (i.e. contaminated with surface OM)

# An analogous version with the "wide" data, mostly useful if you want to compare
#   different elements within the same sample
widedats4=group_by(widedats,stand,year,depth,rep,site,LU) %>%
  summarise_if(is.numeric, mean,na.rm=T)

# Calculate pseudo-tau relative to Zr (or other element) in 2004
dats4=group_by(dats4,stand,rep,depth) %>%
  mutate(Zrrat=ifelse(sum(!is.na(repval[element=='Zr'& year=='04']))>0 &
                        sum(!is.na(repval[element=='Zr'& year=='16']))>0,
                      repval[element=='Zr'&year=='16']/
                        repval[element=='Zr'&year=='04'],NA))
dats4=group_by(dats4,stand,rep,depth,element)%>%
  mutate(tau=ifelse(!is.na(Zrrat),
                    ((repval[year=='16']/repval[year=='04'])/Zrrat)-1,
                    NA))


##### Further analysis of different depths
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
shorttests=droplevels(ttests[ttests$element %in% 
                               c('C','N','P','K','Ca','Mg'),])

# For ease of visualization, get a single average value per stand, depth, 
#   year, and element
datsmean=group_by(dats4[!is.element(dats4$unit,c('dl','err','int')),],
                  stand,depth,year,element,unit,stdonly,site,LU,biome,avgBD,BDsd) %>%
  summarise(mn=mean(repval,na.rm=T),sd=sd(repval,na.rm=T),
            taumn=mean(tau,na.rm=T),sdtau=sd(tau,na.rm=T))
datsmnok=droplevels(datsmean[datsmean$site!='TM' &
                               datsmean$site!='Cr'& datsmean$LU!='A',])

# Use the lattice functions to visualize changes at different depths
# Compare to the t-test tables to see what differences between years 
#   are statistically significant
datsmnok=datsmnok[order(datsmnok$depth),]

abdatsmn=droplevels(datsmean[is.element(datsmean$site, c('TM','Cr','JP'))&
                               !is.element(datsmean$stand,c('JP.P'))&
                               datsmean$year==16,])

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
                 site,LU,biome,avgBD,oldBD,newBDsd) %>%
  mutate(inc=inc_to-inc_from,  repval=mean(value,na.rm=T),repvar=var(value,na.rm=T),
         #stock=ifelse(unit=='pct',repval*avgBD*inc,repval*avgBD*inc*.0001),
         # Alternative: use 2016 values for everything
         stock=ifelse(unit=='pct',repval*BD16*inc,repval*BD16*inc*.0001),
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

# Added 9-25-18
datsstk5=group_by(datsstk,stand,year,element) %>%
  mutate(maxrep=max(rep[inc_to>=60],na.rm=T),
         pitstock=ifelse(maxrep==5,
                         stock[rep==5 & inc_to==60]+stock[rep==5 & inc_to==100],
                         NA))
dats2deps=group_by(datsstk5[!is.element(datsstk5$site,c('TM','Cr'))&
                             datsstk5$stand!='JP.A',],stand,year,rep,
                   element,site,LU,biome,unit,stockunit) %>% 
  summarise(ndeps20=length(unique(inc_to[inc_to<=20])),
            stock20=ifelse(ndeps20==2,sum(stock[inc_to<=20]),NA),
            stocklo20=ifelse(ndeps20==2,sum(stocklo[inc_to<=20]),NA),
            stockhi20=ifelse(ndeps20==2,sum(stockhi[inc_to<=20]),NA),
            oldBDstock20=ifelse(ndeps20==2,sum(oldBDstock[inc_to<=20]),NA),
            #conc20=ifelse(ndeps20==2,
            #              sum(repval[inc_to<=20]*inc[inc_to<=20]*avgBD[inc_to<=20])/
            #                sum(inc[inc_to<=20]*avgBD[inc_to<=20]),NA),
            conc20=ifelse(ndeps20==2,
                          sum(repval[inc_to<=20]*inc[inc_to<=20]*BD16[inc_to<=20])/
                            sum(inc[inc_to<=20]*BD16[inc_to<=20]),NA),
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
                            #ifelse(maxrep==5 & year=='16',
                            #       sum(stock[inc_to<=40])+sum(stock),NA)),
            
            stocklo100=ifelse(ndeps100==5,sum(stocklo),NA),
            stockhi100=ifelse(ndeps100==5,sum(stockhi),NA),
            oldBDstock100=ifelse(ndeps100==5,sum(oldBDstock),NA),
            #conc100=ifelse(ndeps100==5,
            #               sum(repval*inc*avgBD)/sum(inc*avgBD),NA),
            conc100=ifelse(ndeps100==5,
                          sum(repval*inc*BD16)/sum(inc*BD16),NA),
            conc60to100=mean(repval[inc_to==100],na.rm=T),
            BD100=ifelse(ndeps100==5,sum(avgBD*inc)/sum(inc),NA),
            #BDsd100=ifelse(ndeps100==5,sqrt(sum(BDsd^2)),NA),
            # weight by inc
            suminc=sum(inc),
            BDsd100=ifelse(ndeps100==5,sqrt(sum((BDsd^2)*inc)/sum(inc)),NA),
            stocksd100=ifelse(ndeps100==5,sqrt(sum(stockvar*inc)/sum(inc)),NA),
            stockratio=stock20/stock100,concratio=conc20/conc100,
            concrat2=conc20/conc60to100)


euc2deps=droplevels(dats2deps[dats2deps$LU=='E',])
test2deps=dats2deps[-which(dats2deps$stand %in% c('It.E1','It.N')),]
euc2deps2=droplevels(test2deps[test2deps$LU=='E',]) 

# More t-tests, but on stocks rather than each depth
# creating a lot of columns again
tstock=group_by(dats2deps,stand,element,unit,stockunit,site,LU,biome) %>%
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
            sd20_16=sd(stock20[year=='16'],na.rm=T),
            sd20_04=sd(stock20[year=='04'],na.rm=T),
            # within-profile variance not relevant here
            #sd20_16=mean(stocksd20[year=='16'],na.rm=T),
            #sd20_04=mean(stocksd20[year=='04'],na.rm=T),
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

shorttstk=droplevels(tstock[tstock$element %in% 
                              c('C','N','P2','K','Ca2','Mg2','Fe','S','Al',
                                'Cl','Nb','Zr'),])# &
#  !is.element(tstock$stand,
#             c('Eu.N','Eu.E1','Eu.E2','JP.A')),])

shorttstk$element=as.character(shorttstk$element)
shorttstk$element[shorttstk$element=='P2']='P'
shorttstk$element[shorttstk$element=='Ca2']='Ca'
shorttstk$element[shorttstk$element=='Mg2']='Mg'

#budgets=read.csv('nutrient_budget_summary.csv')
budgets=read.csv('nutrient_budgets_linked_bark.csv')
otherconcs=data.frame(Egrandconc=c(0.00118814,0.000030795,0.00037971,0.001193941,
                               0.000128281,0.00000841287,0.000077024, NA, NA),
                      Plconc=c(0.000599674, 0.0000700706,0.000845645,
                               0.001126018,0.000192921,0.00000573683,
                               0.000256473,0.00000234058,0.00000458012),
                      SantanaMG=c(0.001365204, 0.000110556, 0.000899214,
                                0.001577597, 0.000291741, NA, NA, NA, NA),
                  Nutrient=c('N','P','K','Ca','Mg','B','S','Cu','Zn'))
# E. grandis from "Pagano 2013", 
# Pl from Plantar data (NUTREEcalc 2015 for Itacambira)
# Missing bark values (BO and Eu) filled with Plantar bark values
budgets=merge(budgets,otherconcs,by='Nutrient')

budgets=mutate(budgets,
               grandconcbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                              Egrandconc*511)/1000,
               modconcbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                                SantanaMG*511)/1000,
               plconcbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                                Plconc*511)/1000,
               denserbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                             Concentration*560)/1000,
               lessdensebudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                                Concentration*460)/1000,
               lessrotbudg=(In_kgha_1+In_kgha_2*Past_harvests_known-
                              (Wood_m3_1+Wood_m3_2*Past_harvests_known)*
                              Concentration*511)/1000,
               lessharvbudg=(In_kgha_1+In_kgha_2*Past_harvests_known-
                               (Wood_m3_1*.5+Wood_m3_2*Past_harvests_known)*
                               Concentration*511)/1000,
               moreharvbudg=(In_kgha_1+In_kgha_2*Past_harvests_known-
                               (Wood_m3_1*1.5+Wood_m3_2*Past_harvests_known)*
                               Concentration*511)/1000,
               woodonlybudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                             Conc_wood*511)/1000,
               bark20budg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                               (Conc_wood*.8+Conc_bark*.2)*511)/1000,
               bark5budg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                            (Conc_wood*.95+Conc_bark*.05)*511)/1000,
               # Pick more probable values based on plots
               budget=ifelse(Stand=='Bp.E1',lessharvbudg,
                             ifelse(Stand=='BO.E',bark5budg,
                                    (In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                         Concentration*511)/1000))
)
shorterstk=merge(shorttstk,budgets,by.x=c('stand','element'),
                 by.y=c('Stand','Nutrient'))

stkchgs=group_by(droplevels(shorterstk),stand,element,biome)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            sdchg20=sqrt(sd20_04^2+sd20_16^2)/2, 
            # = sqrt(sd04^2/n04+sd16^2/n16), assuming n=4
            sdchg100=sqrt(sd100_04^2+sd100_16^2)/2,
            stk20_04=stock20_04,efs20=log((stk20_04+budget)/stk20_04),
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04),
            #budget=Budget/1000,
            budget=budget,
            grandconcbudg=grandconcbudg, plconcbudg=plconcbudg,
            denserbudg=denserbudg, lessrotbudg=lessrotbudg,
            lessdensebudg=lessdensebudg,conc=Concentration,
            lessharvbudg=lessharvbudg,moreharvbudg=moreharvbudg,
            woodonlybudg=woodonlybudg,bark5budg=bark5budg,
            bark20budg=bark20budg,modconcbudg=modconcbudg,
            minbudg=min(grandconcbudg,plconcbudg,denserbudg,lessdensebudg,
                        modconcbudg,lessrotbudg,lessharvbudg,moreharvbudg,
                        woodonlybudg,bark5budg,bark20budg,budget,na.rm=T),
            maxbudg=max(grandconcbudg,plconcbudg,denserbudg,lessdensebudg,
                        modconcbudg,lessrotbudg,lessharvbudg,moreharvbudg,
                        woodonlybudg,bark5budg,bark20budg,budget,na.rm=T),
            minbudgconc=min(grandconcbudg,plconcbudg,denserbudg,lessdensebudg,
                        modconcbudg,#lessrotbudg,lessharvbudg,moreharvbudg,
                        woodonlybudg,bark5budg,bark20budg,budget,na.rm=T),
            maxbudgconc=max(grandconcbudg,plconcbudg,denserbudg,lessdensebudg,
                        modconcbudg,#lessrotbudg,lessharvbudg,moreharvbudg,
                        woodonlybudg,bark5budg,bark20budg,budget,na.rm=T))#,
           # whichminconc=names(stkchgs)[which.min(c(grandconcbudg,plconcbudg,denserbudg,
           #                                lessdensebudg,modconcbudg,woodonlybudg,
           #                                bark5budg,bark20budg,budget))])
# I want to know which column is the min, but that doesn't work
stkchgs2=stkchgs[stkchgs$stand!='It.E1',]

yrdiffstockplot100_LU=function(sub_ttests){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  plot(stock100_16~stock100_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-100 cm',
                      sep=' '))
  # Replace std devs by std errors
  segments(sub_ttests$stock100_04,sub_ttests$stock100_16-2*sub_ttests$sd100_16/sqrt(3),
           sub_ttests$stock100_04,sub_ttests$stock100_16+2*sub_ttests$sd100_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$stock100_04-2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           sub_ttests$stock100_04+2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           col='gray60',lwd=2)
  points(stock100_16~stock100_04,data=sub_ttests,col=LU,
         pch=as.numeric(site)+19,cex=2,bg=LU)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
  par(mar=c(4,4,2,2))
}

yrdiffstockplot20_LU=function(sub_ttests){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  plot(stock20_16~stock20_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-20 cm',
                      sep=' '))
  # Replace std devs by std errors
  segments(sub_ttests$stock20_04,sub_ttests$stock20_16-2*sub_ttests$sd20_16/sqrt(3),
           sub_ttests$stock20_04,sub_ttests$stock20_16+2*sub_ttests$sd20_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$stock20_04-2*sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           sub_ttests$stock20_04+2*sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           col='gray60',lwd=2)
  points(stock20_16~stock20_04,data=sub_ttests,col=LU,
         pch=as.numeric(site)+19,cex=2,bg=LU)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
  par(mar=c(4,4,2,2))
}


yrdiffratplot_LU=function(sub_ttests){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  plot(rat_16~rat_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$rat_04,sub_ttests$rat_16),na.rm=T)*.9,
              max(c(sub_ttests$rat_04,sub_ttests$rat_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$rat_04,sub_ttests$rat_16),na.rm=T)*.9,
              max(c(sub_ttests$rat_04,sub_ttests$rat_16),na.rm=T)*1.05))
  abline(0,1)
  abline(h=0.2,lty=3)
  abline(v=0.2,lty=3)
  legend('topleft',bty='n',cex=1.5,
         legend=unique(sub_ttests$element))
  # Replace std devs by std errors
  segments(sub_ttests$rat_04,sub_ttests$rat_16-sub_ttests$sdrat_16/sqrt(3),
           sub_ttests$rat_04,sub_ttests$rat_16+sub_ttests$sdrat_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$rat_04-sub_ttests$sdrat_04/sqrt(3),sub_ttests$rat_16,
           sub_ttests$rat_04+sub_ttests$sdrat_04/sqrt(3),sub_ttests$rat_16,
           col='gray60',lwd=2)
  points(rat_16~rat_04,data=sub_ttests,col=LU,
         pch=as.numeric(site)+19,cex=2,bg=LU)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
  par(mar=c(4,4,2,2))
}

# Needs work:
yrdiffstockplot20_bmLU=function(sub_ttests){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$biome=factor(sub_ttests$biome,levels=c('AF','Cer'))
  plot(bmstock20_16~bmstock20_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$bmstock20_04,sub_ttests$bmstock20_16),na.rm=T)*.9,
              max(c(sub_ttests$bmstock20_04,sub_ttests$bmstock20_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$bmstock20_04,sub_ttests$bmstock20_16),na.rm=T)*.9,
              max(c(sub_ttests$bmstock20_04,sub_ttests$bmstock20_16),na.rm=T)*1.05))
  abline(0,1)
  legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-20 cm',
                      sep=' '))
  # Replace std devs by std errors
  segments(sub_ttests$bmstock20_04,sub_ttests$bmstock20_16-sub_ttests$sd20_16/sqrt(sub_ttests$ngrp-1),
           sub_ttests$bmstock20_04,sub_ttests$bmstock20_16+sub_ttests$sd20_16/sqrt(sub_ttests$ngrp-1),
           col='gray60',lwd=2)
  segments(sub_ttests$bmstock20_04-sub_ttests$sd20_04/sqrt(sub_ttests$ngrp-1),sub_ttests$bmstock20_16,
           sub_ttests$bmstock20_04+sub_ttests$sd20_04/sqrt(sub_ttests$ngrp-1),sub_ttests$bmstock20_16,
           col='gray60',lwd=2)
  points(bmstock20_16~bmstock20_04,data=sub_ttests,col=LU,
         pch=as.numeric(biome)+14,cex=2,bg=LU)
  #legend('bottomright',pch=c(16,17,15),bty='n',
  #       legend=c('Eucalyptus','Native vegetation','Pasture'))
  palette('default')
  par(mar=c(4,4,2,2))
}

yrdiffstockplot20_bmLUall=function(sub_ttests,label=T,fulllegend=F){
  #par(mar=c(5,5,2,2))
  par(mar=c(4,4,.5,.5))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$biome=factor(sub_ttests$biome,levels=c('Cer','AF'))
  plot(stock20_16~stock20_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05))
  abline(0,1)
  if(label==T){legend('topleft',bty='n',#cex=1.8,
         legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-20 cm',
                      sep=' '))}
  # Replace std devs by std errors
  segments(sub_ttests$stock20_04,sub_ttests$stock20_16-2*sub_ttests$sd20_16/sqrt(3),
           sub_ttests$stock20_04,sub_ttests$stock20_16+2*sub_ttests$sd20_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$stock20_04-2*sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           sub_ttests$stock20_04+2*sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           col='gray60',lwd=2)
  points(stock20_16~stock20_04,data=sub_ttests,col=LU,
         pch=as.numeric(biome)+15,cex=2,bg=LU)
  if(fulllegend==T){legend('bottomright',pch=c(15,15,15,17,16),bty='n',#cex=1.6,
         col=c('blue3','springgreen','darkgoldenrod1','gray50','gray50'),
         legend=c('Eucalyptus','Native vegetation','Pasture',
                  'Atlantic Forest','Cerrado'))}
  palette('default')
  par(mar=c(4,4,2,2))
}

yrdiffstockplot100_bmLUall=function(sub_ttests,label=T,fulllegend=F){
  par(mar=c(5,5,2,2))
  palette(c('darkgoldenrod1','blue3','springgreen'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$biome=factor(sub_ttests$biome,levels=c('Cer','AF'))
  plot(stock100_16~stock100_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*.9,
              max(c(sub_ttests$stock100_04,sub_ttests$stock100_16),na.rm=T)*1.05))
  abline(0,1)
  if(label==T){legend('topleft',bty='n',#cex=1.8,
                      legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-100 cm',
                                   sep=' '))}
  # Replace std devs by std errors
  segments(sub_ttests$stock100_04,sub_ttests$stock100_16-2*sub_ttests$sd100_16/sqrt(3),
           sub_ttests$stock100_04,sub_ttests$stock100_16+2*sub_ttests$sd100_16/sqrt(3),
           col='gray60',lwd=2)
  segments(sub_ttests$stock100_04-2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           sub_ttests$stock100_04+2*sub_ttests$sd100_04/sqrt(3),sub_ttests$stock100_16,
           col='gray60',lwd=2)
  points(stock100_16~stock100_04,data=sub_ttests,col=LU,
         pch=as.numeric(biome)+15,cex=2,bg=LU)
  if(fulllegend==T){legend('bottomright',pch=c(15,15,15,17,16),bty='n',#cex=1.6,
                           col=c('blue3','springgreen','darkgoldenrod1','gray50','gray50'),
                           legend=c('Eucalyptus','Native vegetation','Pasture',
                                    'Atlantic Forest','Cerrado'))}
  palette('default')
  par(mar=c(4,4,2,2))
}


simple20=droplevels(dats2deps[dats2deps$stand %in% 
                                c('BO.E','BO.P','Vg.E','Vg.N','Eu.E2','Eu.N',
                                  'JP.E2','JP.N','It.E1','It.N'),])
simple20=mutate(simple20,LU2=ifelse(LU=='E','E','O'))
# for this analysis, O = "other"

simple20_2=dats2deps
simple20_2$site=as.character(simple20_2$site)
simple20_2=mutate(simple20_2,site2=ifelse(stand=='JP.E2'|stand=='JP.N','JP2',site),
                  LU2=ifelse(LU=='E','E','O'))
simple20_2=droplevels(simple20_2[simple20_2$stand %in% 
                                   c('BO.E','BO.P','Vg.E','Vg.N','Eu.E2','Eu.N',
                                     'JP.E1','JP.N','JP.E2','JP.P','It.E1','It.N'),])
simple20_2=mutate(simple20_2,LUlong=ifelse(LU=='E','Eucalyptus',
                                           ifelse(LU=='P','Pasture',
                                                  'Native vegetation')),
                  biomelong=ifelse(biome=='Cer','Cerrado','Atlantic Forest'),
                  LU2long=ifelse(LU2=='E','Eucalyptus','Other vegetation'))
simple20_2$LUlong=factor(simple20_2$LUlong,
                         levels=c('Eucalyptus','Native vegetation','Pasture'))
simple20_2$LU2long=factor(simple20_2$LU2long,
                          levels=c('Eucalyptus','Other vegetation'))
simple20_2$biomelong=factor(simple20_2$biomelong,
                            levels=c('Atlantic Forest','Cerrado'))

simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]
shorttstk$site=as.character(shorttstk$site)
shorttstk=mutate(shorttstk,
                 site2=ifelse(stand=='JP.E2'|stand=='JP.N','JP2',site),
                 LU2=ifelse(LU=='E','E','O'))
mrstkchgs=group_by(droplevels(shorttstk[shorttstk$element!='Cl',]),
                   stand,site,LU,site2,LU2,element,biome)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04),
            sig20=ifelse(pval20<0.05,2,1),
            sig100=ifelse(pval100<0.05,2,1))
mrstkchgs$element=factor(mrstkchgs$element,levels=
                           c('C','N','K','P','S','Ca','Mg',
                             'Al','Fe','Nb','Zr'))
mrstkchgs$LU=factor(mrstkchgs$LU,levels=c('E','N','P'))
mrstkchgs=mutate(mrstkchgs,LUlong=ifelse(LU=='E','Eucalyptus',
                                         ifelse(LU=='P','Pasture',
                                                'Native vegetation')),
                 biomelong=ifelse(biome=='Cer','Cerrado','Atlantic Forest'),
                 LU2long=ifelse(LU2=='E','Eucalyptus','Other vegetation'),
                 LUlongish=ifelse(LU=='E','Euc',ifelse(LU=='P','Past','Nat')))
mrstkchgs$LUlong=factor(mrstkchgs$LUlong,
                        levels=c('Eucalyptus','Native vegetation','Pasture'))
mrstkchgs$LUlongish=factor(mrstkchgs$LUlongish,
                           levels=c('Euc','Nat','Past'))
mrstkchgs$LU2long=factor(mrstkchgs$LU2long,
                         levels=c('Eucalyptus','Other vegetation'))
mrstkchgs$biomelong=factor(mrstkchgs$biomelong,
                           levels=c('Atlantic Forest','Cerrado'))
# in longer script, tried paired t-tests between euc and noneuc
#   comparing (log) changes in stocks, n=6 sites
# Decided to use lmes instead


L_to_pct=function(x){
  (exp(x)-1)*100
}
pct_to_L=function(l){
  log((l/100)+1)
}
panel.mean <- function(x, y, ...) {
  tmp <- tapply(y, x, FUN = mean)#; print(tmp)
  panel.points(y=tmp, x=seq_along(tmp), ...)
}
# to reflect simp.lme3 analysis as in text:

# Set up data frame of change by stand, for visualization of 
#   simp.lme3 analyses (i.e. change in stock)
# Is this misleading? Didn't actually analyze change in stock
#   as a response variable
mrstkchgslim=mrstkchgs[mrstkchgs$element %in% c('C','N','K','P',
                                                'Ca'),]
#mrstkchgslim=mrstkchgs[mrstkchgs$element %in% c('C','N','K','P',
#                                               'Ca','Al','Nb'),]
#mrmelt=melt.data.frame(mrstkchgslim,measure.vars=c("chgln20","chgln100"))
# "Names do not match previous names", still
mrd=mutate(mrstkchgslim,depth=rep('0-100'))
mrd=mrd[,-which(names(mrd) %in% c('stk20_16', 'chg20', 'chgln20',
                                  'chgrt20', 'sig20'))]
mru=mutate(mrstkchgslim,depth=rep('0-20'))
mru=mru[,-which(names(mru) %in% c('stk100_16', 'chg100', 'chgln100',
                                  'chgrt100', 'sig100'))]
names(mrd)[which(names(mrd)=='stk100_16')]='stock_16'
names(mrd)[which(names(mrd)=='chg100')]='change'
names(mrd)[which(names(mrd)=='chgln100')]='chgln'
names(mrd)[which(names(mrd)=='chgrt100')]='chgrt'
names(mrd)[which(names(mrd)=='sig100')]='sig'
names(mru)[which(names(mru)=='stk20_16')]='stock_16'
names(mru)[which(names(mru)=='chg20')]='change'
names(mru)[which(names(mru)=='chgln20')]='chgln'
names(mru)[which(names(mru)=='chgrt20')]='chgrt'
names(mru)[which(names(mru)=='sig20')]='sig'
mr2=rbind(mru,mrd)

mr2=mr2[mr2$stand %in% c('BO.E','BO.P','Vg.E','Vg.N','Eu.E2','Eu.N',
                         'JP.E1','JP.N','JP.E2','JP.P','It.E1','It.N'),]
#mr2$chgln[mr2$biome=='Cer'&mr2$LU=='N'&mr2$depth=='0-100']=NA
# remove cerrado natveg? this messes up the plot
# only if varwidth=T: can't handle different numbers of obs 
#   in paired boxes
# remove pairs together?
mr2$chgln[mr2$site2=='JP2'&mr2$depth=='0-100']=NA
mr2$chgln[mr2$site2=='It'&mr2$depth=='0-100']=NA
mr2=group_by(mr2,element,LU,depth) %>%
  mutate(chglnmn = mean(chgln,na.rm=T),newnobs=n())
mr2=group_by(mr2,element) %>%
  mutate(maxchgln=max(chgln,na.rm=T),minchgln=min(chgln,na.rm=T))

mr2=mutate(mr2,sigyr=rep(NA),sigveg=rep(NA))
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU=='E']='+'
# Not with different slopes for different groups
# Also, lme doesn't say if changes significant within N and P
# Just focus on if euc changed and if N and P are diff from euc
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU=='P']='-'
mr2$sigyr[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='E']='+'
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='P']='+'
mr2$sigveg[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='N']='*'
mr2$sigveg[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU!='E']='*'
#mr2$sigyr[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU!='E']='-'
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='N' & mr2$LU=='E']='+'
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='K' & mr2$LU=='N']='*'
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='E']='+'
#mr2$sigveg[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='N']='*' nope
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU=='P']='*'

mr2$depth=factor(mr2$depth,levels=c('0-20','0-100'))

simple20_2$LU=factor(simple20_2$LU,levels=c('E','N','P'))
simple20_2lim=simple20_2[simple20_2$element %in% c('C','N','K','P','Ca'),]


simple20_2lim=simple20_2[simple20_2$element %in% c('C','N','K','P','Ca'),]
# Make a new column to show which euc stands are paired with native vs pasture
simple20_2lim=mutate(simple20_2lim,LU3=as.character(LU))
simple20_2lim$LU3[simple20_2lim$stand=='JP.E1']='EP'
simple20_2lim$LU3[simple20_2lim$LU3=='E']='EN'
simple20_2lim$LU3=factor(simple20_2lim$LU3,levels=c('EN','N','EP','P'))

# Simpler version with less data?
simple20_2limd=group_by(simple20_2lim,element,year,LU3) %>%
  mutate(mn20=mean(stock20,na.rm=T),nobs=n(),
         se20=sd(stock20,na.rm=T)/sqrt(nobs),
         min20=min(stock20,na.rm=T),
         max20=max(stock20,na.rm=T))
simple20_2limd=distinct(simple20_2limd,element,year,LU,LU3,.keep_all=T)
