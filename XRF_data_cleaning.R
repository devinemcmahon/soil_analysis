# XRF data visualization: July 2, 2018

# set your working directory (or open a project containing all the files)
setwd('C:\\Users\\Devin\\Documents\\Soil data')

# useful packages:
library(dplyr)
library(lattice)
library(nlme)

# some functions to make the data usable (see other script)
source('xrf_formatting_functions.R')

# read in the data
xrf01=myreader('std_Mendoza_6-27-16')
xrf02=myreader('std_Mendoza_6-28-2017')
xrf03=myreader('std_McMahon_testEu')
xrf04=myreader('std_McMahon_VgN')
xrf05=myreader('std_McMahon_VgE_2')
xrf06=myreader('std_McMahon_VgE')
xrf07=myreader('std_McMahon_May2018')
xrf08=myreader('std_McMahon_JPE2_plus')
# 9 was a weird duplicate
xrf10=myreader('std_McMahon_ItE1_hot')
xrf11=myreader('std_McMahon_ItE1')
xrf12=myreader('std_McMahon_BoE_It')
xrf13=myreader('std_McMahon_BoE')
xrf14=myreader('std_McMahon_Bo')
xrf15=myreader('std_McMahon_6-5-18')
# 16 was a weird duplicate
xrf17=myreader('std_McMahon_6-13-17_3')
xrf18=myreader('std_McMahon_6-13-17_2_fixed')
xrf19=myreader('std_McMahon_6-13-17_1_fixed')
xrf20=myreader('std_McMahon_6-12-17_3')
xrf21=myreader('std_McMahon_6-12-17_2')
xrf22=myreader('std_McMahon_6-12-17')
xrf23=myreader('std_McMahon_5-23_withyear')
xrf24=myreader('std_McMahon_6-27_fixedtypos')
xrf25=myreader('std_Mendoza6-27-17')
xrf26=myreader('McMahon_6-27_7-9-18_fixedtypos')
xrf27=myreader('McMahon_7-13-18_fixed') # rep 5 = pit
xrf28=myreader('McMahon_7-16_19-18_no_spike')
#xrf28=myreader('McM_7-16-18_no_spike') # rep 5 = pit
#xrf29=myreader('McM_7-17-18_no_spike') # rep 5 = pit
xrf29=myreader('JPE2_xelemts_7-20-18')
xrf29=xrf29[,which(names(xrf29) %in% names(xrf28))]
xrf30=myreader('JPNE2_BpE1_xlmts_7-24-18')
xrf30=xrf30[,which(names(xrf30) %in% names(xrf28))]
xrf31=myreader('McMahon_7-26-18_all_xlmts_typosfixed')
xrf31=xrf31[,which(names(xrf31) %in% names(xrf28))]
xrf32=myreader('JPP_EuNplus_7-3-18')
xrf32=xrf32[,which(names(xrf32) %in% names(xrf28))]
xrf33=myreader('EuE1JPA7-31_8-2')
xrf33=xrf33[,which(names(xrf33) %in% names(xrf28))]
xrf34=myreader('Eus_XRF_8-8-18')
xrf34=xrf34[,which(names(xrf34) %in% names(xrf28))]
xrf35=myreader('XRF_9-4-18')
xrf35=xrf35[,which(names(xrf35) %in% names(xrf28))]
xrf36=myreader('McMahon_aband_all')
xrf36=xrf36[,which(names(xrf36) %in% names(xrf28))]
xrf37=myreader('EuE2_8-22-18')
xrf37=xrf37[,which(names(xrf37) %in% names(xrf28))]
#xrf38=myreader('XRF_10-2-18_last3handground')
xrf38=myreader('XRF_10-4-18') # adding 1 
xrf38=xrf38[,which(names(xrf38) %in% names(xrf28))]

sofar=rbind(xrf01,xrf02,xrf03,xrf04,xrf05,xrf06,xrf07,xrf08,
            xrf10,xrf11,xrf12,xrf13,xrf14,xrf15,
            xrf17,xrf18,xrf19,xrf20,xrf21,xrf22,xrf23,xrf24,
            xrf25,xrf26,xrf27,xrf28,xrf29,xrf30,xrf31,xrf32,
            xrf33,xrf34,xrf35,xrf36,xrf37,xrf38)
#head(sofar)


sofar=cleanuptxt(sofar)
#head(sofar)

xrf=xrffun(sofar)


# manual exclusion of other soils
# Skip Cr and TM sites for now (data from just one year)
# now including those too

xrf=droplevels(xrf[-which(xrf$site %in% c('loamy_sand','sandy_clay_loam','clay',
                                          'NIST','silica','Val_mix_0',
                                          'Val_mix_1','Val_mix_2','Val_mix_3',
                                          'Mix1','Mix2')),]) #,'Cr','TM'

unique(xrf$site)
unique(xrf$year)
unique(xrf$depth)
xrf$ID[!is.element(xrf$year,c('04','16'))]
xrf$evaldate[!is.element(xrf$year,c('04','16'))]
xrf$ID=gsub('_2','',xrf$ID)
xrf$ID=gsub('_rerun','',xrf$ID)
xrf$ID[xrf$ID=='Eu.E1.Ee1.4.-0-10.16']='Eu.E1.Ee1.4.0-10.16'

# Add position to stands that were only sampled in one position in 2004
xrf$ID[xrf$stand=='It.E2'&
         xrf$year=='04']=gsub('It.E2','It.E2.L',
                              xrf$ID[xrf$stand=='It.E2'&xrf$year=='04'])
xrf$ID[xrf$stand=='It.E1'&
         xrf$year=='04']=gsub('It.E1','It.E1.T',
                              xrf$ID[xrf$stand=='It.E1'&xrf$year=='04'])
xrf$ID[xrf$stand=='JP.E1'&
         xrf$year=='04']=gsub('JP.E1','JP.E1.E',
                              xrf$ID[xrf$stand=='JP.E1'&xrf$year=='04'])
xrf$ID[xrf$stand=='JP.E2'&
         xrf$year=='04']=gsub('JP.E2','JP.E2.E',
                              xrf$ID[xrf$stand=='JP.E2'&xrf$year=='04'])
xrf$ID[xrf$stand=='Bp.E2'&
         xrf$year=='04']=gsub('Bp.E2','Bp.E2.E',
                              xrf$ID[xrf$stand=='Bp.E2'&xrf$year=='04'])

# looked at difference in composites with and without L (apparent OM contam)
# very similar, just leave the one with all four positions
#   for consistency with other depths
xrf$ID[xrf$ID=='Eu.E2.1.40-60.16_wL']='Eu.E2.1.40-60.16'
xrf=xrf[-which(xrf$ID=='Eu.E2.1.40-60.16_noL'),]

xrf=strfun(xrf)
xrf=droplevels(xrf)

#masses1=read.csv('')
# A little extra data cleaning
#XRF=XRF[XRF$Mass_g>1,] # remove one messed-up sample

#xrf=xrf[!is.na(xrf$ID),]
xrf=xrf[-which(xrf$K==0),] # exported before it was analyzed, maybe?
xrf=xrf[-which(xrf$stand=='BO.P' & xrf$evaldate==as.Date("2018-07-05") &
                 xrf$year=='16'),]
# those ones weren't ground
#xrf=xrf[-which(xrf$ID=='Vg.E.4.40-60.04'),] # Mistyped mass; results off
xrf=xrf[-which(xrf$ID=='Vg.E.4.40-60.04' & xrf$evaldate=="2017-07-13"),]
# Fixed now
xrf=xrf[-which(xrf$Method!='TQ Powders'),] # some analyzed by other method
# saveRDS(xrf,'xrf_through_Jul30_withduplicates.Rds')

# average duplicates
dupIDs=xrf$ID[duplicated(xrf$ID)==T]
dups=xrf[xrf$ID %in% dupIDs,]
dups =dups %>% mutate_if(is.factor, as.character) 
#sapply(dups,class)
#View(dups)
dups1=group_by(dups,ID)%>%
  summarise_all(function(x){x[1]})
dups2=group_by(dups,ID)%>%
  summarise_all(function(x){x[2]})
#dups3=group_by(dups,ID)%>%
#  summarise_all(function(x){x[3]})
#dups3=dups3[!is.na(dups3$Al),]
#xrf$K[xrf$ID %in% dups3$ID]
# just 2 samples, third evaluation similar to others; drop third one
length(dups1$ID)/length(xrf$ID) # 3.7%
summary(abs(dups1$P-dups2$P)*2/(dups1$P+dups2$P)) 
# median 15%, min 1%, max 86%
summary(abs(dups1$Ca-dups2$Ca)*2/(dups1$Ca+dups2$Ca)) 
# med 8% mean 14%
summary(abs(dups1$Mg-dups2$Mg)*2/(dups1$Mg+dups2$Mg)) #13%
summary(abs(dups1$K-dups2$K)*2/(dups1$K+dups2$K)) 
# median just 2.5%, nice
summary(abs(dups1$Fe-dups2$Fe)*2/(dups1$Fe+dups2$Fe)) 
# median just 2.8

dups1$ID
dups1$evaldate
dups2$evaldate
dups2$P/dups1$P # those are really different...
dups2$K/dups1$K # up to 12% different
# oh, it's a typo: two consecutive JP.E1.3.0-10.04s (one is rep 4)
# also two consecutive BO.E.3.20-40.16s on July 9th
# fixed now, but two Eu.E2s and a JP.E2 still >10% different; most within 4%
dups2$Fe/dups1$Fe
# It.E1s: before and after overheating (3% different, not in a consistent direction)
# Rerun on same date: very consistent
# Reanalyzed in November vs July 2018: quite different for Mg and Ca
# Just because those numbers are very small
# Re-run of JP.N: K doesn't change, but Fe does, weirdly
(dups2$Ca-dups1$Ca)*10000
dups2$Fe[dups2$stand=='Eu.N'&dups2$year=='04']/
  dups1$Fe[dups1$stand=='Eu.N'&dups1$year=='04']
plot(I(dups1$Ca/dups2$Ca)~I(dups1$Zr/dups2$Zr))
dups1$ID[dups1$Ca/dups2$Ca>1.2|dups1$Ca/dups2$Ca<.8]

recentdups=group_by(dups,ID)%>% mutate(latestdate=max(evaldate))
recentdups=recentdups[recentdups$latestdate>as.Date('2018-07-20'),]
rec1=dups1[dups1$ID %in% recentdups$ID,]
rec2=dups2[dups2$ID %in% recentdups$ID,]
recdups=merge(rec1,rec2,by='ID')
recdups[,c('ID','P.x','P.y','K.x','K.y','Ca.x','Ca.y','Mg.x','Mg.y',
           'Fe.x','Fe.y')]
recdups[,c('ID','evaldate.x','evaldate.y')]
# re-runs of same sample on 7-24 are pretty close-- 
#   except that P and Mg about halve for JP.E2.2.60-100, and Fe increases
# for redos of Eu.N.04s, some very close, others change a lot (Ca in 1.20-40)
# save these and quantify the variability between duplicates

dups$evaldate=as.character(dups$evaldate)
dupsavg=group_by(dups,ID)%>%
  summarise_all(function(x){ifelse(is.numeric(x)==T,
                                   mean(x,na.rm=T),x[1])})
dupsavg$evaldate=as.Date(dupsavg$evaldate)

XRF=xrf[-which(xrf$ID %in% dupIDs),]

XRF=rbind(XRF,dupsavg)
XRF=mutate(XRF,xrfavgd=as.numeric(is.element(ID,dupIDs)))
# ok

# Mark which ones were analyzed this summer
XRF$recent=as.numeric(XRF$evaldate>as.Date('2018-06-20'))

#write.csv(XRF,'XRFdata_through_7-16-18.csv')
#write.csv(XRF,'XRFdata_through_7-20-18.csv')
#saveRDS(XRF,'XRFdata_through_7-20-18.Rds')
#saveRDS(XRF,'XRFdata_through_7-25-18.Rds')
#saveRDS(XRF,'XRFdata_through_7-26-18.Rds')
#saveRDS(XRF,'XRFdata_through_7-27-18.Rds')
#saveRDS(XRF,'XRFdata_through_7-30-18.Rds')
#saveRDS(XRF,'XRFdata_through_8-9-18.Rds')
#saveRDS(XRF,'XRFdata_through_9-11-18.Rds')
#saveRDS(XRF,'XRFdata_through_10-4-18.Rds')
