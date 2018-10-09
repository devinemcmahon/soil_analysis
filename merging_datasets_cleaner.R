setwd('C:\\Users\\Devin\\Documents\\Soil data')
library(dplyr)
library(lattice)
library(nlme)

# some functions to make the data usable (see other script)
source('xrf_formatting_functions.R')

ea17=read.csv('2017_EA_data.csv')
head(ea17)
length(ea17$ID)

#ea18=read.csv('CN18_to7-27.csv')
#ea18=read.csv('CN18_to8-17.csv')
#ea18=read.csv('CN18_to9-13.csv')
ea18=read.csv('CN18_to_10-2.csv')

ea17=ea17[,-which(!is.element(names(ea17),names(ea18)))]
ea=rbind(ea17,ea18)
unique(ea$stand)
#ea=droplevels(ea[-which(ea$site%in%c('TM','Cr')),]) # come back to these later
ea=ea[,-which(names(ea) %in% c('C.mg','N.mg'))]

length(ea$ID)
length(unique(ea$ID))

# average duplicates
dupIDs=ea$ID[duplicated(ea$ID)==T]
dups=ea[ea$ID %in% dupIDs,]
dups = dups %>% mutate_if(is.factor, as.character) 
dups1=group_by(dups,ID)%>%
  summarise_all(function(x){x[1]})
dups2=group_by(dups,ID)%>%
  summarise_all(function(x){x[2]})

plot(dups1$C,dups2$C/dups1$C) 
dups1$ID[dups2$C/dups1$C>1.2]
dups2$C[dups2$C/dups1$C>1.2]
dups1$C[dups2$C/dups1$C>1.2]
dups2$CNevalday[dups2$C/dups1$C>1.2]
dups1$CNevalday[dups2$C/dups1$C>1.2]
plot(dups1$C,dups2$C)
abline(0,1) # mostly pretty ok, some scatter where C very low 
# maybe 2 points way off, dups2 << dups1
dups1$ID[dups2$C/dups1$C<.8]
dups1$C[dups2$C/dups1$C<.8]
dups2$C[dups2$C/dups1$C<.8]
plot(C~position,data=ea[ea$stand=='Eu.E1'&ea$depth==5,])
plot(dups1$N,dups2$N)
abline(0,1) 
# same deal; dups2 more likely to be < dups1 
# does that mean N was volatilizing in the oven? Probably not
#   some are run on the same day
dups1$ID[dups2$N/dups1$N<.8]
# take the more recent values for JP.A because those were foiled
#   a really long time before running and some may have leaked?
plot(dups1$N,dups2$N,col=as.numeric(as.factor(dups1$site)))
legend('topleft',col=as.numeric(as.factor(unique(dups1$site))),
       legend=unique(dups1$site),bty='n',pch=16)
abline(0,1) # pretty ok
# It and JP cause the most issues
plot(dups1$C,dups2$C,col=as.numeric(as.factor(dups1$site)))
legend('topleft',col=as.numeric(as.factor(unique(dups1$site))),
       legend=unique(dups1$site),bty='n',pch=16)
abline(0,1) # nice
# just average them together.

dupsavg=group_by(dups,ID)%>%
  summarise_all(function(x){ifelse(is.numeric(x)==T,
                                   mean(x,na.rm=T),x[1])})

ea=ea[-which(ea$ID %in% dupIDs),]
ea=rbind(ea,dupsavg)
ea=mutate(ea,CNavgd=as.numeric(is.element(ID,dupIDs)))
length(unique(ea$ID))==length(ea$ID)

CN04=read.csv('2004_CN_data.csv')
names(CN04)=c('ID','N','C','BD')

CN16=strfun(ea) # remake year column

CN04str=strfun(CN04)
# correct between methods of CN analysis (EA vs Walkley-Black)

eawbs=merge(CN16,CN04str,by=c('stand','depth','rep','year'),
            suffixes=c('.16','.04'))
head(eawbs)
length(unique(eawbs$ID.16))

# When position is specified in 2004, match it to 2016 value
eawbs$position.04=as.character(eawbs$position.04)
eawbs$position.16=as.character(eawbs$position.16)
eawbs=eawbs[eawbs$position.04==eawbs$position.16 | 
               is.na(eawbs$position.16)==1,]
length(unique(eawbs$ID.16))==length(eawbs$ID.16)


# Plots
#######
palette(rainbow(length(levels(eawbs$CNevalday))))
plot(Nadj~N.04,data=eawbs,col=CNevalday,pch=16)
abline(0,1)
legend('bottomright',col=as.numeric(as.factor(levels(eawbs$CNevalday))),
       pch=16,legend=levels(eawbs$CNevalday),ncol=2)
palette('default')

plot(C.16~C.04,data=eawbs)
abline(0,1) # pretty close to right on; 16 a little higher?
plot(C.16~C.04,data=eawbs[eawbs$C.16<1,]) #ok
#######

Cmethod.lm=lm(C.16~C.04,data=eawbs)
summary(Cmethod.lm) # R-sq .93, slope = 1.09
Cbetas=summary(Cmethod.lm)$coef[,1]


Nmethod.lm=lm(Nadj~N.04,data=eawbs)
summary(Nmethod.lm) # slope is 1.3 instead of 1.1 for C
Nbetas=summary(Nmethod.lm)$coef[,1]


## More visualization (can skip this)
#########
plot(Nadj~N.04,data=eawbs)
abline(0,1)
plot(Nadj~predict(Nmethod.lm),data=droplevels(eawbs),col=stand,pch=16)
abline(0,1)
legend('bottomright',col=as.numeric(as.factor(
  levels(droplevels(eawbs$stand)))),
       pch=16,legend=levels(droplevels(eawbs$stand)),ncol=2)

# JP.E2 all fall below the 1:1 line (predicted too low)

plot(Nadj~N.04,data=eawbs[eawbs$site.04=='JP',],col=stand,pch=16)
abline(0,1)
plot(Nadj~predict(Nmethod.lm,newdata=eawbs[eawbs$site.04=='JP',]),
     data=eawbs[eawbs$site.04=='JP',],col=stand,pch=16)
abline(0,1) # correction makes it worse for JP.E2, better for the others
# apply same correction to all for now, look at results from JP.E2
#   and check later if correcting them differently changes anything.
#########

CN04$measured_in=rep('04')
CN16$measured_in=rep('16')
CN04$CNavgd=rep(0)
CN04$CNevalday=rep(NA)
CN04$run.date=rep(NA)
CN04$N.old=CN04$N
CN04$N=CN04$N.old*Nbetas[2]+Nbetas[1]
CN04$C.WB=CN04$C
CN04$C=CN04$C.WB*Cbetas[2]+Cbetas[1]

names(CN16)[names(CN16)=='N']='N.old'
names(CN16)[names(CN16)=='Nadj']='N'
CN16$C.WB=rep(NA)

CN16_2=CN16[,which(is.element(names(CN16),names(CN04)))]
CN04_2=CN04[,which(is.element(names(CN04),names(CN16_2)))]
head(CN16_2)
allCN=data.frame(rbind(CN16_2,CN04_2))
head(allCN)
allCN=group_by(allCN,ID) %>% mutate(newest=max(
  as.numeric(as.character(measured_in))))
allCN$oldmeas=as.numeric(as.character(allCN$measured_in))
allCN=allCN[allCN$oldmeas==allCN$newest | is.na(allCN$oldmeas),] 
allCN=ungroup(allCN)

#saveRDS(allCN,'allCN_through_7-22-18.Rds') 
#saveRDS(allCN,'allCN_through_9-13-18.Rds') 
#  If you want to come back to just CN data
#allCN=readRDS('allCN_through_7-27-18.Rds')
# saveRDS(allCN,'allCN_through_10-2-18.Rds')

#xrf=readRDS('xrf_through_8-9-18.Rds')
#xrf=readRDS('XRFdata_through_8-9-18.Rds')
#xrf=readRDS('XRFdata_through_9-11-18.Rds')
xrf=readRDS('XRFdata_through_10-4-18.Rds')

dats=merge(allCN,xrf,by='ID',all=T)
dats=strfun(dats)
dats=eltfun(dats)
head(dats)

BD=readRDS('allBD_9-29-18.Rds')

dats=merge(dats,BD,by=c('stand','depth','year'),all.x=T,all.y=F)
head(dats)
#saveRDS(dats,'all_data_8-17-18.Rds')
#saveRDS(dats,'all_data_9-13-18.Rds')
#saveRDS(dats,'all_data_9-29-18.Rds')
#saveRDS(dats,'all_data_10-4-18.Rds')
