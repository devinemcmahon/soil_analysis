## A critical look at bulk density changes
setwd('C:\\Users\\Devin\\Documents\\Soil data')


# Bulk density:
#BD=read.csv('revised_BD_8-16-17.csv')
BD=read.csv('revised_BD_7-19-18.csv')
head(BD)
BD$ID=paste(BD$Sample,'16',sep='.') 
#	note that this column will not be accurate in the final dataset
BD=strfun(BD)
BD=eltfun(BD)
summary(BD$elt)
head(BD)
BD$Density_gm.cm3[which(is.na(BD$Density_gm.cm3))]=
  BD$Interp_density[which(is.na(BD$Density_gm.cm3))]
BD$Density_gm.cm3[BD$ID=='Vg.N.B.60-100.16']=
  BD$Interp_density[BD$ID=='Vg.N.B.60-100.16'] 
# soil fell out of this one but the other sample at this depth was 
# 	exceptionally dense; this weights the dense intact one 75%
# 	because the interpolated value is .5*N.A + .5*N.B


plot(Density_gm.cm3~as.numeric(stand),data=BD[BD$depth==5,],
     col=elt,pch=16,las=1)
palette(rainbow(24))
plot(Density_gm.cm3~as.numeric(elt),data=BD[BD$depth==5,],
     col=stand,pch=16,xaxt='n',las=1)
axis(side=1,labels=levels(BD$elt),at=seq(1,6))
# looks like density generally less in L and more in E for 0-10 cm
# just as much difference between A and B, maybe
legend('bottomright',ncol=5,pch=15,cex=.8,
       legend=levels(BD$stand),col=as.numeric(as.factor(
         levels(BD$stand))))
plot(Density_2004_uncorr~as.numeric(Position_2004),data=BD[BD$depth==5,],
     col=stand,pch=16,xaxt='n',las=1)
# measured in E in 2004: JP.E1 and E2, Eu.E2, Bp.E2
# In L: It.E1, It.E2
# In both: Bp.E1 (again, L is less dense)
# average BD in 2016 > 04: Eu (including in N), JP except E1, It.E
# 2016 < 04 in BO.E, Bp.E1 (E in 2016 is low), Vg.E a little
# Oh, Eu values for 2004 are average of 2 stands
# Likely reasons for changes: only sampling in L in It in 2004
# Wrong values in Eu in 2004
# That guy in Bp.E1 doing a bad job in 2016? 
#     Or just having incomplete sampling in 2016; missing values for E
#  T is also low at some depths (also low for It at 80 cm)

BDaov=aov(Density_gm.cm3~elt,
          data=droplevels(BD[BD$elt %in% c('E','L','T'),]))
summary(BDaov) # p=.08
TukeyHSD(BDaov) # no difference; E marginally denser than L and T
BD5aov=aov(Density_gm.cm3~elt,
          data=droplevels(BD[BD$elt %in% c('E','L','T') & 
                               BD$depth==5,]))
TukeyHSD(BD5aov) # also no difference
# E tends to be densest though, L least dense

# To fix this: just use the 2016 values for Eu both years? 
# Just the L values for It?
# Use 2016 values for everything? Should be more representative
#   due to more reps
# Except in BO and Vg? Changes there could be real
#   or some sampling weirdness in 2004

palette('default')

BD=BD[order(BD$depth),]
xyplot(depth~Density_gm.cm3|stand,groups=elt,type='l',ylab='Depth (cm)',
       data=droplevels(BD[BD$elt %in% c('E','L','T'),]),
       ylim=c(90,0),as.table=T,
       par.settings = list(superpose.line = list(col = c(1,2,3),lwd = 2)),
       auto.key=list(space='top', columns=3,lines=TRUE, points=FALSE))

# BO and Vg: no pattern, braided lines
# Bp.E2, Eu.E2, JPs: lower density in L, especially near surface
# It positions pretty similar, but E denser in top 20-30 cm
# L less dense in JP.A and TM.E, ok
# Cr.E: L > T, huh
xyplot(depth~Density_gm.cm3|stand,groups=elt,type='l',ylab='Depth (cm)',
       data=BD, ylim=c(90,0),as.table=T,
       par.settings = list(superpose.line = list(col = c(4,5,6,1,2,3),lwd = 2)),
       auto.key=list(space='top', columns=3,lines=TRUE, points=FALSE))

plot(Density_gm.cm3~Density_2004_uncorr,data=BD,pch=as.character(elt),col=stand)
abline(0,1)

plot(Density_gm.cm3~Density_2004_uncorr,data=BD,
     pch=as.character(elt),col=Position_2004)
# Lower-density sites were sampled in L in 2004, more likely to increase
# L samples from 2016 still above the line (denser than 2004) in those stands
#   but E especially likely to be denser
# E also more likely to be denser than 2004 in stands sampled in E originally
# but not necessarily (everything less dense in BO)

BDchgaov=aov(log(Density_gm.cm3/Density_2004_uncorr)~elt,
             data=droplevels(BD[BD$elt %in% c('E','L','T'),]))
summary(BDchgaov) # not signif diff
BDchgaov2=aov(log(Density_gm.cm3/Density_2004_uncorr)~Position_2004,
             data=droplevels(BD[BD$Position_2004 %in% c('E','L','T'),]))
summary(BDchgaov2) # this one is different
# but is that an effect of position in 2004, or low BD in 2004?


# average the different positions and quantify variability
BDsum16=group_by(BD,stand,depth) %>% summarise(
  avgBD=mean(Density_gm.cm3),BDsd=sd(Density_gm.cm3))
BDsum16$year=rep('16')
BDsum04=group_by(BD,stand,depth) %>% summarise(
  avgBD=mean(Density_2004_uncorr,na.rm=T),
  BDsd=sd(Density_2004_uncorr,na.rm=T))
BDsum04$year=rep('04')
BDboth=data.frame(rbind(BDsum04,BDsum16))
BDboth$BDsd[BDboth$BDsd==0]=NA

