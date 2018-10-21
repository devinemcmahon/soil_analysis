# Soil stock analysis

#setwd('C:\\Users\\Devin\\Documents\\Soil data')
source('soil_data_reader.R')

# Plots by year and depth in each stand
trellis.par.set(strip.background = list(col = 'grey80'),
                par.strip.text=list(cex=.8))

xyplot(depth~mn|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='C' ,],
       ylim=c(90,0),xlab='C (g / 100 g)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~mn*10|stand,groups=year,type='l',ylab='Depth (cm)',
        data=datsmnok[datsmnok$element=='N' ,],
        ylim=c(90,0),xlab='N (g / kg)',as.table=T,
        par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
        auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
# With separate correction for JP.E2, N increases at depth there as well
# Don't use that

xyplot(depth~I(mn/1000)|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca2' ,],
       ylim=c(90,0),xlim=c(-.2,1.2),xlab='Ca (g / kg)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~I(mn)|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca2' ,],
       ylim=c(90,0),xlim=c(-80,500),xlab='Ca (mg / kg)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~I(mn/1000)|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='Ca2' &
                       datsmnok$site %in% c('Eu','JP') &
                       datsmnok$stand!='JP.N',],
       ylim=c(90,0),xlab='Ca (g / kg)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~mn/1000|stand,groups=year,type='l',ylab='Depth (cm)',
        data=datsmnok[datsmnok$element=='Zr' ,],
        ylim=c(90,0),xlab='Zr (g / kg)',as.table=T,
        par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
        auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

xyplot(depth~value*10|stand,groups=year,type='p',ylab='Depth (cm)',
       data=dats[dats$element=='N' ,],
       ylim=c(90,0),xlab='N (g / kg)',as.table=T,
       par.settings = list(superpose.points = list(col = c(2,4))),
       auto.key=list(space='top', columns=2,lines=F, points=T))
xyplot(depth~value|stand,groups=year,type='p',ylab='Depth (cm)',
       data=dats[dats$element=='K' & !is.element(dats$site,
                                                 c('Bp','TM','Cr')),],
       ylim=c(90,0),xlab='K (mg / kg)',as.table=T,
       par.settings = list(superpose.points = list(col = c(2,4))),
       auto.key=list(space='top', columns=2,lines=F, points=T))
xyplot(depth~value/1000|stand,groups=rep,type='p',ylab='Depth (cm)',
       data=dats[dats$element=='K' & is.element(dats$site,
                                                 c('Bp','TM','Cr')) &
                   dats$year=='16',],
       ylim=c(90,0),xlab='K (g / kg) 2016',as.table=T,
       par.settings = list(superpose.points = list(col = c(2,4))),
       auto.key=list(space='top', columns=2,lines=F, points=T))
# Huge spatial variation in K in the high-K sites
# also very slightly depleted in top 10 cm?

xyplot(depth~value|stand,groups=year,type='p',ylab='Depth (cm)',
       data=dats[dats$element=='S' ,],
       ylim=c(90,0),xlab='S (mg / kg)',as.table=T,
       par.settings = list(superpose.points = list(col = c(2,4))),
       auto.key=list(space='top', columns=2,lines=F, points=T))

xyplot(depth~mn/1000|stand,groups=year,type='l',ylab='Depth (cm)',
       data=datsmnok[datsmnok$element=='S' ,],
       ylim=c(90,0),xlab='S (g / kg)',as.table=T,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2)),
       auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))

plot(P~S,data=widedats4) # no apparent relationship

# Where did a given element change? T-tests by stand and depth, n=4 per year
gen_ttable(ttests[ttests$element=='N',]$depth,
           ttests[ttests$element=='N',]$stand,
           ttests[ttests$element=='N',]$pval,
           ttests[ttests$element=='N',]$tstat,0.05)

# Do the unlikely changes in JP.N come out significant?
gen_ttable(ttests[ttests$stand=='JP.N',]$depth,
           ttests[ttests$stand=='JP.N',]$element,
           ttests[ttests$stand=='JP.N',]$pval,
           ttests[ttests$stand=='JP.N',]$tstat,0.05)
# Yes, especially N, S, Fe, Zr, Nb, Mo (things that shouldn't change)


# Across all sites, does C accumulate over time?
allC20.lme=lme(stock20~year,random=~1|site/stand,
               data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20.lme) # yes, increases (marginally, adding Eu)
# with latest data, not even marginally
plot(resid(allC20.lme)~dats2deps$stock20[dats2deps$element=='C' &
                                           !is.na(dats2deps$stock20)])
# residuals still highest at highest C
plot(resid(allC20.lme)~dats2deps$year[dats2deps$element=='C' &
                                        !is.na(dats2deps$stock20)])
# more spread in 04? nah
qqnorm(resid(allC20.lme))
qqline(resid(allC20.lme)) # upper tail off but probably ok


allC20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20LU.lme) # yes, increases, but not in pasture
# now marginal increase (p=0.065), no difference in pasture (itrxn p=0.10)
# with fixed BD, increase overall (.052), decrease in pasture, p = 0.024
qqnorm(resid(allC20LU.lme))
qqline(resid(allC20LU.lme)) # probably ok? upper tail quite off
plot(resid(allC20LU.lme)~dats2deps$year[dats2deps$element=='C' &
                                          !is.na(dats2deps$stock20)]) #ok
plot(resid(allC20LU.lme)~dats2deps$LU[dats2deps$element=='C' &
                                        !is.na(dats2deps$stock20)]) 
# less ok, tiny spread for pasture
# which way is correct?
# more similar with 20 cm and fixed BD



allC100LUsite.lme=lme(stock100~year*LU,random=~1|site,
                      data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC100LUsite.lme) # p for increase =.06, year:pasture = .06 also
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

allC100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                     data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC100LU.lme) 
AIC(allC100LUsite.lme,allC100LUstd.lme,allC100LU.lme)
# both = lowest, but almost identical to site alone
# fix this to take out It.E1 and It.N

#allN100LU.lme=lme(log(stock100)~year*LU,random=~1|site/stand,
#                  data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
allN100LU.lme=lme(log(stock100)~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='N',],na.action = na.omit)
summary(allN100LU.lme) 
allN20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
qqr(allN20LU.lme) #ok
summary(allN20LU.lme)
# 20 cm: increases only in natveg (signif intrxn);
#   within site not stand, both natveg and interaction signif
# no change with constant bulk density
# 100 cm within site, only intrxn signif again
# with Eu, now increase is significant and larger in N 
# should these be proportional increases?
#  now still signif increase, but decrease in P, no change in N,
#   with or without It.E1 and N
qqr(allN100LU.lme)# tails way off, even with log

plot(resid(allN20LU.lme)~dats2deps$year[dats2deps$element=='N' &
                                          !is.na(dats2deps$stock20)]) #ok
plot(resid(allN20LU.lme)~dats2deps$LU[dats2deps$element=='N' &
                                        !is.na(dats2deps$stock20)]) #eh
plot(resid(allN100LU.lme)~test2deps$year[test2deps$element=='N' &
                                          !is.na(test2deps$stock100)])#ok
# different means, similar spreads
plot(resid(allN100LU.lme)~test2deps$LU[test2deps$element=='N' &
                                        !is.na(test2deps$stock100)]) 
# how to get around the difference in the unbalanced land uses?
plot(resid(allN100LU.lme)~dats2deps$stand[dats2deps$element=='N' &
                                            !is.na(dats2deps$stock100)])
plot(stock100~stand,data=dats2deps[dats2deps$element=='N',])
# variance differs within groups; is that a problem?

# To 20 cm, K decreases at p=0.07 if Bp included, increases if excluded
# P decreases only in native veg (JP.N weirdness); same result for P2
# no change in noisy Ca, but increase at p=.062 if using Ca2
# no model has the full contingent of 128 observations--still missing data

allP100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='P2',],na.action = na.omit)
summary(allP100LU.lme) 
# now p=0.1 for increase between years, still decreases in native though
#   that one crazy value for JP.N
# residuals ugly, log no help

allK100LU.lme=lme(log(stock100)~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='K' &
                                   test2deps$site!='Bp' &
                                   test2deps$stand!='JP.N',],na.action = na.omit)
summary(allK100LU.lme) # with Bp: decreases in 2016, 
# but nearly increases in native veg in 2016
# without: increases overall (fixed BD: no change), more (only) in native veg
# residuals not as bad as some others, but still a big outlier
# lower AIC with stand within site vs just stand
# no change without JP.N (down to 11 stands)

allCa100LU.lme=lme(log(stock100)~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='Ca2',],na.action = na.omit)
summary(allCa100LU.lme) # increases, but not in natveg


allNcratLU.lme=lme(concratio~year*LU,random=~1|site/stand,
                   data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allNcratLU.lme) # concratio decreases with time
# residuals pretty normal
allNcratLUstd.lme=lme(concratio~year*LU,random=~1|stand,
                      data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
AIC(allNcratLU.lme,allNcratLUstd.lme) # lower for just stand

allCacratLU.lme=lme(concratio~year*LU,random=~1|stand,
                    data=dats2deps[dats2deps$element=='Ca',],na.action = na.omit)
summary(allCacratLU.lme) # increases, nice; by more in euc, much less in nat
# more in pasture to start with 
# K: no change over time for euc, but signif year*natveg intrxn 
#   (distribution gets deeper in nat, because of JP.N weirdness)

allNcratLU2.lme=lme(concrat2~year*LU,random=~1|stand,
                    data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allNcratLU2.lme) # barely decreases, p=.054 (now it's .011, why?)


Nrat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
#                  data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
  data=test2deps[test2deps$element=='N',],na.action = na.omit)
summary(Nrat100LU.lme) # residuals normalish
# proportion in top 20 cm decreases between years (p=.047); no LU effects
# no change without It.E1 and N

Crat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='C',],na.action = na.omit)
summary(Crat100LU.lme) # gets shallower in pasture, maybe
Prat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='P',],na.action = na.omit)
summary(Prat100LU.lme) # no change
Krat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='K',],na.action = na.omit)
summary(Krat100LU.lme)# K: gets shallower in pasture, deeper in natveg 
Krat100LU2.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='K' &
                                   dats2deps$site!='Bp',],na.action = na.omit)
summary(Krat100LU2.lme) # same deal

#dats2deps=mutate(dats2deps,biome=ifelse(site %in% c('BO','Vg','Eu'),'AF','Cer'))


allC20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
#allC20bm.lme=lme(stock20~year*LU*biome,random=~1|site/stand,
             data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
summary(allC20bm.lme) # decreases in AF, increases in Cerrado
# probably not enough data for LU*biome to be very credible
# difference in what vegetation is in each biome
# Just keep eucs and nats?

qqr(allC20bm.lme)
plot(resid(allC20bm.lme)~dats2deps$year[dats2deps$element=='C' &
                                          !is.na(dats2deps$stock20)]) 
# bigger spread in 16
plot(resid(allC20bm.lme)~dats2deps$biome[dats2deps$element=='C' &
                                        !is.na(dats2deps$stock20)])
allC20bmen.lme=lme(stock20~year*biome*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='C' &
                                  dats2deps$LU!='P',],na.action = na.omit)
summary(allC20bmen.lme) # increases in Cerrado, esp in native

allC20bmen2.lme=lme(log(stock20)~year*biome*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='C' &
                                    test2deps$LU!='P',],na.action = na.omit)
summary(allC20bmen2.lme) 
# still increases in 16 in Cerrado, but nothing else is significant
# log transform with the residuals
plot(resid(allC20bmen2.lme)~as.factor(test2deps$biome[test2deps$element=='C' &
                                              test2deps$LU!='P'&
                                           !is.na(test2deps$stock20)])) 
# No patterns by biome, land use, or C stock

# Group tstock by biome to get one obs per stand and calc SE off that
tstock=mutate(tstock,LU2=ifelse(LU=='E','E','O'))

bmstock=group_by(tstock,element,LU,biome) %>%
  summarise(bmstock100_16=mean(stock100_16),
            bmstock100_04=mean(stock100_04),
            bmstock20_16=mean(stock20_16),
            bmstock20_04=mean(stock20_04),
            sd100_16=sd(stock100_16),
            sd100_04=sd(stock100_04),
            sd20_16=sd(stock20_16),
            sd20_04=sd(stock20_04),
            ngrp=n())
# standard errors are huge

shortbmstk=droplevels(bmstock[bmstock$element %in% 
                              c('C','N','P2','K','Ca2','Mg2','Fe','S','Al',
                                'Cl','Nb','Zr'),])
yrdiffstockplot20_bmLU(bmstock[bmstock$element=='N',])
# no that's no good. standard errors huge.
yrdiffstockplot20_bmLUall(tstock[tstock$element=='C',],fulllegend = T)
yrdiffstockplot100_bmLUall(tstock[tstock$element=='C',],fulllegend = T)

yrdiffstockplot20_bmLUall(tstock[tstock$element=='Ca2',],label=F,fulllegend = T)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')
yrdiffstockplot20_bmLUall(tstock[tstock$element=='Ca2' & tstock$stock20_16<1,],label=F)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')

yrdiffstockplot20_bmLUall(tstock[tstock$element=='Zn' & tstock$stock20_16<.1,])
tstock$stand[tstock$element=='Zn' & tstock$stock20_16>.1] # Bp.E2 has huge Zn incr


# Ratios for analyzing change in different vegetation types


# Repeat analysis with just eucs
# Subset data frames now make in soil_data_reader
#euc2deps=droplevels(dats2deps[dats2deps$LU=='E',])
#test2deps=dats2deps[-which(dats2deps$stand %in% c('It.E1','It.N')),]
#euc2deps2=droplevels(test2deps[test2deps$LU=='E',]) 

eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                    data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme) # increase in Cerrado only (also if using euc2deps2)
qqr(eucC20bm.lme) # again, log tranform helps a bit, increases significance

eucN20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) # no change
qqr(eucN20bm.lme)

eucCa20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='Ca2',],na.action = na.omit)
qqr(eucCa20bm.lme)
summary(eucCa20bm.lme) # maybe increase in AF (p=.07), increase in Cer (.0005)

eucK20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='K'&
                                  euc2deps$site!='Bp',],na.action = na.omit)
qqr(eucK20bm.lme)
summary(eucK20bm.lme) # increases in AF only

eucS20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='S',],na.action = na.omit)
qqr(eucS20bm.lme) # tails way off; bimodal distribution of S?
summary(eucS20bm.lme) # Cerrado starts lower, increases (p=.053)
eucS100bm.lme=lme(log(stock100)~year*biome,random=~1|site/stand,
                  data=euc2deps2[euc2deps2$element=='S',],na.action = na.omit)
# tails still off, esp. upper; same result but lower p-vals as for 20 cm

eucZn20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='Zn',],na.action = na.omit)
qqr(eucZn20bm.lme) # upper tail off
summary(eucZn20bm.lme) # increases in Cerrado (also to 100 cm)

eucCu20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Cu',],na.action = na.omit)
qqr(eucCu20bm.lme) # lower tail off
summary(eucCu20bm.lme) # no change

eucMn20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Mn',],na.action = na.omit)
qqr(eucMn20bm.lme) # both tails off, log possibly worse
summary(eucMn20bm.lme) # increases in Cerrado (Bp stands, also in JP.N)
#tstock$stand[tstock$element=='Mn' & tstock$stock20_16>.3]

eucC100.lme=lme(stock100~year,data=euc2deps2[euc2deps2$element=='C',],
                random=~1|site/stand,na.action=na.omit)
summary(eucC100.lme) # increases (but only with It.E1)
# now also increases without It.E1 (p=.03)
qqnorm(resid(eucC100.lme))
qqline(resid(eucC100.lme)) # upper tail still off (not bad without It.E1)
plot(resid(eucC100.lme)~euc2deps2$year[euc2deps2$element=='C' &
                                           !is.na(euc2deps2$stock100)])
# smaller IQR in 2016; just 40 obs/yr # ok without It.E1

eucN100.lme=lme(stock100~year,data=euc2deps2[euc2deps2$element=='N',],
                random=~1|site/stand,na.action=na.omit)
summary(eucN100.lme) # increases quite a bit, p = 0.0041 # nope
# With 2016 BD in both years and no It.E1, N increases (p=.0003. now .0004?)
qqnorm(resid(eucN100.lme))
qqline(resid(eucN100.lme)) # tails way off; ok without It.E1 and fixed BD?
# something changed when I copied everything over for version control...

# P: no change, p=0.10 (increasing tendency, when using P2)
#     but that's the crazy outlier? residuals horrible
# K decreases with Bp in, increases without (p = .027, increase is small)
# No change in P, or K without Bp
# no consistent trend in Ca; residual distrib has crazy tails
eucCa100.lme=lme(log(stock100)~year,data=euc2deps2[euc2deps2$element=='Ca2',],
                random=~1|site/stand,na.action=na.omit)
summary(eucCa100.lme) # increases at p < 0.0001
# residuals normal with log transform
# with fixed bulk denisity, no change at 20 cm, but increase if Eu removed
eucK100bm.lme=lme(log(stock100)~year*biome,
                data=euc2deps2[euc2deps2$element=='K'&euc2deps2$site!='Bp',],
                 random=~1|site/stand,na.action=na.omit)
qqr(eucK100bm.lme)
summary(eucK100bm.lme) # increases in AF, maybe decreases in Cerrado

eucC100bm.lme=lme(stock100~year*biome,data=euc2deps2[euc2deps2$element=='C',],
                random=~1|site/stand,na.action=na.omit)
qqr(eucC100bm.lme)
summary(eucC100bm.lme) # increases in Cerrado only (p=.035)

eucCa100bm.lme=lme(log(stock100)~year*biome,
                   data=euc2deps2[euc2deps2$element=='Ca2',],
                 random=~1|site/stand,na.action=na.omit)
summary(eucCa100bm.lme) # increases at p < 0.0001 in Cerrado only

eucN100bm.lme=lme(stock100~year*biome,
                   data=euc2deps2[euc2deps2$element=='N',],
                   random=~1|site/stand,na.action=na.omit)
summary(eucN100bm.lme) # increases in AF, less in Cerrado



# just-euc ratios
Crat100euc.lme=lme(stockratio~year,random=~1|site/stand,
                  data=euc2deps2[euc2deps2$element=='C',],na.action = na.omit)
summary(Crat100euc.lme) # no change
qqr(Crat100euc.lme) # a few outliers, esp at upper tail (log no help)
Nrat100euc.lme=lme(stockratio~year,random=~1|site/stand,
                   data=euc2deps2[euc2deps2$element=='N',],na.action = na.omit)
summary(Nrat100euc.lme) # decrease, p=.0304 #now no change
# residuals less bad
Krat100euc.lme=lme(stockratio~year,random=~1|site/stand,
                   data=euc2deps2[euc2deps2$element=='K',],na.action = na.omit)
# marginally shallower, p=.0612; tails also far off; log no help
# with 16 densities, p=.074 increase
# now it's .066, with version control and redone N content. Why the change???
# no change with euc2deps2?
Krat100euc2.lme=lme(stockratio~year,random=~1|site/stand,
                   data=euc2deps[euc2deps$element=='K'&
                                   euc2deps$site!='Bp',],na.action = na.omit)
# this is p=.03

# Ca: no change but residuals look pretty normal
# When using Ca2, large increase (makes sense), p=0.0005
# no change in P; weird residuals even when log-transformed

eucC20.lme=lme(stock20~year,data=euc2deps[euc2deps$element=='C',],
                random=~1|site/stand,na.action=na.omit)
summary(eucC20.lme) 
# when It.E1 included, marginal increase (p=.071); without, p=.678 
# increase at surface would not be due to OM falling down core--real?
# drive by surface bulk density change--is that real?
# Without bulk density change, significant increase (p=.007, now .03)
#   but .18 without It.E1
qqr(eucC20.lme) # upper tail still off
# No change in N at all
# P = weird residual with or without log
# Ca increases with log
eucK20.lme=lme(log(stock20)~year,data=euc2deps[euc2deps$element=='K' &
                                            euc2deps$site!='Bp',],
               random=~1|site/stand,na.action=na.omit)
# nice residual distrib with log; increase at p=.035

allN20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                 data=dats2deps[dats2deps$element=='N' & 
                               dats2deps$stand!='JP.N',],na.action = na.omit)
summary(allN20LU.lme) # no change with fixed BD
# excluding JP.N, natveg N term p = .1, but still no change 
qqr(allN20LU.lme) # pretty good


allCa20LU.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                   data=dats2deps[dats2deps$element=='Ca2',],
                  na.action = na.omit)
qqr(allCa20LU.lme) # really normal with log
summary(allCa20LU.lme)
# significant increase with time in euc/ overall 
# smaller increase in pasture (which was more to start at p=.05) at p=.08
# no change in native (p <.0001, coefficient opposite to yr16 coef)
# same results (qualitatively) if JP.N removed
allCa20LU2.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                   data=dats2deps[dats2deps$element=='Ca2' &
                        dats2deps$site!='Eu',],na.action = na.omit)
summary(allCa20LU2.lme) # same deal, but lower p-values

allP20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='P2',],# &
                                #   dats2deps$stand!='JP.N',],
                  na.action = na.omit)
qqr(allP20LU.lme) # no good with or without log
# no change with or without JP.N

allK20LU.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='K' &
                                    #  dats2deps$stand!='JP.N'&
                                   dats2deps$site!='Bp',],
                  na.action = na.omit)
# without Bp, upper tail off a bit but not so bad with log
summary(allK20LU.lme) # increase overall, some decrease in natveg

# I do want to compare change in ratio in diff land uses 
Nrat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='N',],na.action = na.omit)
summary(Nrat100LU.lme) 
qqr(Nrat100LU.lme) # residuals normalish
# proportion in top 20 cm decreases between years; no LU effects
# Without It.E1 and It.N, nothing changes

Crat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='C',],na.action = na.omit)
summary(Crat100LU.lme) # gets shallower in pasture, maybe
# now gets deeper in native (marginally)
Prat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='P2',],na.action = na.omit)
summary(Prat100LU.lme) # no change (deeper in pasture? Not if using P2)
Krat100LU.lme=lme(log(stockratio)~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='K',],na.action = na.omit)
summary(Krat100LU.lme)# K: gets shallower in pasture, deeper in natveg 
# now just deeper in natveg (JP.N?)
Krat100LU2.lme=lme(stockratio~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='K' &
                                    test2deps$site!='Bp',],na.action = na.omit)
summary(Krat100LU2.lme) # same deal
qqr(Krat100LU.lme) # way off, even with log
Carat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=test2deps[test2deps$element=='Ca2',],na.action = na.omit)
summary(Carat100LU.lme) 
# gets shallower overall, started shallower and gets deeper in pasture,
#   no change in native (signif, coeffs cancel w yr, smaller p with Ca2)
Zrrat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='Zr',],na.action = na.omit)
summary(Zrrat100LU.lme) # gets deeper in N because of weird JP probably


qqr(Carat100LU.lme) # nice
# Wait, these are proportion data and don't come from normal distributions
# Try generalized mixed models
#library(MASS)
Krat100LU2.pql=glmmPQL(stockratio~year*LU,random=~1|site/stand,
                       data=test2deps[test2deps$element=='K' &
                                        test2deps$site!='Bp',],
                       na.action = na.omit,family='quasibinomial')
summary(Krat100LU2.pql) # same deal as lme
Krat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site2,
                       data=simple20[simple20$element=='K',],
                       na.action = na.omit,family='quasibinomial')
Crat100LU.pql=glmmPQL(stockratio~year*LU,random=~1|site/stand,
                       data=test2deps[test2deps$element=='C',],
                       na.action = na.omit,family='quasibinomial')
# parameter estimates differ from lme, but p-values very similar


Carat100LU.pql=glmmPQL(stockratio~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='Ca2',],
                   na.action = na.omit,family='quasibinomial')
summary(Carat100LU.pql) 
qqr(Carat100LU.pql) # not quite as good as lme but pretty similar

Krateuc2.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                       data=euc2deps2[euc2deps2$element=='K' &
                                        euc2deps2$site!='Bp',],
                       na.action = na.omit,family='quasibinomial')
summary(Krateuc2.pql) # same deal as lme, increases at p=.07
# with Bp included, same result, somewhat smaller increase

Nrateuc2.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='N',],
                     na.action = na.omit,family='quasibinomial')
summary(Nrateuc2.pql)

Crat100bm.aov=aov(stockratio~biome*LU,
                      data=test2deps[test2deps$element=='C',],na.action = na.omit)

# variance vs mean
plot(I(sdrat_16^2)~rat_16,data=tstock[tstock$element=='C',])
# No real patterns except for P and P2

yrdiffratplot_LU(tstock[tstock$element=='C' & 
                          tstock$stand!='It.E1' & tstock$stand!='It.N',])
hist(test2deps$stockratio[test2deps$element=='C'])
densityplot(~stockratio|year,auto.key=T,
            data=test2deps[test2deps$element=='C',],groups=LU)


# what if we take out any changes in It.E1 and It.N that could be due to
#   surface contamination?
# start with just taking out It.E1 and N altogether

test2deps=dats2deps[-which(dats2deps$stand %in% c('It.E1','It.N')),]
mostC100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='C',],na.action = na.omit)
summary(mostC100LU.lme) # no change except marginal decrease in pasture
# with 2016 BD, pasture decreases, so does N (p=.027)
qqr(mostC100LU.lme) # ok except 1 outlier
mostN100LU.lme=lme(stock100~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='N',],na.action = na.omit)
summary(mostN100LU.lme) # no change except marginal decrease in pasture
# also no change except increase in N still (JP.N)
# No change in Ca, K (except increase in N), P (ex decr N)
# With 2016 bulk density for everything, N decreases, P decreases more

CN100LU.lme=lme(conc100~year*LU,random=~1|site/stand,
                   data=test2deps[test2deps$element=='CN',],na.action = na.omit)
qqr(CN100LU.lme) # a few huge outliers
plot(conc100~stand,data=test2deps[test2deps$element=='CN',],las=2)

# does bulk density change significantly?
BD100LU.lme=lme(BD100~year*LU,random=~1|site/stand,
                data=dats2deps[dats2deps$element=='K',],na.action = na.omit)
summary(BD100LU.lme) # increases in native veg, where it was ~lower in 04
BD20LU.lme=lme(BD20~year*LU,random=~1|site/stand,
                data=dats2deps[dats2deps$element=='K',],na.action = na.omit)
summary(BD20LU.lme) # lower in N, increases in P, ok



yrdiffstockplot100_LU(tstock[tstock$element=='N',])
#yrdiffstockplot100_LU(tstock[tstock$element=='N' & 
#                               !is.element(tstock$stand,c('It.E1','It.N')),])
legend('bottomright',pch=16,bty='n',#cex=1.6,
       col=c('blue3','springgreen','darkgoldenrod1'),
       legend=c('Eucalyptus','Native vegetation','Pasture'))
legend('right',pch=seq(20,25), bty='n',#cex=1.6,
       legend=c('BO','Bp','It','JP','Eu','Vg'),pt.bg=1)
yrdiffstockplot100_LU(tstock[tstock$element=='C',])
yrdiffstockplot100_LU(tstock[tstock$element=='P2',])
legend('topleft',bty='n',cex=1.8,
       legend='P (Mg / ha)\n0-100 cm')

yrdiffstockplot100_LU(tstock[tstock$element=='K'&tstock$stand!='Bp.E1',])

plot(value~depth,data=dats[dats$stand=='Bp.E2'&dats$element=='K',],
     col=rep,pch=as.numeric(as.factor(year)))
# only reps 2 and 3 have the really high values
# 1 and 4 very similar to the 2004 values
unique(widedats$ID[widedats$stand=='Bp.E2'&
                     widedats$depth==50&!is.na(widedats$K)])
# for rep 1, two row positions (L and Ee2) as well as the composite
# no rep 3 
# values for the two positions pretty close to those for composite, good
# some of the CN data are missing: reps 2-4 of 40-60 and rep 1 of 20-40
#   the four samples for 40-60 are all the positions of rep 1, oops
# fixed as of October 4
# Also missing rep 1 of 2004 for 40-60 and rep 4 for 60-100
# reps 1 and 4 have way more C than 2 and 3, which are close to 2004 values


# Spatial heterogeneity
plot(value~depth,data=dats[dats$stand=='Bp.E1'&dats$element=='K',],
     col=rep,pch=as.numeric(as.factor(elt)))
# reps, not row positions, are variable (row positions similar within a rep)
# also true for E2
plot(value~as.numeric(elt),data=dats[dats$depth==5&dats$element=='K'&
                                       dats$elt %in% c('E','L','T'),],
     col=stand,pch=rep)
plot(value~as.numeric(elt),data=dats[dats$depth==5&dats$element=='K'&
                                       dats$elt %in% c('E','L','T') &
                                       dats$site!='Bp',],
     col=stand,pch=rep)
xyplot(value~as.numeric(elt)|stand,groups=rep,
       data=dats[dats$depth==5&dats$element=='K'& dats$year=='16' &
                   dats$elt %in% c('E','L','T') &dats$site!='Bp',],
       ylim=c(0,600),pch=15)
xyplot(value~as.numeric(elt)|stand,groups=rep,
       data=dats[dats$depth==5&dats$element=='P'& dats$year=='16' &
                   dats$elt %in% c('E','L','T'),],ylim=c(0,600),pch=19)
elt.lme=lme(log(value)~elt, random=~1|element/stand,
            data=dats[dats$depth==5 & dats$year=='16' &
                        dats$elt %in% c('E','L','T') &
                        dats$element %in% c('C','N','P2','Ca2','K'),],na.action=na.omit)
# No way that distribution of residuals will be normal
Celt.lme=lme(log(value)~elt, random=~1|stand,
            data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                        dats$element =='C',],na.action=na.omit)
qqr(Celt.lme) # Nice (with log) 
summary(Celt.lme) # L significantly less than E
xyplot(value~as.numeric(elt)|stand,groups=rep,
       data=dats[dats$depth==5&dats$element=='C'& dats$year=='16' &
                   dats$elt %in% c('E','L','T'),],pch=19)
Celtbm.lme=lme(log(value)~elt*biome, random=~1|stand,
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='C',],na.action=na.omit)
summary(Celtbm.lme) # difference with L is only in Cerrado; I think driven by It.E1
# where L is substantially < T for 2 of 4 reps
Celtbm15.lme=lme(log(value)~elt*biome, random=~1|stand,
                 data=dats[dats$depth==15 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                             dats$element =='C',],na.action=na.omit)
summary(Celtbm15.lme) # from 10-20 cm, L > E in AF only (represented only by Eu.E2)
# Cerrado*L term is significant and opposite 
# no general pattern
xyplot(value~as.numeric(elt)|stand,groups=rep,
       data=dats[dats$depth==15&dats$element=='C'& dats$year=='16' &
                   dats$elt %in% c('E','L','T'),],pch=19)
Cdepeltbm.lme=lme(log(value)~elt*biome, random=~1|stand/depth,
               data=dats[dats$year=='16' & dats$elt %in% c('E','L','T') &
                           dats$element =='C',],na.action=na.omit)
qqr(Cdepeltbm.lme) # not as good
summary(Cdepeltbm.lme) # still difference is just in L and Cerrado; driven by top 5 cm
# best to just look at that layer
xyplot(value~depth|stand,groups=elt,
       data=dats[dats$element=='C'& dats$year=='16' &
                   dats$elt %in% c('E','L','T'),],pch=19)
plot(value~depth,data=dats[dats$element=='C'& dats$year=='16' & dats$stand=='It.E1' &
                             dats$elt %in% c('E','L','T') &dats$depth<30,],
     col=as.numeric(elt)-2,pch=rep+14)

plot(value~depth,data=dats[dats$element=='Ca2'& dats$year=='16' & dats$stand=='It.E1' &
                             dats$elt %in% c('E','L','T') &dats$depth<30,],
     col=as.numeric(elt)-2,pch=rep+14) 
legend('topright',bty='n',pch=15,col=as.numeric(unique(
  dats$elt[dats$elt %in% c('E','L','T')]))-2,
  legend=unique(dats$elt[dats$elt %in% c('E','L','T')])) # more Ca in linha, less in toco
# In Bp.E1, more Ca in E--one super high value at surface
# Bp.E2 has more in T; one rep has very high Ca in T at surface
# Eu.E2 mixed, maybe lower for E, smaller spread for L and some high T
# Eu.E1? maybe higher in L 
plot(value~depth,data=dats[dats$element=='Ca2'& dats$year=='04' & dats$stand=='Eu.E1',],
      col=rep,pch=18) # 1 and 3 are ~E (30 cm from planting line) and
#   2 and 4 are ~T (60 cm); 4 has highest value at surface and 2 second
#   below 20 cm, though, 4 has highest and 3 second-highest
# For C, 1 and 2 are similar, 4 has highest values at surface
# 

Nelt.lme=lme(value~elt, random=~1|stand,
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='N',],na.action=na.omit) # log makes qqr worse
summary(Nelt.lme) # L lower again at p=.06
Neltbm.lme=lme(value~elt*biome, random=~1|stand,
               data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                           dats$element =='N',],na.action=na.omit)
summary(Neltbm.lme) # L < E for Cerrado at p=.028; T < E at .072


Kelt.lme=lme(log(value)~elt, random=~1|stand,
              data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                          dats$element =='K',],na.action=na.omit)
qqr(Kelt.lme)
summary(Kelt.lme) # no difference
Pelt.lme=lme(value~elt, random=~1|stand,
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='P2',],na.action=na.omit)
qqr(Pelt.lme) # still not very good (log worse)
summary(Pelt.lme) # no difference

Caelt.lme=lme(log(value)~elt, random=~1|stand,
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='Ca2',],na.action=na.omit)
qqr(Caelt.lme)
summary(Caelt.lme) # no difference
# Spatial heterogeneity overall:
depcvs=group_by(dats4,stand,LU,biome,depth,element,year) %>%
  summarise(CV=sd(repval,na.rm=T)/mean(repval,na.rm=T))
tapply(depcvs[depcvs$element=='C',]$CV,depcvs[depcvs$element=='C',]$LU,summary)
tapply(depcvs[depcvs$element=='C',]$CV,depcvs[depcvs$element=='C',]$depth,summary)
tapply(droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])$CV,
       droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])$element,summary)
dep16cvs=group_by(droplevels(dats[dats$position %in% c('Ee1','Ee2','L','T') &
                                    dats$LU!='A' & dats$year=='16',]),
                  stand,LU,biome,depth,element) %>%
  summarise(CV=sd(value,na.rm=T)/mean(value,na.rm=T))
tapply(droplevels(dep16cvs[dep16cvs$element %in% c('C','N','P2','Ca2','K'),])$CV,
       droplevels(dep16cvs[dep16cvs$element %in% c('C','N','P2','Ca2','K'),])$element,summary)

plot(CV~LU,data=depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])
plot(CV~LU,data=depcvs[depcvs$element=='C',])
plot(CV~depth,data=depcvs[depcvs$element=='C' & depcvs$LU!='A',],col=LU,pch=5,cex=2)

bwplot(CV~LU|element,data=depcvs[depcvs$element %in% c('C','N','K') &
                                   depcvs$LU!='A',])
bwplot(CV~as.factor(depth)|LU,data=depcvs[depcvs$element %in% c('C','N','K') &
                                            depcvs$stand!='It.E1'&depcvs$stand!='It.N'&
                                            depcvs$LU!='A',])
# CV doesn't change a lot with depth; IQR actually larger at deeper depths
#   maybe due to influence of pit samples? Yes, less of an effect without It.E1
cvaov=aov(CV~LU,data=depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])
summary(cvaov)
qqr(cvaov) # nooo

yrdiffstockplot20_LU(tstock[tstock$element=='K'&tstock$site!='Bp',])
yrdiffstockplot20_LU(tstock[tstock$element=='N',])
legend('bottomright',pch=c(20,24,25,21,22,23), bty='n',ncol=2,
       legend=c('Atlan-','tic','Forest','Cer-','ra-','do'),pt.bg=1)
yrdiffstockplot20_LU(tstock[tstock$element=='Ca2' &
                              !is.element(tstock$stand,c('Eu.E1','Eu.E2')),])
yrdiffstockplot20_LU(tstock[tstock$element=='Ca2' &tstock$stock20_04<.5 &
                              tstock$stock20_16<.5,])
yrdiffstockplot20_LU(tstock[tstock$element=='S',])


# Mean stocks
tapply(shorttstk$stock100_16,shorttstk$element,
       function(x){mean(x,na.rm=T)})
tapply(shorttstk$stock100_16,shorttstk$element,
       function(x){sd(x,na.rm=T)})
mean(shorttstk$BD100_16[shorttstk$element=='C'])
sd(shorttstk$BD100_16[shorttstk$element=='C'])

shortE = shorttstk[shorttstk$LU=='E' &
                     shorttstk$stand!='It.E1',]
tapply(shortE$stock100_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock100_16,shortE$element,
       function(x){sd(x,na.rm=T)})
tapply(shortE$rat_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$rat_16,shortE$element,
       function(x){sd(x,na.rm=T)})
mean(shortE$BD100_16[shortE$element=='C'])
sd(shortE$BD100_16[shortE$element=='C'])

shorttstk$site=as.character(shorttstk$site)
shorttstk=mutate(shorttstk,
                 site2=ifelse(stand=='JP.E2'|stand=='JP.N','JP2',site),
                  LU2=ifelse(LU=='E','E','O'))
mrstkchgs=group_by(droplevels(shorttstk[shorttstk$element %in% shortE$element,]),
                 stand,site,LU,site2,LU2,element,biome)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04))
mrstkE=mrstkchgs[mrstkchgs$LU=='E' & mrstkchgs$stand!='It.E1',]
mrstkO=mrstkchgs[mrstkchgs$LU!='E' & mrstkchgs$stand!='It.N'&
                   mrstkchgs$stand!='JP.N',]
mrstkcomp=merge(mrstkE,mrstkO,by=c('biome','site2','element'),
                suffixes = c('_E','_O'),all=F)
mrstkE2=mrstkchgs[mrstkchgs$LU=='E',]
mrstkO2=mrstkchgs[mrstkchgs$LU!='E',]
mrstkcomp2=merge(mrstkE2,mrstkO2,by=c('biome','site2','element'),
                suffixes = c('_E','_O'),all=F)

t.test(mrstkcomp$chgln20_E[mrstkcomp$element=='C'],
       mrstkcomp$chgln20_O[mrstkcomp$element=='C'],
       paired=T) # p = .019, the log ratios are not the same
t.test(mrstkcomp2$chgln20_E[mrstkcomp2$element=='C'],
       mrstkcomp2$chgln20_O[mrstkcomp2$element=='C'],
       paired=T) # p = .235, the log ratios are the same
hist(mrstkchgs$chgln20[mrstkchgs$element=='C']) 
# There isn't a lot of data here
# Probably best to focus on the pairs
# but maybe do the lmes because that way you can see 
#   increase or decrease for each vegetation type
# or:
t.test(mrstkcomp2$chgln20_E[mrstkcomp2$element=='C'])
t.test(mrstkcomp2$chgln20_O[mrstkcomp2$element=='C'])
# no significant change between years in either

t.test(mrstkcomp$chgln20_E[mrstkcomp$element=='P'],
       mrstkcomp$chgln20_O[mrstkcomp$element=='P'],
       paired=T) # no difference for N, P, K, yes for Ca
t.test(mrstkcomp2$chgln20_E[mrstkcomp2$element=='P'],
       mrstkcomp2$chgln20_O[mrstkcomp2$element=='P'],
       paired=T) # no: N, P, K yes: Ca
hist(mrstkchgs$chgln20[mrstkchgs$element=='P']) 
qqnorm(mrstkchgs$chgln20[mrstkchgs$element=='P']) 
# ok for Ca, 2 weird K outliers, P still weird

t.test(mrstkcomp2$chgln20_E[mrstkcomp2$element=='P']) 
#p=.08 for N, .04 for Ca, no for K and P
t.test(mrstkcomp2$chgln20_O[mrstkcomp2$element=='P'])
# no change in N, Ca, K, P

tapply(mrstkE$chgrt100,mrstkE$element,mean)            
tapply(mrstkE$chgrt100,mrstkE$element,median)            
tapply(mrstkE$chgrt100,mrstkE$element,sd) 

tapply(mrstkO$chgrt100,mrstkO$element,mean)            
tapply(mrstkO$chgrt100,mrstkO$element,median)            
tapply(mrstkO$chgrt100,mrstkO$element,sd) 
plot(chgrt100~LU,data=mrstkchgs[mrstkchgs$element=='C',])

simpstkchgs=mrstkchgs[mrstkchgs$stand %in% simple20_2$stand,]
simpstkchgs$site=as.character(simpstkchgs$site)
simpstkchgs=mutate(simpstkchgs,site2=ifelse(stand=='JP.E2'|stand=='JP.N','JP2',site),
                  LU2=ifelse(LU=='E','E','O'))
Cyrrat.lme=lme(chgln20~LU2,random=~1|site2,
               data=simpstkchgs[simpstkchgs$element=='C',])
qqr(Cyrrat.lme)
summary(Cyrrat.lme) # no signif effect
# just reformat the data and do a paired t-test



tapply(shorttstk$stock100_16[shorttstk$LU=='N'&
                               shorttstk$site!='It'],
       shortE$element[shorttstk$LU=='N'&
                        shorttstk$site!='It'],
       function(x){mean(x,na.rm=T)})
# This isn't really helpful because there are unpaired eucs: site effect
tapply(shorttstk$rat_16[shorttstk$site!='It'],
       shorttstk$element[shorttstk$site!='It'],
       function(x){mean(x,na.rm=T)})

# Difference between euc and not
vegrats=group_by(shorttstk, site, element) %>%
  mutate(natrat100_16=ifelse(sum(LU=='N')>0,
                             stock100_16/stock100_16[LU=='N'],NA),
         natrat100_04=ifelse(sum(LU=='N')>0,
                             stock100_04/stock100_04[LU=='N'],NA),
         natratrat=natrat100_16/natrat100_04)
summary(vegrats$natrat100_04[vegrats$element=='C'&vegrats$LU=='E'])
summary(vegrats$natrat100_16[vegrats$element=='C'&vegrats$LU=='E'])
# generally similar C in euc and native
t.test(log(vegrats$natrat100_04[vegrats$element=='C'&vegrats$LU=='E'])) 
# p = .88
t.test(log(vegrats$natratrat[vegrats$element=='C'&vegrats$LU=='E']))
# p = .10, maybe increase relative to native between 04 and 16 
#   (CI -.01 to .12)

# is this different from the lmes? 
t.test(tstock$stock100_04[tstock$element=='C'],
       tstock$stock100_16[tstock$element=='C'],paired=T)
# p = .41, but .027 for N
# different because just 1 sample per stand, not 4
# is this more correct? probably not, doesn't include variance within site

summary(vegrats$natratrat[vegrats$element=='C'&vegrats$LU=='E'])
# median 1.07, mean 1.06, max 1.15
summary(vegrats$natratrat[vegrats$element=='N'&vegrats$LU=='E'])
# median 1.29, mean 1.08


tstock$stock100_16[tstock$stand=='Eu.E2'&tstock$element=='N']
# E1 is the one that increased a bunch
# Should have had a loss of 115 kg ha-1; instead gained some 3500 kg
#   well, that isn't right
# Increases at depth = big stock changes when multiplied out
# Eu.E2 was where I was cutting off the fallen OM, I think
# Bioturbation?? Leaching in lame stands?

# looking at the data to compare to budgets
#budgets=mutate(budgets,inoutrat=In_kgha/Out_kgha)
# usually close except for P and Mg
# a change in either input or output should affect obs-predicted agreement
plot(Concentration~Nutrient,data=budgets)

pagano=data.frame(Egrandconc=c(0.00118814,0.000030795,0.00037971,
                    0.001193941,0.000128281,0.00000841287,0.000077024),
                  Nutrient=c('N','P','K','Ca','Mg','B','S'))
budgets=merge(budgets,pagano,by='Nutrient')
budgets=mutate(budgets,budget=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                 Concentration*511)/1000,
                 stdconcbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                                Egrandconc*511)/1000,
                 denserbudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                              Egrandconc*560)/1000,
                 lessdensebudg=(In_kgha_1+In_kgha_2-(Wood_m3_1+Wood_m3_2)*
                             Egrandconc*460)/1000,
                 lessrotbudg=(In_kgha_1+In_kgha_2*Past_harvests_known-
                                (Wood_m3_1+Wood_m3_2*Past_harvests_known)*
                                Egrandconc*511)/1000
               )

plot(Concentration~Egrandconc,data=budgets,pch=as.character(Nutrient),col=Nutrient)
abline(0,1)
# K and N vary a lot; one very low Ca value = BO.E (no bark measured)
# My measured N is almost always lower than Pagano estimate
# So if using their estimate, even larger N decreases would be predicted
# Plantar data has much less N, more K than my estimates
budgets$Stand[budgets$Concentration<.0008 & budgets$Nutrient=='Ca']
plot(budget~stdconcbudg,data=budgets,pch=as.character(Nutrient),col=Nutrient)
abline(0,1)
# N and sometimes K influence the budget a lot.
plot(denserbudg~lessdensebudg,data=budgets,
     pch=as.character(Nutrient),col=Nutrient)
abline(0,1) # pretty close; only matters for C and N

# variation in percent bark could also be important



shorterstk=merge(shorttstk,budgets,by.x=c('stand','element'),
                 by.y=c('Stand','Nutrient'))

stkchgs=group_by(droplevels(shorterstk),stand,element,biome)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04),
            #budget=Budget/1000,
            budget=budget,
            stdconcbudg=stdconcbudg, denserbudg=denserbudg,
            lessdensebudg=lessdensebudg,conc=Concentration)
stkchgs2=stkchgs[stkchgs$stand!='It.E1',]
t.test(stkchgs$chg20,stkchgs$budget,paired=T) # for all elements, p=.069
t.test(stkchgs$chg20[stkchgs$element=='N'],
       stkchgs$budget[stkchgs$element=='N'],paired=T)
# Ca sort of differs at p=.1, NPK don't 

palette(rainbow(9))
plot(chg100~budget,data=stkchgs[stkchgs$element=='N',],
     pch=16,col=stand)
legend('bottomright',pch=15,col=as.factor(levels(stkchgs$stand)),
       legend=levels(stkchgs$stand),bty='n')
abline(0,1) # yeah, not even related

knowns=stkchgs[-which(stkchgs$stand %in% 
                        c('Bp.E1','BO.E','JP.E1','JP.E2')),]
plot(chg100~budget,data=knowns[knowns$element=='Ca',],
     pch=16,col=stand)

plot(chg100~budget,data=stkchgs,type='n', 
     xlab='Fertilizer - harvest, Mg ha-1',
     ylab='Observed change in stocks to 100 cm',las=1)
rect(xleft=-1, ybottom=-1, xright=2, ytop=2,border='gray50')
abline(0,1)
abline(h=0,lty=3)
abline(v=0,lty=3)
text(stkchgs$budget,stkchgs$chg100,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))


plot(chg20~budget,data=stkchgs,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',
     #ylim=c(-1,1.8),
     las=1)
rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
segments(x0=stkchgs$minbudg,x1=stkchgs$maxbudg,y0=stkchgs$chg20,
         col=as.numeric(stkchgs$stand))
segments(x0=stkchgs$budget,y0=stkchgs$chg20-stkchgs$sdchg20/sqrt(3),
         y1=stkchgs$chg20+stkchgs$sdchg20/sqrt(3),
         col=as.numeric(stkchgs$stand))
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
#text(stkchgs$stdconcbudg,stkchgs$chg20,labels=stkchgs$element,
#text(stkchgs$lessrotbudg,stkchgs$chg20,labels=stkchgs$element,
      #cex=stkchgs$conc*1000, 
      col=as.numeric(stkchgs$stand))
legend('bottomright',pch=15,col=as.factor(levels(stkchgs$stand)),
       legend=levels(stkchgs$stand),bty='n',ncol=2)
# Concentration matters? Higher estimated wood associated with larger N increases
#   than expected (lower concs would make expected losses less in Vg.E and Eu.E2)
# And lower N conc in BO.E could be associated with underestimated expected losses?
#   Number of harvests/harvested volume more important
# but large observed Ca losses in Eu.E2 associated with high Ca concs 
#       (i.e. estimated losses should be large, too)
# High Ca concs in Bp.E1 and It.E1 also should offset large expected increases,
#     but expected still >> observed
# Low Ca concs in JP associated with larger-than-predicted Ca increases
#   so again concentration doesn't seem to be the main driver. Harvested biomass, yes?
# Or maybe they left the bark onsite in JP so increases > expected?
# Changing the concentrations to the literature values I cited for E. grandis
#   improves Ca in BO.E, but doesn't help much with other nutrients

# taking out the second rotation helps for Ca, K in BO.E and Bp.E1, not for JPs or N

plot(chg20~budget,data=stkchgs,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1 (lessrot)',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',las=1,
     xlim=c(-.2,.5),ylim=c(-.2,.5))
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=stkchgs$minbudg,x1=stkchgs$maxbudg,y0=stkchgs$chg20,
         col=as.factor(stkchgs$stand))
segments(x0=stkchgs$budget,y0=stkchgs$chg20-stkchgs$sdchg20/sqrt(3),
         y1=stkchgs$chg20+stkchgs$sdchg20/sqrt(3),
         col=as.factor(stkchgs$stand))
text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
#text(stkchgs$lessrotbudg,stkchgs$chg20,labels=stkchgs$element,
          cex=stkchgs$conc*2000,
     col=as.numeric(stkchgs$stand))

palette('default')
plot(chg20~budget,data=stkchgs,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',las=1)
rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
     col=as.numeric(as.factor(stkchgs$biome))+2)
legend('bottomright',pch=15,col=c(3,4),
       legend=c('Atlantic Forest','Cerrado'),bty='n')



stkchgs2=stkchgs[stkchgs$stand!='It.E1',]
plot(chg100~budget,data=stkchgs,type='n', 
     xlab='Fertilizer - harvest, Mg ha-1',
     ylab='Observed change in stocks to 100 cm',las=1,
     ylim=c(-1,2),xlim=c(-1,2))
text(stkchgs$budget,stkchgs$chg100,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))
legend('bottomright',pch=15,col=as.factor(levels(stkchgs$stand)),
       legend=levels(stkchgs$stand),bty='n',ncol=2)


# Plots showing average change for different nutrients?
tapply(stkchgs$chgrt100,stkchgs$element,mean) # each stand gets weight of 1
tapply(stkchgs$chgrt20,stkchgs$element,mean)
barplot(tapply(stkchgs2$chgln20,stkchgs2$element,mean))
sds20=tapply(stkchgs2$chgln20,stkchgs2$element,function(x){sd(x,na.rm=T)})

plot(unique(as.numeric(as.factor(stkchgs2$element))),
     tapply(stkchgs2$chgln20,stkchgs2$element,mean), ylim=c(-.4,2),
     las=1,xaxt='n',xlab='Element',ylab='Change in 20 cm stock (log ratio)')
segments(x0=seq(1,length(unique(stkchgs2$element))),
         y0=tapply(stkchgs2$chgln20,stkchgs2$element,mean)-
           tapply(stkchgs2$chgln20,stkchgs2$element,sefun),
         y1=tapply(stkchgs2$chgln20,stkchgs2$element,mean)+
           tapply(stkchgs2$chgln20,stkchgs2$element,sefun))
abline(h=0,lty=2)
axis(side=1,at=seq(1,length(unique(stkchgs2$element))),
     labels = unique(stkchgs2$element))

sefun=function(x){sd(x)/sqrt(length(x)-1)}

plot(unique(as.numeric(as.factor(stkchgs2$element))),
     tapply(stkchgs2$chgln100,stkchgs2$element,mean), ylim=c(-0.5,1.5),
     las=1,xaxt='n',xlab='Element',ylab='Change in 100 cm stock (log ratio)')
segments(x0=seq(1,length(unique(stkchgs2$element))),
         y0=tapply(stkchgs2$chgln100,stkchgs2$element,mean)-
           tapply(stkchgs2$chgln100,stkchgs2$element,sefun),
         y1=tapply(stkchgs2$chgln100,stkchgs2$element,mean)+
           tapply(stkchgs2$chgln100,stkchgs2$element,sefun))
abline(h=0,lty=2)
axis(side=1,at=seq(1,length(unique(stkchgs2$element))),
     labels = unique(stkchgs2$element))

bwplot(chgln100~element,data=stkchgs2)

simp2deps=droplevels(dats2deps[dats2deps$stand %in% 
                                c('BO.E','BO.P','Vg.E','Vg.N',
                                  'JP.E1','JP.N','It.E1','It.N'),])
simp2deps=mutate(simp2deps,LU2=ifelse(LU=='E','E','O'))

simpstks=droplevels(shorttstk[shorttstk$stand %in% 
                                 c('BO.E','BO.P','Vg.E','Vg.N',
                                   'JP.E1','JP.N','It.E1','It.N'),])
simpstks=mutate(simpstks,LU2=ifelse(LU=='E','E','O'))


simpchgs=group_by(simpstks,stand,LU2,element)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04))

bwplot(chgln20~element|LU2,data=simpchgs)
bwplot(chgln20~LU2|element,data=simpchgs) # neither very informative



# Power tests
summary(tstock$sd100_16[tstock$element=='C']/
          tstock$stock100_16[tstock$element=='C']) 
# 1.5% to 18%, mean and median 8.7%

summary(tstock$sd100_16[tstock$element=='N']/
          tstock$stock100_16[tstock$element=='N']) 
# median 7.2%, mean 8.6% (1.3 to 25%)

# Most elements have med/mean around 8% (more for K if Bp.E2 included)
# Ca: median is 36%

# Anovas: Concentration in top 20 cm
Ccaov=aov(conc20~year*LU,data=dats2deps)
summary(Ccaov) # no, maybe close to a LU effect
library(car)
Anova(Ccaov, type = "III") # not close to LU effect any more
Ccaov_stand=aov(conc20~year*stand,data=dats2deps)
summary(Ccaov_stand) # significant stand effect of course
Ccaov_siteLU=aov(conc20~year*site*LU,data=dats2deps)
summary(Ccaov_siteLU) # only site is significant
# I want mixed effects

# Conc in top 20 cm

Cc20.lme=lme(log(conc20)~year*LU,random=~1|site/stand,
               data=dats2deps[dats2deps$element=='C',],
             na.action = na.omit)
summary(Cc20.lme) # concentration increases significantly
# greatest increase in euc: decreases in pasture, small decrease in N
qqr(Cc20.lme) # nice with log transform
plot(resid(Cc20.lme)~dats2deps$conc20[dats2deps$element=='C' &
                                        !is.na(dats2deps$conc20)]) #ok

plot(resid(Cc20.lme)~dats2deps$site[dats2deps$element=='C' &
                                        !is.na(dats2deps$conc20)]) 
# pretty ok
euc2deps=droplevels(dats2deps[dats2deps$LU=='E',])
test2deps=dats2deps[-which(dats2deps$stand %in% c('It.E1','It.N')),]
euc2deps2=droplevels(test2deps[test2deps$LU=='E',]) 

Cc20euc.lme=lme(log(conc20)~year,random=~1|site/stand,
                data=euc2deps[euc2deps$element=='C',],
                na.action = na.omit)
summary(Cc20euc.lme) # yes, signif increase
qqr(Cc20euc.lme) # what are the outliers? still there without It.E1
qqnorm(resid(Cc20euc.lme),pch=as.character(
  dats2deps$site[dats2deps$element=='C' &
                   !is.na(dats2deps$conc20)]))
# some I and B: not sure

Nc20.lme=lme(log(conc20)~year*LU,random=~1|site/stand,
             data=dats2deps[dats2deps$element=='N',],
             na.action = na.omit)
summary(Nc20.lme) # no change but maybe small increase in N (p=.07)
qqr(Nc20.lme)

Nc20euc.lme=lme(log(conc20)~year,random=~1|site/stand,
                data=euc2deps[euc2deps$element=='N',],
                na.action = na.omit)
summary(Nc20euc.lme) # marginal increase, p = .08
qqr(Nc20euc.lme) # pretty good
# Looks like the increases are Vg.E, Eu.E1, It.E1 
# Increased litterfall, recent fertilization?

# P: no change; log-transformed residuals still way not normal
# K: no log and no Bp = increase at p=0.041
# Ca: log transform, big increase

# Next up: superficiality metric based on concentration, 
#   without It.E1 and N? Or with pit data substituted in for those
# Also amended stocks using BD from 2016 in all cases

Ccrat.lme=lme(log(concrat2)~year,random=~1|site/stand,
              data=euc2deps2[euc2deps2$element=='C',],
              na.action = na.omit)
qqr(Ccrat.lme) # tails way off, ok with log
summary(Ccrat.lme) # no change; is concratio w/o log best?

# maybe do stock ratios, but with fixed bulk density?
# that's misleading if change in OM really changed BD
#   could be the case in BO, maybe Vg
# I think it's ok

Ceucrat.lme=lme(log(concratio)~year,random=~1|site/stand,
              data=euc2deps2[euc2deps2$element=='C',],
              na.action = na.omit)
summary(Ceucrat.lme) # no change

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

table(simple20$LU2[simple20$element=='C'],simple20$year[simple20$element=='C']) 
# Not balanced--rep 5
Csimp.lme=lme(stock20~year*LU2,random=~1|site,
               data=simple20[simple20$element=='C',], na.action=na.omit)
qqr(Csimp.lme) # wavery, but not so bad
summary(Csimp.lme) # nothing is significant

Csimp.lme2=lme(stock20~year*LU2,random=~1|site2,
              data=simple20_2[simple20_2$element=='C',], na.action=na.omit)
qqr(Csimp.lme2) # tails off, mostly ok?
summary(Csimp.lme2) # C increases between years if JP.E1 and E2 both included

Nsimp.lme=lme(stock20~year*LU2,random=~1|site,
              data=simple20[simple20$element=='N',], na.action=na.omit)
qqr(Nsimp.lme) # nice
summary(Nsimp.lme) # N initially higher in other, doesn't increase without JP.E1
Nsimp.aov=aov(stock20~year*LU2,data=simple20[simple20$element=='N',])
qqr(Nsimp.aov) # tails off
summary(Nsimp.aov) # not correct I think
Nsimp.lme2=lme(stock20~year*LU2,random=~1|site,
              data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
qqr(Nsimp.lme2) # nice
summary(Nsimp.lme2) # same deal

Ksimp.lme=lme(log(stock20)~year*LU2,random=~1|site,
              data=simple20[simple20$element=='K',], na.action=na.omit)
qqr(Ksimp.lme) # some outliers
summary(Ksimp.lme) # no significant terms; increase in euc when JP excluded

xyplot(stock20~year|site,groups=LU2,data=simple20[simple20$element=='K',])
# different trends for euc and non-euc veg, but not consistent among sites
# For C, apparent increases in euc in most sites and other in It,
#   but decreases in JP.N, Vg.N, BO.P, and Eu.N

xyplot(stock20~year|site,groups=LU2,data=simple20[simple20$element=='N',],
       auto.key=list(space='top', columns=2,lines=FALSE, points=TRUE),
       ylab='N stock to 20 cm (Mg ha-1)',xlab='Sampling year')

xyplot(stock20~year|site2,groups=LU2,data=simple20_2[simple20_2$element=='C',],
       auto.key=list(space='top', columns=2,lines=FALSE, points=TRUE),
       ylab='C stock to 20 cm (Mg ha-1)',xlab='Sampling year')

xyplot(stock20~year|site2,groups=LU2,data=simple20_2[simple20_2$element=='K',],
       auto.key=list(space='top', columns=2,lines=FALSE, points=TRUE),
       ylab='K stock to 20 cm (Mg ha-1)',xlab='Sampling year')


Casimp.lme=lme(log(stock20)~year*LU2,random=~1|site,
              data=simple20[simple20$element=='Ca2',], na.action=na.omit)
qqr(Casimp.lme) # good with log
summary(Casimp.lme) # increase in euc, barely in other when JP.P included (p=.047)

Psimp.lme=lme(stock20~year*LU2,random=~1|site,
               data=simple20[simple20$element=='P2',], na.action=na.omit)
qqr(Psimp.lme) # upper tail off
summary(Psimp.lme) # no signif changes as you might expect

Zrsimp.lme=lme(log(stock20)~year*LU2,random=~1|site,
               data=simple20[simple20$element=='Zr',], na.action=na.omit)
qqr(Zrsimp.lme) 
summary(Zrsimp.lme) # starts lower in other, marginally increases there (p=.07)
# increase significant when JP.P included (p=.045)

Ksimp.lme2=lme(log(stock20)~year*LU2,random=~1|site2,
               data=simple20_2[simple20_2$element=='K',], na.action=na.omit)


Casimp.lme2=lme(log(stock20)~year*LU2,random=~1|site2,
               data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
qqr(Casimp.lme2) # good with log
summary(Casimp.lme2) # increase in euc but not in other (opposite sign)


Psimp.lme2=lme(stock20~year*LU2,random=~1|site2,
              data=simple20_2[simple20_2$element=='P2',], na.action=na.omit)
qqr(Psimp.lme2) # upper tail off
summary(Psimp.lme2) # no signif changes as you might expect

Zrsimp.lme2=lme(log(stock20)~year*LU2,random=~1|site2,
               data=simple20_2[simple20_2$element=='Zr',], na.action=na.omit)
qqr(Zrsimp.lme2) 
summary(Zrsimp.lme2) # minor increase in other (p=.052)

# is it appropriate to group pasture with native veg? 
Casimp.lme3=lme(log(stock20)~year*LU,random=~1|site2,
                data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
qqr(Casimp.lme3) # good with log
summary(Casimp.lme3) # ca starts higher under P, increases as under E
# starts same under N, does not increase
# Probably better to keep the three types separate

# But how to show which are different?
Casimp.aov=aov(log(stock20)~year*LU,
               data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
qqr(Casimp.aov)
summary(Casimp.aov)
TukeyHSD(Casimp.aov) # none of these are actually different
# except N < P and 16 > 04 at p=.058
# I don't want all the year/veg comparisons, just year within veg
# Dunnett test?

Csimp.aov=aov(stock20~year*LU,
               data=simple20_2[simple20_2$element=='C',], na.action=na.omit)
qqr(Csimp.aov) # tails way off, opposite direction from usual
summary(Csimp.aov) # marginal LU effect


Csimp.lme3=lme(log(stock20)~year*LU,random=~1|site2,
               data=simple20_2[simple20_2$element=='C',], na.action=na.omit)
qqr(Csimp.lme3) # tails off, mostly ok? Log is better
summary(Csimp.lme3) # C increases between years if JP.E1 and E2 both included
# But decreases in pasture; no change in native?
exp(0.161)
exp(0.161-.191)
library(sjPlot) # download failed
# this makes nice plots of coefficients of models,
#   but would probably not be good for visualizing many elements at once
fakedat=data.frame(site2=rep(unique(simple20_2$site2),2,each=2),
                   year=rep(c('04','16'),12),
                   LU=c(rep('E',12),rep(c('P','N','N','P','N','N'),each=2)))
fakefix=data.frame(year=rep(c('04','16'),3),LU=rep(c('E','P','N'),each=2))
Csimp20fix=summary(Csimp.lme3)$tTable
Csimp20preds=predict(Csimp.lme3,newdata=fakedat)
fakedat$C20=predict(Csimp.lme3,newdata=fakedat)
fakedatgr=group_by(fakedat,year,LU)%>%mutate(mnC20=mean(C20))
plot(exp(mnC20)~as.numeric(as.factor(LU)),data=fakedatgr,
     col=as.numeric(as.factor(year))*2,pch=16, xaxt='n',xlim=c(.5,3.5),
     las=1,ylab='Predicted C stock to 20 cm',xlab='')
axis(side=1,at=seq(1,3),labels=c('Eucalyptus','Native','Pasture'))
text(seq(1,3),rep(35,3), c('*','','#'),cex=c(1.5,1,1))
#points(seq(1,3),rep(35,3), pch=c(3,NA,4),cex=c(1.5,1,1))
legend('topright',pch=15,col=c(2,4),legend=c('2004','2016'))

simple20_2$LU=factor(simple20_2$LU,levels=c('E','N','P'))
plot(stock20~as.numeric(as.factor(LU)),data=simple20_2[simple20_2$element=='C',],
     col=as.numeric(as.factor(year))*2, xaxt='n',xlim=c(.5,3.5),
     las=1,ylab='Predicted C stock to 20 cm',xlab='')
axis(side=1,at=seq(1,3),labels=c('Eucalyptus','Native','Pasture'))
text(seq(1,3),rep(35,3), c('*','','#'),cex=c(1.5,1,1))
points(exp(mnC20)~as.numeric(as.factor(LU)),data=fakedatgr,
     col=as.numeric(as.factor(year))*2,pch=18,cex=2)
     
#points(seq(1,3),rep(3.5,3), pch=c(3,NA,4),cex=c(1.5,1,1))
legend('topright',pch=15,col=c(2,4),legend=c('2004','2016'))


Nsimp.lme3=lme(stock20~year*LU,random=~1|site,
               data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
qqr(Nsimp.lme3) # nice
summary(Nsimp.lme3) # with all sites included, starts higher in N, no change

Ksimp.lme3=lme(log(stock20)~year*LU,random=~1|site2,
               data=simple20_2[simple20_2$element=='K',], na.action=na.omit)
# no differences between types or years

Psimp.lme3=lme(stock20~year*LU,random=~1|site2,
               data=simple20_2[simple20_2$element=='P2',], na.action=na.omit)
qqr(Psimp.lme3) # both tails off
summary(Psimp.lme3) # no signif changes as you might expect

Zrsimp.lme3=lme(log(stock20)~year*LU,random=~1|site2,
                data=simple20_2[simple20_2$element=='Zr',], na.action=na.omit)
qqr(Zrsimp.lme3) 
summary(Zrsimp.lme3) # minor increase in other (p=.052)
# when separated by veg type, P starts with more Zr, N with less; 
# increases in native 


simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]

Ksimp100.lme=lme(log(stock100)~year*LU,random=~1|site,
              data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme) # mostly ok
summary(Ksimp100.lme) # excluding JP2 and It, 
#   increase overall and both noneuc start higher than euc, but no intrxn 
Nsimp100.lme=lme(log(stock100)~year*LU,random=~1|site,
                 data=simp100[simp100$element=='N',], na.action=na.omit)
qqr(Nsimp100.lme) # ok
summary(Nsimp100.lme)
Csimp100.lme=lme(log(stock100)~year*LU,random=~1|site,
                  data=simp100[simp100$element=='C',], na.action=na.omit)
qqr(Csimp100.lme)
summary(Csimp100.lme) # now increase in euc isn't signif, 
#   but decreases in native and pasture are! 
# don't do log transform for P (in which nothing changes)


Krat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='K',],
                        na.action = na.omit,family='quasibinomial')
qqr(Krat100simp.pql) # ok? upper tail off
summary(Krat100simp.pql) # ratio starts bigger in noneuc (p=.053)
#   and decreases (p=.027)
# without It and JP2 (i.e. native Cerrado), no significant changes
# Native (just AF) almost higher with 3 separate LUs (p=.097)

Crat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='C',],
                        na.action = na.omit,family='quasibinomial')
qqr(Crat100simp.pql) # ok
summary(Crat100simp.pql) # no change

Nrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='N',],
                        na.action = na.omit,family='quasibinomial')
qqr(Nrat100simp.pql) # tails quite off; ok without Cerr nat
summary(Nrat100simp.pql) # with JP2 and It, decreases overall
# without, no change

Carat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='Ca',],
                        na.action = na.omit,family='quasibinomial')
qqr(Carat100simp.pql) 
summary(Carat100simp.pql) # marginal increase in euc, decrease in non
# starts higher in non
# without Cerr nat, same deal, but changes are signif only for pasture
# (marginally higher starting value p=.083 and decrease p=.067 in nat)

Prat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                         data=simp100[simp100$element=='P2',],
                         na.action = na.omit,family='quasibinomial')
qqr(Prat100simp.pql) 
summary(Prat100simp.pql) # shallower under non-euc (if using P not P2), no change

Alrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='Al',],
                        na.action = na.omit,family='quasibinomial')
qqr(Alrat100simp.pql) 
summary(Alrat100simp.pql) # gets deeper in euc but not in native, shallower in pasture
