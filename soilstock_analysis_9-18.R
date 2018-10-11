# Soil stock analysis

#setwd('C:\\Users\\Devin\\Documents\\Soil data')
source('soil_data_reader.R')

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

gen_ttable(ttests[ttests$element=='N',]$depth,
           ttests[ttests$element=='N',]$stand,
           ttests[ttests$element=='N',]$pval,
           ttests[ttests$element=='N',]$tstat,0.05)


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
# 20 cm: increases only in natveg (signif intrxn);
#   within site not stand, both natveg and interaction signif
# 100 cm within site, only intrxn signif again
# with Eu, now increase is significant and larger in N 
# should these be proportional increases?
#  now still signif increase, but decrease in P, no change in N,
#   with or without It.E1 and N
qqnorm(resid(allN20LU.lme))
qqline(resid(allN20LU.lme)) # nice
qqnorm(resid(allN100LU.lme))
qqline(resid(allN100LU.lme)) # tails way off, one crazy outlier
# log isn't better
plot(resid(allN20LU.lme)~dats2deps$year[dats2deps$element=='N' &
                                          !is.na(dats2deps$stock20)]) #ok
plot(resid(allN20LU.lme)~dats2deps$LU[dats2deps$element=='N' &
                                        !is.na(dats2deps$stock20)]) #eh
plot(resid(allN100LU.lme)~dats2deps$year[dats2deps$element=='N' &
                                          !is.na(dats2deps$stock100)]) 
# different means, similar spreads
plot(resid(allN100LU.lme)~dats2deps$LU[dats2deps$element=='N' &
                                        !is.na(dats2deps$stock100)]) 
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
                  data=dats2deps[dats2deps$element=='K' &
                                   dats2deps$site!='Bp',],na.action = na.omit)
summary(allK100LU.lme) # with Bp: decreases in 2016, 
# but nearly increases in native veg in 2016
# without: increases overall (fixed BD: no change), more (only) in native veg
# residuals not as bad as some others, but still a big outlier
# lower AIC with stand within site vs just stand

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
                  data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(Nrat100LU.lme) # residuals normalish
# proportion in top 20 cm decreases between years (p=.047); no LU effects

Crat100LU.lme=lme(stockratio~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='C',],na.action = na.omit)
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
Krat100euc2.lme=lme(stockratio~year,random=~1|site/stand,
                   data=euc2deps[euc2deps$element=='K'&
                                   euc2deps$site!='Bp',],na.action = na.omit)

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
                 data=dats2deps[dats2deps$element=='N',],na.action = na.omit)
summary(allN20LU.lme) # no change with fixed BD
qqr(allN20LU.lme) # pretty good


allCa20LU.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                   data=dats2deps[dats2deps$element=='Ca2',],
                  na.action = na.omit)
qqr(allCa20LU.lme) # really normal with log
summary(allCa20LU.lme)
# significant increase with time in euc/ overall 
# smaller increase in pasture (which was more to start at p=.05) at p=.08
# decrease in native 
allCa20LU2.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                   data=dats2deps[dats2deps$element=='Ca2' &
                        dats2deps$site!='Eu',],na.action = na.omit)
summary(allCa20LU2.lme) # same deal, but lower p-values

allP20LU.lme=lme(stock20~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='P2',],
                  na.action = na.omit)
qqr(allP20LU.lme) # no good with or without log

allK20LU.lme=lme(log(stock20)~year*LU,random=~1|site/stand,
                  data=dats2deps[dats2deps$element=='K' &
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
library(MASS)
Krat100LU2.pql=glmmPQL(stockratio~year*LU,random=~1|site/stand,
                       data=test2deps[test2deps$element=='K' &
                                        test2deps$site!='Bp',],
                       na.action = na.omit,family='quasibinomial')
summary(Krat100LU2.pql) # same deal as lme
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
# Also missing rep 1 of 2004 for 40-60 and rep 4 for 60-100
# reps 1 and 4 have way more C than 2 and 3, which are close to 2004 values

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
mean(shortE$BD100_16[shortE$element=='C'])
sd(shortE$BD100_16[shortE$element=='C'])
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
stkchgs=group_by(droplevels(shorterstk),stand,element)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04),
            budget=Budget/1000)
stkchgs2=stkchgs[stkchgs$stand!='It.E1',]

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
     xlab='Fertilizer - harvest, Mg ha-1',
     ylab='Observed change in stocks to 20 cm',las=1)
rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))

plot(chg20~budget,data=stkchgs,type='n', 
     xlab='Fertilizer - harvest, Mg ha-1',
     ylab='Observed change in stocks to 20 cm',las=1,
     xlim=c(-.2,.5),ylim=c(-.2,.5))
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))
legend('bottomright',pch=15,col=as.factor(levels(stkchgs$stand)),
       legend=levels(stkchgs$stand),bty='n',ncol=2)

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

