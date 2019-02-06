#source('soil_data_reader.R')

# Year trends between stocks in all eucalyptus stands together, to 20 cm
eucC20.lme=lme(log(stock20)~year,data=euc2deps[euc2deps$element=='C',],
               random=~1|site/stand,na.action=na.omit)
summary(eucC20.lme) 
qqr(eucC20.lme) # upper tail still off
summary(euc2deps$stock20[euc2deps$element=='C'&euc2deps$year=='04'])
summary(euc2deps$stock20[euc2deps$element=='C'&euc2deps$year=='16'])
exp(3.7818) # estimate for 2004
exp(3.7818)*exp(.1014) # estimate for 2016
# these are slightly lower than just taking the mean

eucN20.lme=lme(log(stock20)~year,random=~1|site/stand,
               data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
qqr(eucN20.lme) 
summary(eucN20.lme)
summary(euc2deps$stock20[euc2deps$element=='N'&euc2deps$year=='04'])

eucP20.lme=lme(log(stock20)~year,random=~1|site/stand,
               data=euc2deps[euc2deps$element=='P2',],na.action = na.omit)
qqr(eucP20.lme) # tails off with or without log
summary(eucP20.lme)
eucK20.lme=lme(log(stock20)~year,data=euc2deps[euc2deps$element=='K' &
                                                 euc2deps$site!='Bp',],
               random=~1|site/stand,na.action=na.omit)
summary(eucK20.lme) 
qqr(eucK20.lme) 
eucCa20.lme=lme(log(stock20)~year,random=~1|site/stand,
                data=euc2deps[euc2deps$element=='Ca2',],na.action = na.omit)
qqr(eucCa20.lme) # tails a bit off
summary(eucCa20.lme)

# Year trends broken out by biome
# Stock in 2004 in Atlantic Forest = exp(intercept +/- se)
# stock in 2016 in AF = exp(intercept)*exp(year coefficient) [both +/- respective SEs]
# stock in 2004 in Cerrado = exp(intercept)*exp(biome coef)
# stock in 2016 in Cerrado = exp(intercept)*exp(biome coef)*exp(year coef)*exp(intrxn coef)

# Increase in AF = (exp(year coef)-1)*100%
# Increase in Cerrado = (exp(year coef)*exp(intrxn coef)-1)*100%

eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme) 
C20sum=summary(eucC20bm.lme)
C20sum$coefficients
intervals(eucC20bm.lme)


eucC20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme2) 
C20sum2=summary(eucC20bm.lme2)
C20sum2$coefficients


eucN20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) 
qqr(eucN20bm.lme) 

eucP20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
summary(eucP20bm.lme)  # keep log to be consistent among nutrients?
qqr(eucP20bm.lme)

eucCa20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Ca2',],
                  na.action = na.omit)
summary(eucCa20bm.lme) 

eucK20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='K'&
                                 euc2deps$site!='Bp',],na.action = na.omit)
summary(eucK20bm.lme) 

eucCN20bm.lme=lme(conc20~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='CN',],na.action = na.omit)
summary(eucCN20bm.lme) 

lmts=c('C','N','P','K','Ca')
lmts2=c('C','N','P','K','Ca','CN')
allT20s=data.frame(rbind(round(summary(eucC20bm.lme)$tTable,3),
                       round(summary(eucN20bm.lme)$tTable,3),
                       round(summary(eucP20bm.lme)$tTable,3),
                       round(summary(eucK20bm.lme)$tTable,3),
                       round(summary(eucCa20bm.lme)$tTable,3)))
allT20s$Element=rep(lmts,each=4)
allT20s$Term=rep(c('Intercept','Year_2016','Biome_Cerrado','2016_x_Cerrado'),5)
allT20s=allT20s[,c(6,7,1,2,3,4,5)]
allT20s
write.csv(allT20s,'allT20s.csv',row.names = F)


# 0-100 cm

eucCN100bm.lme=lme(conc100~year*biome,random=~1|site/stand,
                   data=euc2deps2[euc2deps2$element=='CN',],na.action = na.omit)
summary(eucCN100bm.lme) 

eucK100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='K'&euc2deps2$site!='Bp',],
                  random=~1|site/stand,na.action=na.omit)
summary(eucK100bm.lme) 

eucC100bm.lme=lme(log(stock100)~year*biome,data=euc2deps2[euc2deps2$element=='C',],
                  random=~1|site/stand,na.action=na.omit)
summary(eucC100bm.lme) 

eucCa100bm.lme=lme(log(stock100)~year*biome,
                   data=euc2deps2[euc2deps2$element=='Ca2',],
                   random=~1|site/stand,na.action=na.omit)
summary(eucCa100bm.lme) 

eucN100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='N',],
                  random=~1|site/stand,na.action=na.omit)
summary(eucN100bm.lme) # increases in AF, doesn't change in Cerrado
qqr(eucN100bm.lme)

eucP100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2',],
                  random=~1|site/stand,na.action=na.omit)
summary(eucP100bm.lme) 
qqr(eucP100bm.lme)

allT100s=data.frame(rbind(round(summary(eucC100bm.lme)$tTable,3),
                         round(summary(eucN100bm.lme)$tTable,3),
                         round(summary(eucP100bm.lme)$tTable,3),
                         round(summary(eucK100bm.lme)$tTable,3),
                         round(summary(eucCa100bm.lme)$tTable,3)))
allT100s$Element=rep(lmts,each=4)
allT100s$Term=rep(c('Intercept','Year_2016','Biome_Cerrado','2016_x_Cerrado'),5)
allT100s=allT100s[,c(6,7,1,2,3,4,5)]
allT100s
write.csv(allT100s,'allT100s.csv',row.names = F)


# ratios of stock in top 20 cm to stock in top 100 cm
Krateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='K',],# &
                     #                 euc2deps2$site!='Bp',],
                     na.action = na.omit,family='binomial')
summary(Krateuc2.pql) 

Nrateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='N',],
                     na.action = na.omit,family='binomial')
summary(Nrateuc2.pql) 
qqr(Nrateuc2.pql) #decreases in AF, no change in Cerrado

Crateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='C',],
                     na.action = na.omit,family='binomial')
summary(Crateuc2.pql) # nada
qqr(Crateuc2.pql) # slightly better than with no biome term

Carateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                      data=euc2deps2[euc2deps2$element=='Ca2',],
                      na.action = na.omit,family='binomial')
summary(Carateuc2.pql) # increases a bunch in Cerrado only

Prateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='P2',],
                     na.action = na.omit,family='binomial')
summary(Prateuc2.pql) # no significant terms

allTrats=data.frame(rbind(round(summary(Crateuc2.pql)$tTable,3),
                         round(summary(Nrateuc2.pql)$tTable,3),
                         round(summary(Prateuc2.pql)$tTable,3),
                         round(summary(Krateuc2.pql)$tTable,3),
                         round(summary(Carateuc2.pql)$tTable,3)))
allTrats$Element=rep(lmts,each=4)
allTrats$Term=rep(c('Intercept','Year_2016','Biome_Cerrado','2016_x_Cerrado'),5)
allTrats=allTrats[,c(6,7,1,2,3,4,5)]
allTrats
write.csv(allTrats,'allTrats.csv',row.names = F)

# Vegetation type interactions
# Not enough sites to do year*land use*biome
# Allow different year effects in different sites
Casimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
qqr(Casimp.lme4) # actually better than lme3 (without random slope)
summary(Casimp.lme4) # increase in euc, not signif in pasture

Csimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='C',], 
               na.action=na.omit)
qqr(Csimp.lme4) # tails off,esp lower
summary(Csimp.lme4) # no effect of year on euc (p=.09) 
Nsimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
qqr(Nsimp.lme4) 
summary(Nsimp.lme4) # no change; started with more N in native

nat2deps=droplevels(dats2deps[dats2deps$LU=='N',])
natN20bm.lme=lme(log(stock20)~year*biome,random=~1|stand,
                 data=nat2deps[nat2deps$element=='N',],na.action = na.omit)
qqr(natN20bm.lme) 
summary(natN20bm.lme) # marginal decrease in AF, def increase in cerrado

Ksimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='K',], na.action=na.omit)
summary(Ksimp.lme4) 

Psimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='P2',], na.action=na.omit)
qqr(Psimp.lme4) # worse
summary(Psimp.lme4) # no signif changes as you might expect

Zrsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Zr',], na.action=na.omit)
summary(Zrsimp.lme4) # significant year-native intrxn

Nbsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Nb',], na.action=na.omit,
                control = lmeControl(opt = 'optim'))
summary(Nbsimp.lme4) # no year or interaction effects

simpT20s=data.frame(rbind(round(summary(Csimp.lme4)$tTable,3),
                          round(summary(Nsimp.lme4)$tTable,3),
                          round(summary(Psimp.lme4)$tTable,3),
                          round(summary(Ksimp.lme4)$tTable,3),
                          round(summary(Casimp.lme4)$tTable,3)))
simpT20s$Element=rep(lmts,each=6)
simpT20s$Transform=rep(c('ln','none','none','ln','ln'),each=6)
simpT20s$Term=rep(c('Intercept','Year_2016','Native','Pasture',
                   '2016_x_Native','2016_x_Pasture'),5)
simpT20s=simpT20s[,c(6,7,8,1,2,3,4,5)]
simpT20s
write.csv(simpT20s,'simpT20s.csv',row.names = F)

# 0-100 cm
# Removed pairs with known issues below 20 cm
simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]


Csimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='C',], na.action=na.omit)
qqr(Csimp100.lme)
summary(Csimp100.lme) 
# with random ~1|year+site (more appropriate),
#   only signif thing is decrease in pasture (in native, p=.08)

Nsimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='N',], na.action=na.omit)
summary(Nsimp100.lme)

Ksimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                  data=simp100[simp100$element=='K',], na.action=na.omit)
summary(Ksimp100.lme) # native veg mean change comes out negative
# but when you take into account the site effects, positive in both AF reserves
# when native = default factor level, not significant
Ksimp100.lme$coefficients

Psimp100.lme=lme(stock100~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='P2',], na.action=na.omit)
summary(Psimp100.lme) # don't do log transform for P (in which nothing changes)

Casimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                  data=simp100[simp100$element=='Ca2',], na.action=na.omit)
summary(Casimp100.lme) # increase in E, N and P NOT diff with diff slopes 
# year term for N or P not signif

simpT100s=data.frame(rbind(round(summary(Csimp100.lme)$tTable,3),
                          round(summary(Nsimp100.lme)$tTable,3),
                          round(summary(Psimp100.lme)$tTable,3),
                          round(summary(Ksimp100.lme)$tTable,3),
                          round(summary(Casimp100.lme)$tTable,3)))
simpT100s$Element=rep(lmts,each=6)
simpT100s$Transform=rep(c('ln','ln','none','ln','ln'),each=6)
simpT100s$Term=rep(c('Intercept','Year_2016','Native','Pasture',
                    '2016_x_Native','2016_x_Pasture'),5)
simpT100s=simpT100s[,c(6,7,8,1,2,3,4,5)]
simpT100s
write.csv(simpT100s,'simpT100s.csv',row.names = F)
