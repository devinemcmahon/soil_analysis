source('soil_data_reader.R')

# Year trends between stocks in all eucalyptus stands together, to 20 cm
eucC20.lme=lme(log(stock20)~year,data=euc2deps[euc2deps$element=='C',],
               random=~1|site/stand,na.action=na.omit)
intervals(eucC20.lme)
eucN20.lme=lme(log(stock20)~year,random=~1|site/stand,
               data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
intervals(eucC20.lme)

eucP20.lme=lme(log(stock20)~year,random=~1|site/stand,
               data=euc2deps[euc2deps$element=='P2',],na.action = na.omit)
intervals(eucP20.lme)
eucK20.lme=lme(log(stock20)~year,data=euc2deps[euc2deps$element=='K' &
                                                 euc2deps$site!='Bp',],
               random=~1|site/stand,na.action=na.omit)
intervals(eucK20.lme)
eucCa20.lme=lme(log(stock20)~year,random=~1|site/stand,
                data=euc2deps[euc2deps$element=='Ca2',],na.action = na.omit)
intervals(eucK20.lme)

# Year trends broken out by biome
# Stock in 2004 in Atlantic Forest = exp(intercept +/- se)
# stock in 2016 in AF = exp(intercept)*exp(year coefficient) [both +/- respective SEs]
# stock in 2004 in Cerrado = exp(intercept)*exp(biome coef)
# stock in 2016 in Cerrado = exp(intercept)*exp(biome coef)*exp(year coef)*exp(intrxn coef)

# Increase in AF = (exp(year coef)-1)*100%
# Increase in Cerrado = (exp(year coef)*exp(intrxn coef)-1)*100%

eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
intervals(eucC20bm.lme)

eucN20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
intervals(eucN20bm.lme)

eucP20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
intervals(eucP20bm.lme)

eucK20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='K'&
                                 euc2deps$site!='Bp',],na.action = na.omit)
intervals(eucK20bm.lme)

eucCa20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Ca2',],
                  na.action = na.omit)
intervals(eucCa20bm.lme)

eucCN20bm.lme=lme(conc20~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='CN',],na.action = na.omit)
intervals(eucCN20bm.lme)


# 0-100 cm

eucC100bm.lme=lme(log(stock100)~year*biome,data=euc2deps2[euc2deps2$element=='C',],
                  random=~1|site/stand,na.action=na.omit)
intervals(eucC100bm.lme) 

eucN100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='N',],
                  random=~1|site/stand,na.action=na.omit)
intervals(eucN100bm.lme) # increases in AF, doesn't change in Cerrado

eucP100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2',],
                  random=~1|site/stand,na.action=na.omit)
intervals(eucP100bm.lme) 

eucK100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='K'&euc2deps2$site!='Bp',],
                  random=~1|site/stand,na.action=na.omit)
intervals(eucK100bm.lme) 

eucCa100bm.lme=lme(log(stock100)~year*biome,
                   data=euc2deps2[euc2deps2$element=='Ca2',],
                   random=~1|site/stand,na.action=na.omit)
intervals(eucCa100bm.lme) 

eucCN100bm.lme=lme(conc100~year*biome,random=~1|site/stand,
                   data=euc2deps2[euc2deps2$element=='CN',],na.action = na.omit)
intervals(eucCN100bm.lme) 


# Vegetation type interactions
# Not enough sites to do year*land use*biome
# Allow different year effects in different sites

Csimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='C',], 
               na.action=na.omit)
intervals(Csimp.lme4) # no effect of year on euc (p=.09) 

Nsimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
intervals(Nsimp.lme4) # no change; started with more N in native

nat2deps=droplevels(dats2deps[dats2deps$LU=='N',])
natN20bm.lme=lme(log(stock20)~year*biome,random=~1|stand,
                 data=nat2deps[nat2deps$element=='N',],na.action = na.omit)
qqr(natN20bm.lme) 
intervals(natN20bm.lme) # marginal decrease in AF, def increase in cerrado

Psimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='P2',], na.action=na.omit)
qqr(Psimp.lme4) # worse
intervals(Psimp.lme4) # no signif changes as you might expect

Ksimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='K',], na.action=na.omit)
intervals(Ksimp.lme4) 

Casimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
intervals(Casimp.lme4) 

Zrsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Zr',], na.action=na.omit)
intervals(Zrsimp.lme4) 

Nbsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Nb',], na.action=na.omit,
                control = lmeControl(opt = 'optim'))
intervals(Nbsimp.lme4) 


# 0-100 cm
# Removed pairs with known issues below 20 cm
simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]


Csimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='C',], na.action=na.omit)
intervals(Csimp100.lme) 

Nsimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='N',], na.action=na.omit)
intervals(Nsimp100.lme)

Psimp100.lme=lme(stock100~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='P2',], na.action=na.omit)
intervals(Psimp100.lme) # don't do log transform for P (in which nothing changes)

Ksimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                  data=simp100[simp100$element=='K',], na.action=na.omit)
intervals(Ksimp100.lme) # native veg mean change comes out negative
# but when you take into account the site effects, positive in both AF reserves
# when native = default factor level, not significant

Casimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                  data=simp100[simp100$element=='Ca2',], na.action=na.omit)
intervals(Casimp100.lme) # increase in E, N and P NOT diff with diff slopes 
# year term for N or P not signif

