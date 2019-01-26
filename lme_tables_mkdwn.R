source('soil_data_reader.R')
eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme) 

eucN20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) 
qqr(eucN20bm.lme) 

eucP20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
summary(eucP20bm.lme) 

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

eucP100bm.lme=lme(stock100~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2',],
                  random=~1|site/stand,na.action=na.omit)
summary(eucP100bm.lme) 
