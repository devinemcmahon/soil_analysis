# Soil stock analysis

#setwd('C:\\Users\\Devin\\Documents\\Soil data')
source('soil_data_reader.R')

# clay content 0-10 cm
clay = data.frame(stand=c('EuA','EuB','BO','Vg',
                          'BpA','BpB','ItA','ItB',
                          'JPA','JPB'),
                  pct5=c(15,18,58,74,85,88,77,65,15,16),
                  pct80=c(35,41,72,74,77,87,79,74,18,22))
summary(clay$pct5)
summary(clay$pct80)

# Figure S2
# Highlighting the scale of the spatial heterogeneity in Bps
xyplot(depth~value/1000|stand+year,groups=rep,type='p',ylab='Depth (cm)',
       data=dats[dats$element=='K' & dats$site=='Bp',],
       ylim=c(90,0),xlab='K (g / kg)',as.table=T)


# Figure S3
par(mfrow=c(2,3))
yrdiffstockplot20_bmLUall(tstock[tstock$element=='C',],fulllegend = F)
legend('bottomright',bty='n',legend='a',cex=1.5)
yrdiffstockplot20_bmLUall(tstock[tstock$element=='N',],fulllegend = T)
legend('bottomright',bty='n',legend='b',cex=1.5)
yrdiffstockplot20_bmLUall(tstock[tstock$element=='K' &
                                   tstock$site!='Bp',],fulllegend = F)
legend('bottomright',bty='n',legend='c',cex=1.5)

yrdiffstockplot20_bmLUall(tstock[tstock$element=='P2',],label=F,fulllegend = F)
legend('topleft',bty='n',legend='P (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='d',cex=1.5)

yrdiffstockplot20_bmLUall(tstock[tstock$element=='Ca2',],label=F,fulllegend = F)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='e',cex=1.5)

yrdiffstockplot20_bmLUall(tstock[tstock$element=='Ca2' & tstock$stock20_16<1,],label=F)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='f',cex=1.5)
par(mfrow=c(1,1))

sum(is.na(dats2deps$stock20[dats2deps$element=='P']))
sum(is.na(dats2deps$stock20[dats2deps$element=='P2']))
length(unique(widedats$ID[is.na(widedats$P)&!is.na(widedats$P_dl) &
                            widedats$stand!='JP.A'&
                            !is.element(widedats$site,c('Cr','TM'))]))
length(unique(widedats$ID[is.na(widedats$Ca)&!is.na(widedats$Ca_dl)&
                            widedats$stand!='JP.A'&
                            !is.element(widedats$site,c('Cr','TM'))]))
length(unique(widedats$ID[!is.na(widedats$P)&!is.na(widedats$P_dl)&
                            widedats$stand!='JP.A'&
                            !is.element(widedats$site,c('Cr','TM'))]))
length(unique(widedats$ID[is.na(widedats$Ca)&!is.na(widedats$P)&
                            widedats$stand!='JP.A'&
                            !is.element(widedats$site,c('Cr','TM'))]))
length(unique(widedats$ID[is.na(widedats$P)&!is.na(widedats$Ca)&
                            widedats$stand!='JP.A'&
                            !is.element(widedats$site,c('Cr','TM'))]))
table(widedats$stand[is.na(widedats$P)])
length(unique(widedats$P_dl[is.na(widedats$P)&widedats$stand!='JP.A'&
                     !is.element(widedats$site,c('Cr','TM'))])) # just 6
length(unique(widedats$Ca_dl[is.na(widedats$Ca)&widedats$stand!='JP.A'&
                              !is.element(widedats$site,c('Cr','TM'))])) # just 10
table(widedats$stand[is.na(widedats$P)],widedats$P_dl[is.na(widedats$P)])
# almost always .0003
table(widedats$stand[is.na(widedats$Ca)],widedats$Ca_dl[is.na(widedats$Ca)])
# usually .001 except in Bp.E1? Mostly an issue in Vg and Bp, some It.E2 and BO

# Year trends between stocks in eucalyptus stands
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


# Year trends between stocks in each biome (Figure 2)

# Added a random slope (year) to these models, 12-29-18
# Does that make sense? I do want to know the overall effect of year
#   which is just an offset, not a slope
# Leave it how it was?
# But allow different slopes for the different LU pairs, below?
# That doesn't seem right--the two analyses should be the same
# Adding a random year effect for each site is the more conservative 
#   approach, and probably more appropriate
# But some sites have just one stand
# Puts all the effect of that stand in the error term?
# Leave how it was for now
# Different nesting structure for biome analysis: 
#   Include year effect to allow for any changes within site across veg type
eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                    data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme) # increase in Cerrado only (also if using euc2deps2)
qqr(eucC20bm.lme) # again, log tranform helps a bit, increases significance
eucC20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme2) # increase in Cerrado now barely signif (p=.046)
qqr(eucC20bm.lme2) # slightly better?
anova(eucC20bm.lme,eucC20bm.lme2) # not different
intervals(eucC20bm.lme2) # near-perfect correlations between year and intercept
# for random effects--overfit?

lmesum=summary(eucC20bm.lme)
lmesum$tTable


# N: doesn't converge with a random slope (does with optim)
eucN20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) # no change
qqr(eucN20bm.lme) # mostly ok
intervals(eucN20bm.lme)

#eucN20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
#                  data=euc2deps[euc2deps$element=='N',],na.action = na.omit,
#                  control=lmeControl(opt='optim'))
# ambiguous error message with log transform
summary(eucN20bm.lme2) # no change
qqr(eucN20bm.lme2) # tails off without log
anova(eucN20bm.lme,eucN20bm.lme2)
# 2 has lower AIC but higher BIC; differ at p=.02 (no longer; 2 is worse)

eucP20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
qqr(eucP20bm.lme) # tails off without Eu, but much less bad
summary(eucP20bm.lme) # no change
intervals(eucP20bm.lme)
#eucP20bm.lme2=lme(stock20~year*biome,random=~1|site/stand,
#                 data=euc2deps[euc2deps$element=='P2'&
#                                 euc2deps$site!='Eu',],na.action = na.omit)
#qqr(eucP20bm.lme2) # tails off without Eu, but much less bad
#summary(eucP20bm.lme2) # when excluding Eu, increase (p=.01) due to Vg;
# p=.08 for opposite sign year-Cerrado interaction
# Eu not necessarily messier than other sites; keep all stands


eucCa20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Ca2',],
                  na.action = na.omit)#,control = lmeControl(opt='optim'))
qqr(eucCa20bm.lme)
summary(eucCa20bm.lme) # maybe increase in AF (p=.07), increase in Cer (.0005)
intervals(eucCa20bm.lme)
eucCa20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
                 data=euc2deps[euc2deps$element=='Ca2',],
                 na.action = na.omit,control = lmeControl(opt='optim'))
#qqr(eucCa20bm.lme2)
#summary(eucCa20bm.lme2) 
intervals(eucCa20bm.lme2) # year16 effect and intercept highly correlated
anova(eucCa20bm.lme,eucCa20bm.lme2)
# with a random slope, no signif effects
# but that seems wrong, clearly a lot of Ca was added in the Cerrado
# although not consistently across cerrado stands, hence loss of significance?
# Or is it just not meaningful to allow a different change between
#   the four observations in each site--not enough data to fit this?
# go with just a random intercept?
# I think that is correct, just random intercepts

eucK20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='K'&
                                  euc2deps$site!='Bp',],na.action = na.omit)
qqr(eucK20bm.lme)
summary(eucK20bm.lme) # increases in AF only
#euc2deps$biome=factor(euc2deps$biome,levels=c('Cer','AF')) #yep

#eucMg20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
#                  data=euc2deps[euc2deps$element=='Mg2',],na.action = na.omit)
#qqr(eucMg20bm.lme) # tails way off even with log
#summary(eucMg20bm.lme) # no change

#eucAl20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
#                  data=euc2deps[euc2deps$element=='Al',],na.action = na.omit)
#qqr(eucAl20bm.lme) # not very good; removing log worse
#summary(eucAl20bm.lme) # decreases (p=0.030), no biome effect
# Note that there is more Al in top 20 cm under eucalyptus than native veg
#   and most in pasture

eucCN20bm.lme=lme(conc20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='CN',],na.action = na.omit)
summary(eucCN20bm.lme) # no change in AF, increases in Cer
qqr(eucCN20bm.lme) # ok

# From t-tables in lme_tables_printing: magnitude of signif chg
# C in Cerrado:
exp(-.062)*exp(.272) # 23% increase
# K in AF
exp(.181) #20% 
# Ca in Cerrado:
exp(.512)*exp(1.34) # 530% increase
# Ca in AF


# or:
summary(stkchgs$chgrt20[stkchgs$element=='N']) # not all the stands
stkchgsall=group_by(shorttstk,stand,element,biome,LU)%>%
  summarise(chg100=stock100_16-stock100_04,stk100_16=stock100_16,
            chg20=stock20_16-stock20_04,stk20_16=stock20_16,
            sdchg20=sqrt(sd20_04^2+sd20_16^2)/2, 
            # = sqrt(sd04^2/n04+sd16^2/n16), assuming n=4
            sdchg100=sqrt(sd100_04^2+sd100_16^2)/2,
            stk20_04=stock20_04,
            chgrt100=(stock100_16-stock100_04)/stock100_04,
            chgrt20=(stock20_16-stock20_04)/stock20_04,
            chgln100=log(stock100_16/stock100_04),
            chgln20=log(stock20_16/stock20_04))
summary(stkchgsall$chgrt20[stkchgsall$LU=='E'&
                             stkchgsall$element=='C'&
                             stkchgsall$biome=='Cer'])
# lme estimates are more representative of change?
# If expressing change as a percent, use mean?
# But log transform is better
# should lme estimates be the basis for Figure 2, too?
# Don't think so because that's stocks
exp(summary(stkchgsall$chgln20[stkchgsall$LU=='E'&
                             stkchgsall$element=='C'&
                             stkchgsall$biome=='Cer']))
# this matches the lmes more closely, ok
exp(summary(stkchgsall$chgln20[stkchgsall$LU=='E'&
                                 stkchgsall$element=='Ca'&
                                 stkchgsall$biome=='Cer']))
# not equivalent, though
# use same values as in text--pick one and stick with it
# maybe use the lme outputs because they match the stats
summary(stkchgsall$stk20_16[stkchgsall$LU=='E'&
                              stkchgsall$element=='C'&
                              stkchgsall$biome=='Cer'])
exp(4.017)*exp(-.062)*exp(-.44)*exp(.272) #(lme est)
# not quite the median or mean, but close
# oh, they're different because of when the avg is taken
# diff no. of reps per stand, diff outlier weighting
# ok.
# Ca in AF to 100 cm:
exp(.266) #wait
exp(summary(stkchgsall$chgln100[stkchgsall$LU=='P'&
                              stkchgsall$element=='C']))
# 16% decrease
summary(Csimp100.lme)
exp(-.0847)*exp(-.1394) #20% decrease
# mean 16/mean 04:
summary(simp100$stock100[simp100$element=='C'&simp100$LU=='P'&
                           simp100$year=='16'])/
  summary(simp100$stock100[simp100$element=='C'& simp100$LU=='P'&
                            simp100$year=='04'])
# change in mean is minus 21%

# 0-100 cm
eucCN100bm.lme=lme(conc100~year*biome,random=~1|site/stand,
                  data=euc2deps2[euc2deps2$element=='CN',],na.action = na.omit)
summary(eucCN100bm.lme) # decreases in AF, increases in Cer
qqr(eucCN100bm.lme) # one point off at each tail, upper far off
plot(eucCN100bm.lme)
mean(euc2deps2$conc100[euc2deps2$element=='CN'&
                         euc2deps2$year=='04'&
                         euc2deps2$biome=='AF'],na.rm=T)
# 16.3
# lme estimate: equivalent to mean
# ok. makes sense to call the lme estimates means
# they're just means of all the replicates, not first avgd by stand
# and I guess they're really medians for ln convert-backs?
median(euc2deps2$conc100[euc2deps2$element=='CN'&
                         euc2deps2$year=='04'&
                         euc2deps2$biome=='AF'],na.rm=T)
# 16.8
mean(euc2deps2$conc100[euc2deps2$element=='CN'&
                         euc2deps2$year=='16'&
                         euc2deps2$biome=='AF'],na.rm=T)
mean(euc2deps2$conc100[euc2deps2$element=='CN'&
                         euc2deps2$year=='04'&
                         euc2deps2$biome=='Cer'],na.rm=T)
mean(euc2deps2$conc100[euc2deps2$element=='CN'&
                         euc2deps2$year=='16'&
                         euc2deps2$biome=='Cer'],na.rm=T)



eucK100bm.lme=lme(log(stock100)~year*biome,
                data=euc2deps2[euc2deps2$element=='K'&euc2deps2$site!='Bp',],
                 random=~1|site/stand,na.action=na.omit)
qqr(eucK100bm.lme)
summary(eucK100bm.lme) # increases in AF, maybe decreases in Cerrado (no)
#euc2deps2$biome=factor(euc2deps2$biome,levels=c('Cer','AF')) 

eucC100bm.lme=lme(log(stock100)~year*biome,data=euc2deps2[euc2deps2$element=='C',],
                random=~1|site/stand,na.action=na.omit)
qqr(eucC100bm.lme)
summary(eucC100bm.lme) # increases in Cerrado only (p=.005)

eucCa100bm.lme=lme(log(stock100)~year*biome,
                   data=euc2deps2[euc2deps2$element=='Ca2',],
                 random=~1|site/stand,na.action=na.omit)
summary(eucCa100bm.lme) # increases at p < 0.0001 in Cerrado only
qqr(eucCa100bm.lme) 

eucN100bm.lme=lme(log(stock100)~year*biome,
                   data=euc2deps2[euc2deps2$element=='N',],
                   random=~1|site/stand,na.action=na.omit)
summary(eucN100bm.lme) # increases in AF, doesn't change in Cerrado
qqr(eucN100bm.lme)

eucP100bm.lme=lme(stock100~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2',],
                  random=~1|site/stand,na.action=na.omit)
qqr(eucP100bm.lme) # quite bad
summary(eucP100bm.lme) # now significant! increases in AF
eucP100bm.lme2=lme(stock100~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2'&
                                   euc2deps2$site!='Eu',],
                  random=~1|site/stand,na.action=na.omit)
qqr(eucP100bm.lme2) # one outlier each end; almost signif in AF w/o Eu
summary(eucP100bm.lme2)

# N change to 100 cm in AF
exp(summary(stkchgsall$chgln100[stkchgsall$LU=='E'&
                               stkchgsall$element=='N'&
                               stkchgsall$biome=='AF']))
# +12% median, +18% mean
exp(.161) #+17%
exp(.161-.046) #Minus 1 standard error = +12% 
exp(.161+.046) #+23%
# Ca in Cer
exp(.266)*exp(1.17) # +320%



# significant changes (year effect) at 20 cm: C in Cerrado +
# Ca in AF (-) p=.075, Ca in Cer + (p=.0005)
# C:N increases in Cer, ok
# K increases in AF (excluding Bp)
# no change in N or P

# significant at 100 cm: +N in AF, +C in Cerrado, 
#   CN decreases in AF and increases in Cerr accordingly
# +Ca in Cerr, +K AF

# overall: +C, Ca at both depths in Cerrado
# + K at both depths in AF
# + N in 100 cm (not signif to just 20) in AF, not Cerrado

#Mean stocks (for fig 2)
table(shorttstk$stockunit,shorttstk$element)
nutstk=group_by(dats2deps[dats2deps$element %in% c('C','N','K','P2','Ca2') &
                            dats2deps$stand !='It.E1',],
                element,LU,biome,year)%>%
  mutate(mnrat=mean(stockratio,na.rm=T),
         #serat=sd(stockratio,na.rm=T)/sqrt(sum(!is.na(stockratio))-1))
         serat=sd(stockratio,na.rm=T)/sqrt(sum(!is.na(stockratio))))
nutstk$stock100[nutstk$site=='Bp'&nutstk$element=='K']=NA
nutstk$stock20[nutstk$site=='Bp'&nutstk$element=='K']=NA
nutstk=group_by(nutstk,element,LU,biome,year,mnrat,serat)%>%
  summarise(stk20down=mean(stock100-stock20,na.rm=T),
            stk20up=mean(stock20,na.rm=T),
            stk100=mean(stock100,na.rm=T),
            nobs=n(),
            se100=sd(stock100,na.rm=T)/sqrt(nobs),
            se20=sd(stock20,na.rm=T)/sqrt(nobs))
#library(reshape2)
nutstk=ungroup(nutstk)
#nutmelt=melt.data.frame(nutstk,measure.vars=c("stk20down","stk20up"))
# "Names do not match previous names"--why not??
#nutmelt=melt.data.frame(nutstk,measure.vars=5:6) #likewise
# Fine I'll do it by hand
# oh, this is probably because I loaded both reshape and reshape2
nutstkd=mutate(nutstk,depth=rep('20down'))
nutstkd=nutstkd[,-which(names(nutstkd) %in% c('stk20up','se20'))]
nutstku=mutate(nutstk,depth=rep('20up'))
nutstku=nutstku[,-which(names(nutstku) %in% c('stk20down','se100'))]
names(nutstkd)[which(names(nutstkd)=='stk20down')]='stock'
names(nutstkd)[which(names(nutstkd)=='se100')]='se'
names(nutstku)[which(names(nutstku)=='stk20up')]='stock'
names(nutstku)[which(names(nutstku)=='se20')]='se'
nutstk2=rbind(nutstku,nutstkd)

nutstk2$element=as.character(nutstk2$element)
nutstk2$element[nutstk2$element=='Ca2']='Ca'
nutstk2$element[nutstk2$element=='P2']='P'
nutstk2$element=factor(nutstk2$element,levels=c('C','N','K','P','Ca'))
nutstk2$year=factor(nutstk2$year,levels=c('16','04'))
nutstk2$stock[nutstk2$element=='C']=nutstk2$stock[nutstk2$element=='C']/10
nutstk2$se[nutstk2$element=='C']=nutstk2$se[nutstk2$element=='C']/10
nutstk2$stk100[nutstk2$element=='C']=nutstk2$stk100[nutstk2$element=='C']/10
nutstk2=mutate(nutstk2,hibar=ifelse(depth=='20up',stock+se,stk100+se),
               lobar=ifelse(depth=='20up',stock-se,stk100-se))
#nutstk2=group_by(nutstk2,element,biome)%>%
#               mutate(maxhibar=max(hibar))
nutstk2=group_by(nutstk2,biome)%>%mutate(maxhibar=max(hibar))
nutstkE=nutstk2[nutstk2$LU=='E',]
#plot(stock~element,data=nutstkE)
#library(ggplot2)
nutstkE <- with(nutstkE, nutstkE[order(element,year,depth),])
nutstkE=mutate(nutstkE,sig=rep(NA),sigrat=rep(''))
nutstkE$sig[nutstkE$depth=='20down'&nutstkE$element=='N'&nutstkE$biome=='AF'&
              nutstkE$year=='16']=1
nutstkE$sig[nutstkE$element=='Ca'& nutstkE$biome=='Cer'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig[nutstkE$element=='C'&nutstkE$biome=='Cer'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig[nutstkE$element=='K'&nutstkE$biome=='AF'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig=as.integer(nutstkE$sig)
nutstkE$sigrat[nutstkE$depth=='20down'&nutstkE$element=='N'&
                 nutstkE$biome=='AF'&
                 nutstkE$year=='16']='*'
nutstkE$sigrat[nutstkE$depth=='20down'&nutstkE$element=='K'&
                 nutstkE$biome=='AF'&
                 nutstkE$year=='16']='*'
nutstkE$sigrat[nutstkE$depth=='20down'&nutstkE$element=='Ca'&
                 nutstkE$biome=='Cer'&
                 nutstkE$year=='16']='*'
# Figure 2

labls <- c(AF = "Atlantic Forest", Cer = "Cerrado")

ggplot(data=nutstkE, aes(x=year, y=stock, fill=depth)) +
  geom_bar(stat="identity") + 
  facet_grid(element~biome,labeller = labeller(biome=labls)) + 
  coord_flip() +
  labs(y="Stock (10 Mg/ha for C; Mg/ha for other elements)",
       x="Year", fill="Depth") +
  theme(strip.text.y = element_text(angle = 0),
        legend.position=c(0.75,0.14),
        panel.background = element_rect(fill='white'),
        panel.grid.major.x = element_line(colour='grey80'),
        panel.grid.major.y = element_blank()) +
  geom_errorbar(aes(ymax=hibar,  ymin=lobar), width=0.15) +
  scale_fill_manual(values=c('slategray','lightblue'),
                      labels=c('20-100 cm','0-20 cm'),
                      guide = guide_legend(reverse=TRUE,title=NULL))+
  geom_point(mapping = aes(y = (hibar+1)*(sig>0)),
             shape=18,size=3,show.legend=F)+
  geom_text(mapping = aes(y = maxhibar-7,
                          label = paste(round((mnrat*100),0),' (',
                                        round((serat*100),0),')',
                                        sigrat,sep='')), 
            nudge_x = .3,data=nutstkE[nutstkE$depth=='20down',],
            na.rm=T,parse=F,show.legend=F,size=3.5,hjust=0)



tapply(shorttstk$stock100_16,shorttstk$element,
       function(x){mean(x,na.rm=T)})
tapply(shorttstk$stock100_04,shorttstk$element,
       function(x){mean(x,na.rm=T)})

tapply(shorttstk$stock100_16,shorttstk$element,
       function(x){sd(x,na.rm=T)})
mean(shorttstk$BD100_16[shorttstk$element=='C'])
sd(shorttstk$BD100_16[shorttstk$element=='C'])

tapply(shorttstk$stock20_16,shorttstk$element,
       function(x){mean(x,na.rm=T)})


shortE = shorttstk[shorttstk$LU=='E' &
                     shorttstk$stand!='It.E1',]
tapply(shortE$stock100_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock100_16,shortE$element,
       function(x){sd(x,na.rm=T)})
tapply(shortE$stock20_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock20_16,shortE$element,
       function(x){sd(x,na.rm=T)})

tapply(shorttstk$stock20_16[shorttstk$LU=='E'],
       shorttstk$element[shorttstk$LU=='E'],
       function(x){mean(x,na.rm=T)})
tapply(shorttstk$stock20_04[shorttstk$LU=='E'],
       shorttstk$element[shorttstk$LU=='E'],
       function(x){mean(x,na.rm=T)})

sefun=function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x))-1)}
tapply(shorttstk$stock20_16[shorttstk$LU=='E'],
       shorttstk$element[shorttstk$LU=='E'],sefun)
tapply(shorttstk$stock20_04[shorttstk$LU=='E'],
       shorttstk$element[shorttstk$LU=='E'],sefun)

tapply(shortE$rat_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$rat_16,shortE$element,
       function(x){sd(x,na.rm=T)})
tapply(shortE$rat_16,shortE$element,sefun)
mean(shortE$BD100_16[shortE$element=='C'])
sd(shortE$BD100_16[shortE$element=='C'])

tapply(shortE$stock100_16[shortE$element=='C'],
       shortE$biome[shortE$element=='C'],
       function(x){mean(x,na.rm=T)}) 

mean(shorttstk$stock20_04[shorttstk$LU=='E' & shorttstk$element=='C']) # 46.0
mean(shorttstk$stock20_16[shorttstk$LU=='E' & shorttstk$element=='C']) # 50.1
mean(shorttstk$stock20_04[shorttstk$LU=='E' & shorttstk$element=='K'&
                            shorttstk$site!='Bp']) 
mean(shorttstk$stock20_16[shorttstk$LU=='E' & shorttstk$element=='K'&
                            shorttstk$site!='Bp']) 
sefun(shorttstk$stock20_04[shorttstk$LU=='E' & shorttstk$element=='K'&
                            shorttstk$site!='Bp']) 
sefun(shorttstk$stock20_16[shorttstk$LU=='E' & shorttstk$element=='K'&
                            shorttstk$site!='Bp']) 

tapply(shortE$stock100_16[shortE$biome=='AF'],
       shortE$element[shortE$biome=='AF'],
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock100_16[shortE$biome=='AF'],
       shortE$element[shortE$biome=='AF'],sefun)

tapply(shortE$stock100_04[shortE$biome=='Cer'],
       shortE$element[shortE$biome=='Cer'],
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock100_04[shortE$biome=='Cer'],
       shortE$element[shortE$biome=='Cer'],sefun)
tapply(shortE$stock100_16[shortE$biome=='Cer'],
       shortE$element[shortE$biome=='Cer'],
       function(x){mean(x,na.rm=T)})
tapply(shortE$stock100_16[shortE$biome=='Cer'],
       shortE$element[shortE$biome=='Cer'],sefun)


# just-euc ratios
# These are proportion data and don't come from normal distributions
# Try generalized mixed models
#library(MASS)
euc2deps2$biome=factor(euc2deps2$biome,levels=c('Cer','AF'))
euc2deps2$biome=factor(euc2deps2$biome,levels=c('AF','Cer'))

Krateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='K' ,],
                     na.action = na.omit,family='binomial')
# using quasibinomial = same result, but I don't understand it as well
summary(Krateuc2.pql) # same deal as lme, increases at p=.07
qqr(Krateuc2.pql)
# with Bp included, same result, somewhat smaller increase

Nrateuc2.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='N',],
                     na.action = na.omit,family='binomial')
summary(Nrateuc2.pql) # again, same as lme: decreases, p=.026
qqr(Nrateuc2.pql) #ok
Nrateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='N',],
                     na.action = na.omit,family='binomial')
summary(Nrateuc2.pql) 
qqr(Nrateuc2.pql) #decreases in AF, no change in Cerrado

Crateuc.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                    data=euc2deps2[euc2deps2$element=='C',],
                    na.action = na.omit,family='binomial')
summary(Crateuc.pql) # nada
qqr(Crateuc.pql) # not great
Crateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                    data=euc2deps2[euc2deps2$element=='C',],
                    na.action = na.omit,family='binomial')
summary(Crateuc2.pql) # nada
qqr(Crateuc2.pql) # slightly better
Crateuc2.lme=lme(log(stockratio)~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='C',],
                     na.action = na.omit)
summary(Crateuc2.lme) # nada
qqr(Crateuc2.lme) # slightly better than pql

Carateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='Ca2',],
                     na.action = na.omit,family='binomial')
summary(Carateuc2.pql) # increases a bunch in Cerrado only
qqr(Carateuc2.pql) # ok
Carateuc2.lme=lme(log(stockratio)~year*biome,random=~1|site/stand,
                  data=euc2deps2[euc2deps2$element=='Ca2',],
                  na.action = na.omit)
summary(Carateuc2.lme) # similar
qqr(Carateuc2.lme) # worse than pql

Carateuc3.pql=glmmPQL(stockratio~year*biome,random=~1+year|site/stand,
                      data=euc2deps2[euc2deps2$element=='Ca2',],
                      na.action = na.omit,family='binomial',
                      control=lmeControl(opt = 'optim'))
summary(Carateuc3.pql) # increases a bunch in Cerrado, lower p-value
qqr(Carateuc3.pql) # ok except 1 outlier
Carateuc3.lme=lme(stockratio~year*biome,random=~1+year|site/stand,
                 data=euc2deps2[euc2deps2$element=='Ca2',],
                 na.action = na.omit,
                 control=lmeControl(opt = 'optim'))
qqr(Carateuc3.lme) # bad with log, doesn't converge without


Prateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='P2',],
                     na.action = na.omit,family='binomial')
summary(Prateuc2.pql) # no significant terms
qqr(Prateuc2.pql) # messed up as usual
Prateuc3.pql=glmmPQL(stockratio~year*biome,random=~1+year|site/stand,
                      data=euc2deps2[euc2deps2$element=='P2',],
                      na.action = na.omit,family='binomial',
                      control=lmeControl(opt = 'optim'))
summary(Prateuc3.pql) 
qqr(Prateuc3.pql) # not an improvement

# How do these really work?
summary(Nrateuc2.pql) 
summary(euc2deps2$stockratio[euc2deps2$element=='N'&euc2deps2$year=='04'&
                               euc2deps2$biome=='AF'])
exp(-.642)/(exp(-.642)+1) #the mean, approx
summary(euc2deps2$stockratio[euc2deps2$element=='N'&euc2deps2$year=='16'&
                               euc2deps2$biome=='AF'])
(exp(-.642-.221)/(exp(-.642-.221)+1)) #ok

######## Row position and other heterogeneity things
##################
Celt.lme=lme(log(value)~elt, random=~1|stand,
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='C',],na.action=na.omit)
qqr(Celt.lme) # Nice (with log) 
summary(Celt.lme) # L significantly less than E
xyplot(value~as.numeric(elt)|stand,groups=rep,
       data=dats[dats$depth==5&dats$element=='C'& dats$year=='16' &
                   dats$elt %in% c('E','L','T'),],pch=19)
eltdats=dats[dats$elt %in% c('E','L','T') & dats$LU=='E'&dats$year=='16',]
eltdats=mutate(eltdats,eltlong=ifelse(elt=='E','Inter-',
                                      ifelse(elt=='L','Current','Previous')))
eltdats$eltlong=factor(eltdats$eltlong,levels=c('Current','Previous',
                                                'Inter-'))
eltdats=group_by(eltdats,stand,depth,element,rep,elt)%>%
  mutate(eltmn=mean(value,na.rm=T))

ggplot(data=eltdats[eltdats$element=='P2'&
                      eltdats$stand %in% c('Bp.E1','Eu.E2','It.E1','JP.E2')&
                      eltdats$depth==5,],
       aes(x=eltlong,y=eltmn/1000, color=as.factor(rep)))+
  geom_point(shape=18,size=3,show.legend=F)+
  geom_line(aes(group=as.factor(rep)),show.legend = F)+
  facet_wrap(~stand,ncol=2,scales='free_y')+
  theme(panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank())+
  labs(y='Soil P content (g/kg), 0-10 cm',x='Row position')+
  scale_color_brewer(palette='Dark2')

eltdats=ungroup(eltdats)
eltdats=group_by(eltdats,stand,depth,element,rep)%>%
  mutate(repmn=mean(value,na.rm=T),repvar=var(value,na.rm=T),
         repcv=sqrt(repvar)/repmn)
eltdats=ungroup(eltdats)
eltdats=group_by(eltdats,stand,depth,element,elt)%>%
  mutate(posmn=mean(value,na.rm=T),posvar=var(value,na.rm=T),
         poscv=sqrt(posvar)/posmn)
eltdats=ungroup(eltdats)
eltdis=distinct(eltdats,stand,depth,element,.keep_all = T)
eltsum=group_by(eltdis[eltdis$element %in% c('C','N','P2','K','Ca2')&
                         eltdis$depth==5,],element) %>%
  summarise(#mnposvar=mean(posvar,na.rm=T),mnrepvar=mean(repvar,na.rm=T),
    #nstands=n(), #n=8
    mnposcv=mean(poscv,na.rm=T),mnrepcv=mean(repcv,na.rm=T))
eltsum
#write.csv(eltsum,'eltrepcvs.csv')
#write.csv(eltsum,'eltrepcvs_3-20-19.csv') 
# with updated compositing/averaging procedure, no longer consistently 
#   more variable within position than within rep
t.test(eltdis[eltdis$element %in% c('C','N','P2','K','Ca2')&
                eltdis$depth==5,]$poscv,
       eltdis[eltdis$element %in% c('C','N','P2','K','Ca2')&
                eltdis$depth==5,]$repcv)
# anova probably more appropriate--did I do that already?
# No, because there should not be a consistent effect of rep across sites
# It would have to be n=16 within sites, ok
t.test(eltdis[eltdis$element =='C' & eltdis$depth==5,]$posvar,
       eltdis[eltdis$element =='C' & eltdis$depth==5,]$repvar)
# cv or var not significantly diff between rep and position for any element

# CV within 16 simple samples
eltcvs5=eltdats[eltdats$depth==5 &eltdats$element %in%
                  c('C','N','P2','K','Ca2'),] %>% 
  group_by(element,stand) %>%
  mutate(stdmn=mean(value,na.rm=T),stdvar=var(value,na.rm=T),
         stdcv=sqrt(stdvar)/stdmn,nobs16=n())
eltcvs5=ungroup(eltcvs5)
eltdis5=distinct(eltcvs5,stand,element,.keep_all = T)
table(eltdis5$nobs16[eltdis5$element=='C']) #Eu.E1 one has just 14 obs
eltdis5$stand[eltdis5$nobs16==14] # Eu.E1
eltsum5=group_by(eltdis5,element) %>%
  summarise(mncv=mean(stdcv,na.rm=T),nobs=n())
eltsum5

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
# More plots in longer script

Nelt.lme=lme(value~elt, random=~1|stand,
             #data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
            #             dats$element =='N',],na.action=na.omit) # log makes qqr worse
            data=eltdats[eltdats$depth==5 &eltdats$element=='N',],na.action=na.omit)
summary(Nelt.lme) # L lower again at p=.06 # now 0.09 
# (take out JP.A, also new compositing since this was last run)
Neltbm.lme=lme(value~elt*biome, random=~1|stand,
               data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                           dats$element =='N',],na.action=na.omit)
summary(Neltbm.lme) # L < E for Cerrado at p=.02; T < E at .06


Kelt.lme=lme(log(value)~elt, random=~1|stand,
             data=eltdats[eltdats$depth==5 &eltdats$element =='K',],na.action=na.omit)
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

Cuelt.lme=lme(value~elt, random=~1|stand,
              data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                          dats$element =='Cu',],na.action=na.omit)
qqr(Cuelt.lme) # no good
summary(Caelt.lme)

# Spatial heterogeneity overall:
depcvs=group_by(droplevels(dats4[dats4$site!='TM'&dats4$site!='Cr'&
                                   dats4$LU!='A',]),
                stand,LU,biome,depth,element,year) %>%
  summarise(CV=sd(repval,na.rm=T)/mean(repval,na.rm=T))
tapply(depcvs[depcvs$element=='C',]$CV,depcvs[depcvs$element=='C',]$LU,summary)
tapply(depcvs[depcvs$element=='C',]$CV,depcvs[depcvs$element=='C',]$depth,summary)
tapply(droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K','S','Zn')&
                           depcvs$LU=='E',])$CV,
       droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K','S','Zn')&
                           depcvs$LU=='E',])$element,summary)

xyplot(CV~depth|element,data=depcvs[depcvs$year=='16' &depcvs$element %in% 
                                      c('C','N','P2','Ca2','K','S'),],
       groups=stand,type='l',ylim=c(0,1)) #ugly
depCV.lme=lme(CV~depth,random=~1+depth|element,
              depcvs[depcvs$year=='16' &depcvs$element %in% 
                       c('C','N','P2','Ca2','K','S'),],
              na.action = na.omit)
summary(depCV.lme) # CV increases with depth at p=.06?
qqr(depCV.lme) #no way

# paired sites, by land use
dep16cvs=group_by(droplevels(dats[dats$position %in% c('Ee1','Ee2','L','T') &
                                    dats$LU!='A' & dats$year=='16',]),
                  stand,LU,biome,depth,element) %>%
  summarise(CV=sd(value,na.rm=T)/mean(value,na.rm=T))
tapply(droplevels(dep16cvs[dep16cvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn'),])$CV,
       droplevels(dep16cvs[dep16cvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn'),])$element,summary)

plot(CV~LU,data=depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])
plot(CV~LU,data=depcvs[depcvs$element=='C',])
plot(CV~depth,data=depcvs[depcvs$element=='C' & depcvs$LU!='A',],col=LU,pch=5,cex=2)

bwplot(CV~LU|element,data=depcvs[depcvs$element %in% c('C','N','K') &
                                   depcvs$LU!='A',])
bwplot(CV~as.factor(depth)|LU,data=depcvs[depcvs$element %in% c('C','N','K') &
                                            #depcvs$stand!='It.E1'&depcvs$stand!='It.N'&
                                            depcvs$LU!='A',])
bwplot(CV~as.factor(depth)|LU*element,data=depcvs[depcvs$year=='16' &depcvs$element %in% 
                                                    c('C','N','P2','Ca2','K','S','Zn'),],
       type='l',ylim=c(0,1)) #too much

# CV doesn't change a lot with depth; IQR actually larger at deeper depths
#   maybe due to influence of pit samples? Yes, less of an effect without It.E1
cvaov=aov(CV~LU,data=depcvs[depcvs$element %in% c('C','N','P2','Ca2','K'),])
summary(cvaov)
qqr(cvaov) # nooo
cvaov=aov(log(CV)~LU,data=depcvs[depcvs$element %in% 
                                   c('C','N','P2','Ca2','K') &
                                   depcvs$CV>0,])
summary(cvaov) #p=.081
qqr(cvaov) #tails a bit off
TukeyHSD(cvaov) # N maybe a little > P (p=.086)

Ccvaov=aov(log(CV)~LU,data=depcvs[depcvs$element =='C',])
summary(Ccvaov) #fine with log
qqr(Ccvaov)

cvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element %in% 
                                       c('C','N','P2','Ca2','K') &
                                       depcvs$CV>0,])
qqr(cvdeplm)
summary(cvdeplm) # no effect
Ccvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element=='C',])
qqr(Ccvdeplm)
summary(Ccvdeplm) # +.003, p=.096
Cacvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element=='Ca',])
qqr(Cacvdeplm)
summary(Cacvdeplm) # no effect, nor for P
Kcvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element=='K',])
qqr(Kcvdeplm)
plot(resid(Kcvdeplm)~depcvs[depcvs$element=='K',]$LU)
# some differences, but IQRs similar
# lower CV for pasture, generally
summary(Kcvdeplm)# there is an effect for K, decreasing with depth

Kcvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element=='K'&
                                        depcvs$stand!='Bp.E1'&
                                        depcvs$stand!='Bp.E2',])
qqr(Kcvdeplm)
summary(Kcvdeplm)# without those, same effect size, lower p-value

Ncvdeplm=lm(log(CV)~depth,data=depcvs[depcvs$element=='N' &
                                        depcvs$CV>0,])


# Why is there a zero?
# 2004 values for N are same for all reps in JP.E1.20-40 and JP.P.10-20
qqr(Ncvdeplm) # tail off, non-log bad too
summary(Ncvdeplm) # no effect



Ccvaov=aov(log(CV)~LU,data=depcvs[depcvs$element=='C',])
qqr(Ccvaov) # nice with log
summary(Ccvaov) # no signif effect of LU 
TukeyHSD(Ccvaov)

stockcvs=group_by(dats2deps,stand,LU,biome,element,year) %>%
  summarise(CV20=sd(stock20,na.rm=T)/mean(stock20,na.rm=T),
            CV100=sd(stock100,na.rm=T)/mean(stock100,na.rm=T))
tapply(droplevels(stockcvs[stockcvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn')&
                             stockcvs$LU=='E',])$CV20,
       droplevels(stockcvs[stockcvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn')&
                             stockcvs$LU=='E',])$element,summary)
##################


tstock$stock100_16[tstock$stand=='Eu.E1'&tstock$element=='N']
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

budgets$Nutrient=factor(budgets$Nutrient,levels = 
                          c('N','P','K','Ca','Mg'))
bdgsum=group_by(budgets,Nutrient)%>%
  summarise(Input=mean(ifelse(Wood_m3_2>0,In_kgha_2,In_kgha_1)),
            Harvest=mean(ifelse(Wood_m3_2>0,Wood_m3_2,Wood_m3_1))*
              mean(Concentration)*511,Conc=mean(Concentration),
            Budget=Input-Harvest)
print.data.frame(bdgsum)


summary(stkchgs$minbudg[stkchgs$element=='N']*1000)
# max net N export = 478 kg ha-1 over the 12 years


#tapply(stkchgs$chgrt20,stkchgs$element,mean)
#tapply(stkchgs$chgln20,stkchgs$element,mean)
#tapply(stkchgs$budget,stkchgs$element,mean)
#tapply(stkchgs$efs20,stkchgs$element,mean)
#log(abs(tapply(stkchgs$budget,stkchgs$element,mean)))

t.test(stkchgs2$chg20,stkchgs2$budget,paired=T) 
t.test(stkchgs$chg20[stkchgs$element=='N'],
       stkchgs$budget[stkchgs$element=='N'],paired=T)
t.test(stkchgs$chg20[stkchgs$element=='P'],
       stkchgs$budget[stkchgs$element=='P'],paired=T)
t.test(stkchgs$chg20[stkchgs$element=='K'],
       stkchgs$budget[stkchgs$element=='K'],paired=T)
t.test(stkchgs$chg20[stkchgs$element=='Ca'],
       stkchgs$budget[stkchgs$element=='Ca'],paired=T)
t.test(stkchgs$chg20[stkchgs$element=='Ca'],
       stkchgs$agbbudg[stkchgs$element=='Ca'],paired=T)
# No significant differences
qqnorm(stkchgs$chg20[stkchgs$element=='K'])
qqline(stkchgs$chg20[stkchgs$element=='K'])
qqnorm(stkchgs$budget[stkchgs$element=='K'])
qqline(stkchgs$budget[stkchgs$element=='K']) # not terrible

t.test(stkchgs$chg20,stkchgs$lessrotbudg,paired=T) 
t.test(stkchgs$chg20,stkchgs$woodonlybudg,paired=T) # p=.059
# no longer signif with new, more limited data (no Bp)? 
t.test(stkchgs$bark20budg,stkchgs$woodonlybudg,paired=T) 
# those are different, good


stkchgs3=stkchgs[stkchgs$element!='Mg',]
# Text after Figure 4
summary(I(abs(stkchgs3$maxbudgconc-stkchgs3$minbudgconc)/
            abs(stkchgs3$budget)))
summary(I(abs(stkchgs3$maxbudg-stkchgs3$minbudg)/
            abs(stkchgs3$budget)))
summary(I(abs(stkchgs3$maxagbbudg-stkchgs3$minagbbudg)/
            abs(stkchgs3$agbbudg)))
summary(I((stkchgs3$chg20-stkchgs3$budget)/
            stkchgs3$budget))
summary(I((stkchgs3$chg20-stkchgs3$agbbudg)/
            stkchgs3$agbbudg))

# Where does agb change the sign to match the observations?
stkchgs3$budget[stkchgs3$element=='N']
stkchgs3$agbbudg[stkchgs3$element=='N']
stkchgs3$chg20[stkchgs3$element=='N']
stkchgs3$stand[stkchgs3$element=='N']

# Deprecated Figure 5
#######################
palette('default')

par(mar=c(4,4,0.5,0.5),mfrow=c(1,2))
plot(chg20~budget,data=stkchgs3,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',las=1)
rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=stkchgs$minbudg,x1=stkchgs$maxbudg,y0=stkchgs$chg20,
         col=as.numeric(as.factor(stkchgs$biome))+2)
#segments(x0=stkchgs3$minbudgconc,x1=stkchgs3$maxbudgconc,y0=stkchgs3$chg20,
#         col=as.numeric(as.factor(stkchgs3$biome))+2)
segments(x0=stkchgs3$budget,y0=stkchgs3$chg20-stkchgs3$sdchg20,
         y1=stkchgs3$chg20+stkchgs3$sdchg20,
         col=as.numeric(as.factor(stkchgs3$biome))+2)

text(stkchgs3$budget,stkchgs3$chg20,labels=stkchgs3$element)#,
#     col=as.numeric(as.factor(stkchgs3$biome))+2)
legend('bottomright',pch=15,col=c(3,4),
       legend=c('Atlantic Forest','Cerrado'),bty='n')
legend('bottomleft',legend='a',cex=1.2,bty='n')
plot(chg20~budget,data=stkchgs3,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='', xlim=c(-.2,.5),ylim=c(-.2,.5),las=1)
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=stkchgs3$minbudg,x1=stkchgs3$maxbudg,y0=stkchgs3$chg20,
         col=as.numeric(as.factor(stkchgs3$biome))+2)
#segments(x0=stkchgs3$minbudgconc,x1=stkchgs3$maxbudgconc,y0=stkchgs3$chg20,
#         col=as.numeric(as.factor(stkchgs3$biome))+2)
segments(x0=stkchgs3$budget,y0=stkchgs3$chg20-stkchgs3$sdchg20,
         y1=stkchgs3$chg20+stkchgs3$sdchg20,
         col=as.numeric(as.factor(stkchgs3$biome))+2)
text(stkchgs3$budget,stkchgs3$chg20,labels=stkchgs3$element)
legend('bottomleft',legend='b',cex=1.2,bty='n')
par(mfrow=c(1,1),mar=c(4,4,1,1))

# presentation version?
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")
plot(x = 1:length(cbbPalette), y = rep(1, length(cbbPalette)), pch = 19, cex = 5, 
     col = cbbPalette, yaxt = "n", bty = "n", xaxt = "n", xlab = "", ylab = "")

stkchgs3=stkchgs[stkchgs$element!='Mg',]
palette(cbbPalette)
plot(chg20~budget,data=stkchgs3,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',
#     xlim=c(-.2,.5),ylim=c(-.2,.5),
     las=1)
#rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=stkchgs3$minbudg,x1=stkchgs3$maxbudg,y0=stkchgs3$chg20,
         col=as.numeric(as.factor(stkchgs3$element))*2)
#segments(x0=stkchgs3$minbudgconc,x1=stkchgs3$maxbudgconc,y0=stkchgs3$chg20,
#         col=as.numeric(as.factor(stkchgs3$biome))+2)
segments(x0=stkchgs3$budget,y0=stkchgs3$chg20-stkchgs3$sdchg20,
         y1=stkchgs3$chg20+stkchgs3$sdchg20,
         col=as.numeric(as.factor(stkchgs3$element))*2)
text(stkchgs3$budget,stkchgs3$chg20,labels=stkchgs3$element,
     col=as.numeric(as.factor(stkchgs3$element))*2)
text(stkchgs3$agbbudg,stkchgs3$chg20,labels=stkchgs3$element)

palette('default')
budgplot=function(elmt,limfac=1.01){
  par(mar=c(4,4,1,1))
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg20-datsub$sdchg20)*limfac
  ymax=max(datsub$chg20+datsub$sdchg20)*limfac
  
plot(chg20~budget,data=datsub,#type='n', 
     xlab='Expected budget, Mg/ha',
     ylab='Observed stock change to 20 cm, Mg/ha',
          xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=2,
     las=1,cex.axis=1.3,cex.lab=1.3)
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=datsub$minbudg,x1=datsub$maxbudg,y0=datsub$chg20,
         lwd=1.5)#,
         #col=as.numeric(as.factor(datsub$element))*2)
segments(x0=datsub$budget,y0=datsub$chg20-datsub$sdchg20,
         y1=datsub$chg20+datsub$sdchg20,lwd=1.5)#,
         #col=as.numeric(as.factor(datsub$element))*2)
#text(datsub$budget,datsub$chg20,labels=datsub$element,
#     col=as.numeric(as.factor(datsub$element))*2)
#text(datsub$agbbudg,datsub$chg20,labels=datsub$element)
points(datsub$agbbudg,datsub$chg20,cex=2,col=3,pch=16)
segments(x0=datsub$minagbbudg,x1=datsub$maxagbbudg,y0=datsub$chg20,
         lty=2,col=3,lwd=1.5)#,
}

budgplotsimp=function(elmt,limfac=1.01){
  par(mar=c(4,4,1,1))
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg20-datsub$sdchg20)*limfac
  ymax=max(datsub$chg20+datsub$sdchg20)*limfac
  
  plot(chg20~budget,data=datsub,#type='n', 
       xlab='Expected budget, Mg/ha',
       ylab='Observed stock change to 20 cm, Mg/ha',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=2,
       las=1,cex.axis=1.3,cex.lab=1.3)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  abline(0,1)
  segments(x0=datsub$minbudg,x1=datsub$maxbudg,y0=datsub$chg20,
           lwd=1.5)#,
  #col=as.numeric(as.factor(datsub$element))*2)
  segments(x0=datsub$budget,y0=datsub$chg20-datsub$sdchg20,
           y1=datsub$chg20+datsub$sdchg20,lwd=1.5)#,
}
budgplot100=function(elmt,limfac=1.01){
  par(mar=c(4,4,1,1))
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg100-datsub$sdchg100)*limfac
  ymax=max(datsub$chg100+datsub$sdchg100)*limfac
  
  plot(chg100~budget,data=datsub,#type='n', 
       xlab='Expected budget, Mg/ha',
       ylab='Observed stock change to 100 cm, Mg/ha',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=2,
       las=1,cex.axis=1.3,cex.lab=1.3)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  abline(0,1)
  segments(x0=datsub$minbudg,x1=datsub$maxbudg,y0=datsub$chg100,
           lwd=1.5)#,
  #col=as.numeric(as.factor(datsub$element))*2)
  segments(x0=datsub$budget,y0=datsub$chg100-datsub$sdchg100,
           y1=datsub$chg100+datsub$sdchg100,lwd=1.5)#,
  #col=as.numeric(as.factor(datsub$element))*2)
  #text(datsub$budget,datsub$chg100,labels=datsub$element,
  #     col=as.numeric(as.factor(datsub$element))*2)
  #text(datsub$agbbudg,datsub$chg100,labels=datsub$element)
  points(datsub$agbbudg,datsub$chg100,cex=2,col=3,pch=16)
  segments(x0=datsub$minagbbudg,x1=datsub$maxagbbudg,y0=datsub$chg100,
           lty=2,col=3,lwd=1.5)#,
}

png('Nbudgplot_green.png',width=6,height=6,units='in',res=150)
budgplot('N',1.01)
legend('topleft',legend='Nitrogen',cex=1.5,bty='n')
legend('bottomleft',bty='n',cex=1.2,col=c(1,3),pch=16,#pch=c(16,1),
       legend=c('Fertilizer - Harvest',
'Fertilizer - Harvest +\nInput from biomass change'))
dev.off()

budgplotsimp('K',type='n')
budgplot('K')

png('Nbudgplot_simp.png',width=6,height=6,units='in',res=150)
budgplotsimp('N',1.01)
legend('topleft',legend='Nitrogen',cex=1.5,bty='n')
dev.off()

png('Kbudgplot_simp.png',width=6,height=6,units='in',res=150)
budgplotsimp('K',1.01)
legend('topleft',legend='Potassium',cex=1.5,bty='n')
dev.off()

budgplot100('N')

budgplotnull=function(elmt,limfac=1.01){
  par(mar=c(4,4,1,1))
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg20-datsub$sdchg20)*limfac
  ymax=max(datsub$chg20+datsub$sdchg20)*limfac
  
  plot(chg20~budget,data=datsub,type='n', 
       xlab='Expected budget, Mg/ha',
       ylab='Observed stock change to 20 cm, Mg/ha',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=2,
       las=1,cex.axis=1.3,cex.lab=1.3)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  abline(0,1)
}
budgplotnull('K')

png('Kbudgplot_null.png',width=6,height=6,units='in',res=150)
budgplotnull('K')
legend('topleft',legend='Potassium',cex=1.5,bty='n')
dev.off()

# New Figure 4?

budgplotfig=function(elmt,limfac=1.01){
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg20-datsub$sdchg20)*limfac
  ymax=max(datsub$chg20+datsub$sdchg20)*limfac
  
  plot(chg20~budget,data=datsub,#type='n', 
       xlab='',ylab='',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=1.5,
       las=1)#,cex.axis=1.2)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  abline(0,1)
  segments(x0=datsub$minbudg,x1=datsub$maxbudg,y0=datsub$chg20,
           lwd=1.5)#,
  #col=as.numeric(as.factor(datsub$element))*2)
  segments(x0=datsub$budget,y0=datsub$chg20-datsub$sdchg20,
           y1=datsub$chg20+datsub$sdchg20,lwd=1.5)#,
}


#dev.new(width=6,height=6,units='in',res=72)
png('possible_fig4.png',width=6.5,height=6.5,units='in',res=150)

par(mfrow=c(2,2),mar=c(1.5,4,4,1))
budgplotfig('N')
legend('topleft',legend='N',cex=1.2,bty='n') 
title(ylab='Observed Δ stock to 20 cm, Mg/ha')#,cex.lab=1.2)
budgplotfig('P')
legend('topleft',legend='P',cex=1.2,bty='n')
par(mar=c(4,4,1.5,1))
budgplotfig('K')
legend('topleft',legend='K',cex=1.2,bty='n')
title(xlab='Expected budget (fertilizer - harvest), Mg/ha')#,cex.lab=1.2)
budgplotfig('Ca')
legend('topleft',legend='Ca',cex=1.2,bty='n')
dev.off()

budgplotfig2=function(elmt,limfac=1.01){
  datsub=stkchgs[stkchgs$element==elmt,]
  xmin=min(datsub$minbudg)*limfac
  xmax=max(datsub$maxbudg)*limfac
  ymin=min(datsub$chg20-datsub$sdchg20)*limfac
  ymax=max(datsub$chg20+datsub$sdchg20)*limfac
  
  plot(chg20~budget,data=datsub,#type='n', 
       xlab='',ylab='',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=1.2,
       las=1)#,cex.axis=1.2)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  abline(0,1)
  segments(x0=datsub$minbudg,x1=datsub$maxbudg,y0=datsub$chg20,
           lwd=1.5)#,
  #col=as.numeric(as.factor(datsub$element))*2)
  segments(x0=datsub$budget,y0=datsub$chg20-datsub$sdchg20,
           y1=datsub$chg20+datsub$sdchg20,lwd=1.5)#,
  points(datsub$agbbudg,datsub$chg20,cex=1.2,col=3,pch=16)
  segments(x0=datsub$minagbbudg,x1=datsub$maxagbbudg,y0=datsub$chg20,
           lty=2,col=3,lwd=1.5)
}

# Supplemental figure to accompany?
#dev.new(width=6,height=3,units='in')
png('possible_supplmental_fig4.png',width=6,height=3,units='in',res=150)
par(mfrow=c(1,2),mar=c(4,4,1,1))
budgplotfig2('N',1.03)
legend('bottomleft',legend='N',cex=1.2,adj=c(1,1),bty='n')
title(ylab='Observed Δ stock to 20 cm, Mg/ha',cex.lab=.9)
title(xlab='Expected budget, Mg/ha',cex.lab=.9)
budgplotfig2('K',1.03)
legend('bottomleft',legend='K',adj=c(1,1),cex=1.2,bty='n')
legend('topleft',bty='n',col=c(1,3),adj=c(0,0),pch=16,cex=.8,
       legend=c('Fertilizer - Harvest',
                'Fertilizer - Harvest +\nInput from biomass change'))
dev.off()

par(mfrow=c(2,2),mar=c(1,4,4,1))
budgplotfig2('N')
legend('topleft',legend='N',cex=1.2,bty='n')
title(ylab='Observed stock change to 20 cm, Mg/ha',cex.lab=1.2)
budgplotfig2('P')
legend('topleft',legend='P',cex=1.2,bty='n')
legend('bottomleft',bty='n',cex=1.2,col=c(1,3),pch=16,#pch=c(16,1),
       legend=c('Fertilizer - Harvest',
                'Fertilizer - Harvest +\nInput from biomass change'))

par(mar=c(4,4,1,1))
budgplotfig2('K')
legend('topleft',legend='K',cex=1.2,bty='n')
title(xlab='Expected budget, Mg/ha',cex.lab=1.2)
budgplotfig2('Ca')
legend('topleft',legend='Ca',cex=1.2,bty='n')
# too much


ggplot(aes(x=budget,y=chg20,color=stand),data=stkchgs3)+
  facet_wrap(.~element)+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)
#######################

# most recent rotation
summary(c(budgets$In_kgha_2[budgets$In_kgha_2!=0],
          budgets$In_kgha_1[budgets$In_kgha_2==0])) # by nutrient tho

stkchgs3=stkchgs[stkchgs$element!='Mg',]
stkchgs3$element=factor(stkchgs3$element,levels=c('N','P','K','Ca'))
chgtypes=group_by(stkchgs3,stand,element, biome, chg20,agbbudg,
                  sdchg20,minbudg,minbudgconc,maxbudg,
                  maxbudgconc,budget,minagbbudg,maxagbbudg) %>%
  summarise(standing=agbchg,
            harvest=(Wood_m3_1+Wood_m3_2)*Concentration*-511/1000,
            fertilizer=(In_kgha_1+In_kgha_2)/1000)

# Table 1:
tapply(chgtypes$fertilizer*1000,chgtypes$element,summary)
tapply(chgtypes$harvest*1000,chgtypes$element,summary)
tapply(chgtypes$budget*1000,chgtypes$element,summary)
tapply(stkchgs3$conc*1000,stkchgs3$element,summary)
tapply(stkchgs3$agbchg*-1000,stkchgs3$element,summary)
tapply(chgtypes$agbbudg*1000,chgtypes$element,summary)

summary(c(stkchgs3$Wood_m3_1[stkchgs3$biome=='AF'&stkchgs3$element=='N'],
       stkchgs3$Wood_m3_2[stkchgs3$biome=='AF'&stkchgs3$element=='N'&
                            stkchgs3$Wood_m3_2>0]))
summary(c(stkchgs3$Wood_m3_1[stkchgs3$biome=='Cer'&stkchgs3$element=='N'],
          stkchgs3$Wood_m3_2[stkchgs3$biome=='Cer'&stkchgs3$element=='N'&
                               stkchgs3$Wood_m3_2>0]))
tapply(chgtypes$fertilizer[chgtypes$biome=='AF']*1000,
       chgtypes$element[chgtypes$biome=='AF'],summary)
tapply(chgtypes$fertilizer[chgtypes$biome=='Cer']*1000,
       chgtypes$element[chgtypes$biome=='Cer'],summary)


sensit=group_by(stkchgs3,element,stand)%>%
  summarise(mindivagb=minagbbudg/agbbudg, maxdivagb=maxagbbudg/agbbudg,
         mindiv=minbudg/budget, maxdiv=maxbudg/budget,
         mindivconc=minbudgconc/budget, maxdivconc=maxbudgconc/budget)
sensitsum=group_by(sensit,element)%>%
  summarise_if(is.numeric,list(min,median,mean,max),na.rm=T) %>%
  mutate_if(is.numeric,round,digits=2)
data.frame(t(sensitsum[,order(names(sensitsum))])) #TMI

sensit2=group_by(stkchgs3,element,stand) %>%
  summarise(budgsens=abs(maxbudg-minbudg)/abs(budget),
         concsens=abs(maxbudgconc-minbudgconc)/abs(budget),
         agbsens=abs(maxagbbudg-minagbbudg)/abs(agbbudg))
sensitsum=group_by(sensit2,element)%>%
  summarise_if(is.numeric,list(min,median,mean,max),na.rm=T) %>%
  mutate_if(is.numeric,round,digits=2)
data.frame(t(sensitsum[,order(names(sensitsum))]))

obsdif=group_by(stkchgs3,element,stand) %>%
 # summarise(budgdif=abs(chg20-budget)/budget,
  #          agbdif=abs(chg20-agbbudg)/agbbudg) # add the abs() to show budget sign
  summarise(budgdif=chg20/budget,
            agbdif=chg20/agbbudg) 
  #summarise(budgdif=budget/chg20,
  #          agbdif=agbbudg/chg20) 
  #summarise(budgdif=log(abs(chg20/budget)),
  #          agbdif=log(abs(chg20/agbbudg))) 
obsdifsum=group_by(obsdif,element)%>%
  summarise_if(is.numeric,funs(min,median,mean,max),na.rm=T) %>%
  mutate_if(is.numeric,round,digits=2)
data.frame(t(obsdifsum[,order(names(obsdifsum))]))
data.frame(obsdifsum[,order(names(obsdifsum))])

summary(I((stkchgs3$chg20-stkchgs3$budget)/
            stkchgs3$budget))
summary(I((stkchgs3$chg20-stkchgs3$agbbudg)/
            stkchgs3$agbbudg))

obsdifm=melt(obsdif,id.vars=c('element','stand'),value.name = 'disc')
obsdifm2=merge(obsdifm,obsdifsum,by='element')
obsdifm2a=obsdifm2[obsdifm2$variable=='agbdif',
                   -which(names(obsdifm2) %in% 
                           c('budgdif_min','budgdif_median',
                             'budgdif_mean','budgdif_max'))]
obsdifm2b=obsdifm2[obsdifm2$variable=='budgdif',
                   -which(names(obsdifm2) %in% 
                            c('agbdif_min','agbdif_median',
                              'agbdif_mean','agbdif_max'))]
names(obsdifm2a)
names(obsdifm2b)
names(obsdifm2a)=c('element','stand','budgtype','disc',
                   'min','median','mean','max')
names(obsdifm2b)=c('element','stand','budgtype','disc',
                   'min','median','mean','max')
obsdifm3=data.frame(rbind(obsdifm2b,obsdifm2a))

ggplot(obsdifm,aes(x=element,y=disc,color=variable))+
  geom_boxplot()

# Figure 4
png(filename = 'fig4_newcen_bw.png',width=5,height=5, units='in',res=150)
ggplot(obsdifm3,aes(x=element,y=disc,color=budgtype))+
  scale_color_manual(values=c('grey40','black'),#c('darkblue','darkgreen'),
                      name='Budget',
                      labels=c('Fertilizer - Harvest',
                               'Fertilizer - Harvest +\nInput from biomass change'))+
  guides(color=guide_legend())+
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_point(size=3, shape=1,position=position_dodge(.3),show.legend = F)+
  geom_pointrange(aes(x=element,y=median,colour=budgtype,ymax=max,
                      ymin=min),fatten=6,size=.8,
                  data=distinct(obsdifm3,element,budgtype,.keep_all=T),
                  shape=18,show.legend = T,position=position_dodge2(.3))+
  #labs(y='(Observed change in soil stock - Budget) / Budget', x=NULL) +
  #labs(y='(Measured Δ soil stock - Budget) / Budget', x=NULL) +
  labs(y='Measured Δ soil stock / Budget', x=NULL) +
  theme(legend.position=c(0.52,0.15),#c(0.22,0.9),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_line(colour='grey80'),
        legend.key=element_blank(),
        legend.title=element_text(size=10)) 
dev.off()  
  


chgtypesm=melt(chgtypes,measure.vars = c('standing','harvest','fertilizer'))

tapply(chgtypes$budget,chgtypes$element,summary)
stdlabls <- c(BO.E = "Atlantic Forest example", 
           It.E2 = "Cerrado example") # or Eu.E1 for AF?

ggplot(chgtypesm[chgtypesm$stand %in% c('BO.E','It.E2'),],
       #[chgtypesm$stand %in% c('Eu.E2','It.E2','JP.E2'),],
       aes(x=element,y=value,fill=variable))+
  #geom_bar(stat = "identity")+
  facet_wrap(~stand, labeller=labeller(stand=stdlabls))+
  #facet_grid(rows=vars(stand))+
  #geom_point(aes(y=chg20))+
  #geom_point(aes(y=budget),shape=1)+
  #geom_errorbar(aes(ymin=minbudg,ymax=maxbudg),width=.2)
  #scale_fill_brewer(palette='Dark2',name=NULL,#name='Change in pool',
  #                  labels=c('Transfer to/from biomass',
  #                           'Removal in harvested wood',
  #                           'Addition in fertilizer'))+
  coord_cartesian(ylim=c(-1,1))+
  geom_hline(yintercept=0,color='gray50')+
  geom_point(aes(y=chg20),position=position_nudge(x=.2),
             show.legend = F)+
  geom_point(aes(y=budget,color='Budget'),
             position=position_nudge(x=-.2),
             show.legend = F)+
  geom_errorbar(aes(ymin=minbudg,ymax=maxbudg,color='Budget'),
                width=.2,position=position_nudge(x=-.2),size=1)+
  geom_errorbar(aes(ymin=chg20-sdchg20/2,ymax=chg20+sdchg20/2,
                    color='Observed'),size=1,
                width=.2,position=position_nudge(x=.2))+
  geom_errorbar(aes(ymin=minagbbudg,ymax=maxagbbudg,
                    color='Budget_standing'),width=.2,size=1)+
  geom_point(aes(y=agbbudg,colour='Budget_standing'),
             size=2,show.legend = F)+
  labs(x='',y='Change in soil nutrient stocks (Mg / ha), 0-20 cm')+
  scale_colour_manual(name=NULL,
                      values=c(Budget="deepskyblue", 
                               Observed="darkblue",
                               Budget_standing='lightgreen'),
                      labels=c('Budget (fertilizer - harvest)',
                               'Budget + change in standing biomass',
                               'Observed change in soil stock'))+
  theme(legend.position=c(0.75,0.2),
    #legend.position=c(0.55,0.8),
    legend.spacing.y = unit(.8,'lines'), 
    panel.background = element_rect(fill='white'),
    #panel.grid.major = element_blank(),
    legend.key=element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=11))#+
  #geom_errorbar(aes(ymin=minbudgconc,ymax=maxbudgconc),
  #              color='orchid',
  #              width=.2,position=position_nudge(x=-.1))
  

chgtypes1per=group_by(chgtypes,element) %>%
  mutate(nobs=n(),sebudg=sd(budget)/sqrt(nobs),
         nagb=sum(!is.na(agbbudg)),
         seagbbudg=sd(agbbudg,na.rm=T)/sqrt(nagb),
         seobs=sd(chg20)/sqrt(nobs)) %>%
  summarise_if((is.numeric), 
               function(x){mean(x[is.finite(x)==1],na.rm=T)})
chgtypes1perm=melt(chgtypes1per,
                   measure.vars = c('standing','harvest','fertilizer'))

chgtypes1permed=group_by(chgtypes,element) %>%
  summarise_if(is.numeric, median, na.rm=T)
chgtypes1permedm=melt(chgtypes1permed,
                   measure.vars = c('standing','harvest','fertilizer'))


ggplot(chgtypes1perm,aes(x=element,y=value,fill=variable))+
  geom_bar(stat = "identity")+
  #facet_wrap(~biome)+
  geom_point(aes(y=chg20,colour='Observed'),size=2,
             show.legend = F,position=position_nudge(x=.3))+
  scale_fill_brewer(palette='Dark2',name='Pool')+
  geom_point(aes(y=budget,colour='Budget'),size=2,
             position=position_nudge(x=-.3),show.legend = F)+
  geom_point(aes(y=agbbudg,colour='Budget_standing'),
             size=2,show.legend = F)+
  geom_errorbar(aes(ymin=minbudg,ymax=maxbudg,
                color='Budget'),width=.2,size=1,
                position=position_nudge(x=-.3))+
  geom_errorbar(aes(ymin=minagbbudg,ymax=maxagbbudg,
                    color='Budget_standing'),width=.2,size=1)+
  geom_errorbar(aes(ymin=chg20-sdchg20,ymax=chg20+sdchg20,
                    color='Observed'),width=.2,size=1,
                position=position_nudge(x=.3))+
  #geom_errorbar(aes(ymin=chg20-seobs,ymax=chg20+seobs,
  #                  color='Observed'),width=.2,size=1,
  #              position=position_nudge(x=.3))+
  #geom_errorbar(aes(ymin=budget-sebudg,ymax=budget+sebudg,
  #                  color='Budget'),width=.2,size=1,
  #              position=position_nudge(x=-.3))+
  #geom_errorbar(aes(ymin=agbbudg-seagbbudg,ymax=agbbudg+seagbbudg,
  #                  color='Budget_standing'),width=.2,size=1)+
  labs(x='',y='Change in soil nutrient stocks (Mg/ha), 0-20 cm')+
  scale_colour_manual(name="Change in stock",
                      values=c(Budget="deepskyblue", 
                               Observed="darkblue",
                               Budget_standing='lightgreen'),
                      labels=c('Budget (fertilizer - harvest)',
                               'Budget + change in standing biomass',
                               'Observed change in soil stock'))+
  theme(#legend.position=c(0.2,0.8),
        legend.position=c(0.4,0.7),
        #legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())#+
  geom_errorbar(aes(ymin=minbudgconc,ymax=maxbudgconc),
                color='orchid',
                width=.2,position=position_nudge(x=-.1))



simple20_2=mutate(simple20_2,pairtype=ifelse(site2 %in% c('JP','BO'),'P','N'))
simple20_2$pairtype=factor(simple20_2$pairtype,levels=c('N','P'))

# is it appropriate to group pasture with native veg? 
# Decided that it's not (use LU, not LU2)
# Allow different year effects in different sites
Casimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='Ca2',], na.action=na.omit)
qqr(Casimp.lme4) # actually better than lme3
summary(Casimp.lme4) # increase in euc, not signif in pasture

Csimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='C',], 
               na.action=na.omit)
qqr(Csimp.lme4) # tails off,esp lower
summary(Csimp.lme4) # no effect of year on euc (p=.09) 
# but sig neg intrxns
#   between year and both native and pasture veg
# changing factor levels: year doesn't have sig effect on pasture C or Ca
# pasture-euc intrxn stays signif, ok
median()

# Test significance of year within each veg type by changing factor levels?
#   Seems wrong...
simple20_2$LU=factor(simple20_2$LU,levels=c('P','E','N'))
simple20_2$LU=factor(simple20_2$LU,levels=c('N','E','P'))
simple20_2$LU=factor(simple20_2$LU,levels=c('E','N','P'))
# Significant things to 20 cm: 
#   increase in euc in Ca 
#   year-native intrxn in Ca and C
#   year-pasture intrxn in C

#Nsimp.lme3=lme(stock20~year*LU,random=~1|site2,
#               data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
#qqr(Nsimp.lme3) # nice-ish...tails off; better without log
#summary(Nsimp.lme3) # with all sites included, starts higher in N, no change
Nsimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='N',], na.action=na.omit)
qqr(Nsimp.lme4) 
summary(Nsimp.lme4) # no change; started with more N in native
#anova(Nsimp.lme3,Nsimp.lme4) # 4 better fit

Ksimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='K',], na.action=na.omit)
qqr(Ksimp.lme4) 
summary(Ksimp.lme4) # marginal year-native intrxn (p=.06)
#anova(Ksimp.lme3,Ksimp.lme4) # not different

Psimp.lme4=lme(stock20~year*LU,random=~1+year|site2,
               data=simple20_2[simple20_2$element=='P2',], na.action=na.omit)
qqr(Psimp.lme4) # worse
summary(Psimp.lme4) # no signif changes as you might expect

Zrsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Zr',], na.action=na.omit)
qqr(Zrsimp.lme4) 
summary(Zrsimp.lme4) # significant year-native intrxn

Nbsimp.lme4=lme(log(stock20)~year*LU,random=~1+year|site2,
                data=simple20_2[simple20_2$element=='Nb',], na.action=na.omit,
                control = lmeControl(opt = 'optim'))
qqr(Nbsimp.lme4) 
summary(Nbsimp.lme4) # no year or interaction effects
# but initially lower value in native and higher in pasture
# is that also true within pairs?
plot(conc20~as.numeric(as.factor(site2)),col=as.factor(LU),
     data=simple20_2[simple20_2$element=='Nb',],pch=as.numeric(as.factor(year)))
plot(conc20~as.factor(site2),col=as.factor(LU),
     data=simple20_2[simple20_2$element=='Nb',],pch=as.numeric(as.factor(year)))
# Changes similarly in both in Vg, more in euc in JP.E2 than in JP.N (both yrs),
#   more in JP.E1 than JP.P both years
# Makes sense that JP pairing less accurate: paired stands more spatially separated,
#   lower element concentrations and sandier soil increase measurement error

CNsimp.lme=lme(conc20~year*LU,random=~1+year|site2,
  data=simple20_2[simple20_2$element=='CN',], na.action=na.omit)
summary(CNsimp.lme) # increases in euc at p=.046
# signif intrxn with native and pasture 
#   but the decreases in those veg types aren't significantly nonzero
qqr(CNsimp.lme) # pretty ok

nat2deps=droplevels(dats2deps[dats2deps$LU=='N',])
nat2deps$biome=factor(nat2deps$biome,levels=c('Cer','AF'))
natC20bm.lme=lme(log(stock20)~year*biome,random=~1|stand,
                 data=nat2deps[nat2deps$element=='C',],na.action = na.omit)
qqr(natC20bm.lme) # tails a bit off
summary(natC20bm.lme) # marginal decrease in AF, def increase in cerrado
exp(-.1835)-1 # mean change in AF = 1
exp(-.1835)*exp(.3062)-1 # in Cerrado

natC20.lme=lme(log(stock20)~year,random=~1|stand,
               data=nat2deps[nat2deps$element=='C',],na.action = na.omit)
qqr(natC20.lme) # ok
summary(natC20.lme) # no significant change

natN20bm.lme=lme(log(stock20)~year*biome,random=~1|stand,
                 data=nat2deps[nat2deps$element=='N',],na.action = na.omit)
qqr(natN20bm.lme) # nice
summary(natN20bm.lme) # marginal decrease in AF, def increase in cerrado
# also marginal decrease in AF (p=.076), strong increase in Cerrado

exp(-.1615)-1 # mean change in AF = 15% decrease
exp(-.1615)*exp(.6178)-1 # in Cerrado

natCN20bm.lme=lme(log(conc20)~year*biome,random=~1|stand,
                  data=nat2deps[nat2deps$element=='CN',],na.action = na.omit)
qqr(natCN20bm.lme) # ok?
summary(natCN20bm.lme) # marginal decrease in Cerrado (p=.07)

natCa20.lme=lme(log(stock20)~year*biome,random=~1|stand,
               data=nat2deps[nat2deps$element=='Ca2',],na.action = na.omit)
qqr(natCa20.lme) # ok
summary(natCa20.lme) # no significant change


# Removed pairs with known issues below 20 cm
simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]
# Fix unbalanced design that causes spurious decreases in K
#   because there are fewer measurements of rock-derived nutrients
#   in Vg.N than in Eu.N in 04 (but not in 16) due to one missing sample
# No, that wasn't the problem. Don't do this.
#simp100$stock100[simp100$element=='K'&simp100$year=='04'&
#                   simp100$stand=='Vg.N'&simp100$rep==1]=
#  mean(simp100$stock100[simp100$element=='K'&simp100$year=='04'&
#                          simp100$stand=='Vg.N'&simp100$rep!=1],na.rm=T)
#simp100$stock100[simp100$element=='P2'&simp100$year=='04'&
#                   simp100$stand=='Vg.N'&simp100$rep==1]=
#  mean(simp100$stock100[simp100$element=='P2'&simp100$year=='04'&
#                          simp100$stand=='Vg.N'&simp100$rep!=1],na.rm=T)
#simp100$stock100[simp100$element=='Ca2'&simp100$year=='04'&
#                   simp100$stand=='Vg.N'&simp100$rep==1]=
#  mean(simp100$stock100[simp100$element=='Ca2'&simp100$year=='04'&
#                          simp100$stand=='Vg.N'&simp100$rep!=1],na.rm=T)


Ksimp100.lme=lme(log(stock100)~year*LU,random=~1|site,
                 data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme) # ok
summary(Ksimp100.lme) # with random intercept only, euc increases, others not different
Ksimp100.lme2=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme2) # mostly ok
summary(Ksimp100.lme2) # with random slopes, no change for euc, 
# signif negative native-euc term
# but when N set to default, decrease is not significant
# error on estimates means that when they're added, mean chg comes out wrong
# K in N increases less than in either E or P
# ok
exp(.071)
exp(.071)*exp(-.131)
# so is it misleading to estimate effect size from lme?
###### native decr #### not really a problem? skip exploration below
##################
# Wait, decreases? No, it doesn't
# Increases significantly less than the increase in euc, 
#   but the increase in euc is not significant.
# Ugh, am I misusing these statistics? I don't think so...
# With native veg as the default, NO significant year effect,
#   but significant interactions with other veg types
# pasture: maybe marginal increase, intrxn with native veg is signif
Ksimp100.lme0=lme(log(stock100)~year+LU,random=~1+year|site,
                  data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme0) # mostly ok
summary(Ksimp100.lme0) # no effect of year
anova(Ksimp100.lme0,Ksimp100.lme2)
anova(Ksimp100.lme,Ksimp100.lme2) # 2 is better
Ksimp100.lm=lm(log(stock100)~year*LU*site,
                data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lm)
summary(Ksimp100.lm)
# signif effects of year, LUs
# why is there no LUN:siteVg term?? Just NAs.
simp100$stock100[simp100$element=='K'&simp100$stand=='Vg.N']
# there are two NaNs in there (reps 1 and 4 of 04), does that ruin things?
Ksimp100.lm=lm(log(stock100)~year*LU*site,
               data=simp100[simp100$element=='K'&!is.na(simp100$stock100),])
qqr(Ksimp100.lm)
summary(Ksimp100.lm)# still not there
summary(simp100$stock100[simp100$element=='K'&simp100$year=='04'])
summary(simp100$stock100[simp100$element=='K'&simp100$year=='16'])
# The median decreases, what do you know
#  bunch of NAs in 2004...
summary(simp100$stock100[simp100$element=='K'&
                           simp100$LU=='E'&
                           simp100$year=='04'])
summary(simp100$stock100[simp100$element=='K'&
                           simp100$LU=='E'&
                           simp100$year=='16']) #median decr
summary(simp100$stock100[simp100$element=='K'&
                           simp100$LU=='N'&
                           simp100$year=='04'])
summary(simp100$stock100[simp100$element=='K'&
                           simp100$LU=='N'&
                           simp100$year=='16']) #also

# mean and median thrown off by diff. no. of reps, lme estimates not
summary(mrstkchgs$chgln100[mrstkchgs$element=='K']) 
# median change is positive, change in medians is negative
simp100$stand[simp100$element=='K'&is.na(simp100$stock100)]
# several. What's up with that? Two BO.E, one JP.E1, one Vg.E, two Vg.N
# also missing values there for P, fewer for N
# But the decrease is really down to the big change in Eu.N 0-10 cm
# Does lme get thrown off by different number of reps in Eu vs Vg in 04?
#   because of the missing values in 04? Seems like it shouldn't
summary(simp100$stock100[simp100$element=='K'&simp100$year=='04'&
                           simp100$LU=='N'])
summary(simp100$stock100[simp100$element=='K'&simp100$year=='16'&
                           simp100$LU=='N']) # median and mean decrease
# now median decreases barely, mean increases barely with balanced design

Ksimp100.lme$coefficients
# just down to differences between sites, I think
# the mean effect is less than the site effect
#   which is why the change without site effect is different
# also the euc year estimate is .07 +/- .08
# intrxn estimate -.13 +/- .06
exp(.14)*exp(-.13) # increase
# leave as is?

simp100$stock100[simp100$element=='K'&simp100$year=='04'&
                   simp100$stand=='Eu.N']
simp100$stock100[simp100$element=='K'&simp100$year=='16'&
                   simp100$stand=='Eu.N']
simp100$stock100[simp100$element=='K'&simp100$year=='04'&
                   simp100$stand=='Vg.N']
simp100$stock100[simp100$element=='K'&simp100$year=='16'&
                   simp100$stand=='Vg.N']

# 04 = 4 samples of high-K Eu and 3 of low-Ca Vg
# 16 = 4 samples of each
# so K decreases
# should also be a problem for Ca and P
# just for summary stats, not for lmes
Ksimp100AF.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                  data=simp100[simp100$element=='K'&
                                 simp100$biome=='AF',],
                  na.action=na.omit)
qqr(Ksimp100AF.lme) # mostly ok
summary(Ksimp100AF.lme) #  marginal increase in euc, not in nat
# change in nat still a little negative
# increase in P? just the one pasture, yeah

##################


Nsimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='N',], na.action=na.omit)
qqr(Nsimp100.lme) # ok
summary(Nsimp100.lme)
#With random slope, nothing significant 

Csimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='C',], na.action=na.omit)
qqr(Csimp100.lme)
summary(Csimp100.lme) 
# with random ~1|site, increase in euc isn't signif, 
#   but decreases in native and pasture are! 
# with random ~1|year+site (more appropriate),
#   only signif thing is decrease in pasture (in native, intrxn p=.07)
#   and when native set to default level, no year effect
summary(simp100$stock100[simp100$element=='C'&
                           simp100$LU=='P'&
                           simp100$year=='04'])
summary(simp100$stock100[simp100$element=='C'&
                           simp100$LU=='P'&
                           simp100$year=='16']) # both have 8 reps ok



# don't do log transform for P (in which nothing changes)
Casimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='Ca2',], na.action=na.omit)
qqr(Casimp100.lme) # some way off
summary(Casimp100.lme) # increase in E, N and P NOT diff with diff slopes 
# year term for N or P not signif



mr2=mr2[mr2$stand %in% c('BO.E','BO.P','Vg.E','Vg.N','Eu.E2','Eu.N',
                         'JP.E1','JP.N','JP.E2','JP.P','It.E1','It.N'),]
#mr2$chgln[mr2$biome=='Cer'&mr2$LU=='N'&mr2$depth=='0-100']=NA
# remove cerrado natveg? this messes up the plot
# only if varwidth=T: can't handle different numbers of obs 
#   in paired boxes
# now works with position.dodge2
# remove pairs together?
mr2$chgln[mr2$site2=='JP2'&mr2$depth=='0-100']=NA
mr2$chgln[mr2$site2=='It'&mr2$depth=='0-100']=NA
mr2=group_by(mr2,element,LU,depth) %>%
  mutate(chglnmn = mean(chgln,na.rm=T),chglnmax=max(chgln,na.rm=T),
         chglnmin=min(chgln,na.rm=T), newnobs=n())
mr2=group_by(mr2,element) %>%
  mutate(maxchgln=max(chgln,na.rm=T),minchgln=min(chgln,na.rm=T))

mr2=mutate(mr2,sigyr=rep(NA),sigveg=rep(NA))
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU=='E']='+'
# Not with different slopes for different groups: p= .1
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU=='P']='-'
mr2$sigyr[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='E']='+'
#mr2$sigyr[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='P']='+'
mr2$sigveg[mr2$depth=='0-20'& mr2$element=='Ca' & mr2$LU=='N']='*'
mr2$sigveg[mr2$depth=='0-20'& mr2$element=='C' & mr2$LU!='E']='*'
#mr2$sigyr[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU!='E']='-' nope
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU=='P']='-'
#mr2$sigyr[mr2$depth=='0-100'& mr2$element=='N' & mr2$LU=='E']='+' nope
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='K' & mr2$LU=='N']='*'
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='E']='+'
#mr2$sigveg[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='N']='*' nope
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU=='P']='*'

mr2$depth=factor(mr2$depth,levels=c('0-20','0-100'))

median(exp(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='E']))
mean(exp(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='E']))
exp(median(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='E']))
exp(mean(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='E']))
exp(.159) #estimate for year effect from lme
# super close to exp(mean), ok
exp(mean(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='N']))
exp(mean(mr2$chgln[mr2$element=='C'&mr2$depth=='0-20'&mr2$LU=='P']))
exp(mean(mr2$chgln[mr2$element=='Ca'&mr2$depth=='0-20'&mr2$LU=='P']))
exp(mean(mr2$chgln[mr2$element=='Ca'&mr2$depth=='0-20'&mr2$LU=='E']))
# From lm table instead:
exp(1.029) # C in euc
exp(1.029)*exp(-.775) # in pasture


# ggplot version

ggplot(data=mr2, aes(x=LUlongish, y=chgln,width=newnobs,#y=chgln20,
                     #fill=depth,colour=LUlongish))+
                     fill=LUlongish,colour=depth))+
  #colour=LUlongish,shape=depth))+
  #colour='grey30')) + 
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_boxplot(varwidth =T,size=.8,#show.legend = F,
               position=position_dodge2(.8,preserve='single'),
               outlier.shape = 1,width=.8) +
  #geom_pointrange(size=3,
  #             position=position_dodge2(1,preserve='single')) +
  #geom_point(aes(x=LUlongish,y=chglnmn,colour=LUlongish),
  #           shape=15,size=2,show.legend = F)+
  #geom_point(aes(x=LUlongish,y=chglnmn,
  #               colour=depth),
  #shape=depth),
  #           size=2,show.legend = F,shape=19,
  #           position=position_dodge(.8))+
  scale_colour_manual(values=c('lightblue','steelblue'),
                      name='Depth (cm)')+
  scale_fill_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                    name='Vegetation type',
                    labels=c('Eucalyptus','Native vegetation',
                             'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  scale_y_continuous(sec.axis = 
                       sec_axis(trans=~.,name='Percent change in stock',
                                breaks=pct_to_L(c(-75,-25,25,75,200, 1000)),
                                labels=paste(c('-75', '-25','+25','+75',
                                               '+200', '+1000'),'%',sep='')))+
  #facet_wrap(depth~element,ncol=3,scales='free_y') +
  labs(y='Log change in stock over 12 years', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())+ 
  #scale_fill_discrete(labels=c('20-100 cm','0-20 cm'),
  #                    guide = guide_legend(reverse=TRUE,title=NULL))+
  #geom_point(mapping = aes(y = hibar+1, shape = as.factor(sig)),show.legend=F)
  geom_text(mapping = aes(x=LUlongish, y=maxchgln*0.9*!is.na(sigveg),
                          label = sigveg,colour=depth),
            size=8,na.rm=T,show.legend=F,#colour='black',
            position=position_dodge(.8))+
  geom_text(mapping = aes(x=LUlongish, y=(maxchgln)*!is.na(sigyr), 
                          label = sigyr,colour=depth),
            size=6,na.rm=T,show.legend=F,#colour='black',
            position=position_dodge(.8))

mr2=mr2[order(mr2$LU),]
# plot the points instead of making boxplots:
ggplot(data=mr2, aes(x=depth, y=chgln,#shape=depth,
                     shape=LUlongish,colour=LUlongish))+
  scale_x_discrete()+
  scale_colour_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                      name='Vegetation type',
                      labels=c('Eucalyptus','Native vegetation',
                               'Pasture'))+
  #scale_shape_manual(values=c(1,0),
  #                   name='Depth (cm)')+
  scale_shape_manual(values=c(1,5,0),
                     name='Vegetation type',
                     labels=c('Eucalyptus','Native vegetation',
                              'Pasture'))+
  
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_point(size=2.5, position=position_dodge(.8))+
  #geom_point(aes(x=LUlongish,y=chglnmn,colour=LUlongish),
  #           data=distinct(mr2,LUlongish,element,depth,.keep_all=T),
  #           shape=18,size=3,show.legend = F,position=position_dodge2(.8))+
            #without subsetting data, print all the points to make a sort of bar
            #shape=15,size=1,show.legend = F,position=position_dodge2(.8))+
  geom_pointrange(aes(x=depth,y=chglnmn,colour=LUlongish,ymax=chglnmax,
                      ymin=chglnmin),
             data=distinct(mr2,LUlongish,element,depth,.keep_all=T),
             shape=18,show.legend = F,position=position_dodge2(.8))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  scale_y_continuous(sec.axis = 
                       sec_axis(trans=~.,name='Percent change in stock',
                                breaks=pct_to_L(c(-75,-25,25,75,200, 1000)),
                                labels=paste(c('-75', '-25','+25','+75',
                                               '+200', '+1000'),'%',sep='')))+
  labs(y='ln(stock in 2016 / stock in 2004)', x="Depth (cm)") +
  #guides(colour=guide_legend(),shape=F)+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_line(colour='grey80'),
        legend.key=element_blank(),
        legend.title=element_text(size=10))+ 
  geom_text(mapping = aes( y=maxchgln*1.15*!is.na(sigveg),#x=LUlongish,
                           x=depth,#colour=LUlongish,
                          label = sigveg),
            size=6,na.rm=T,show.legend=F,colour='black',
            position=position_dodge(.8))+
  geom_text(mapping = aes(y=(maxchgln)*1.15*!is.na(sigyr),#x=LUlongish, 
                          x=depth,#colour=LUlongish,
                          label = sigyr),
            size=5,na.rm=T,show.legend=F, colour='black',
            position=position_dodge(.8))
# Too busy?



# old figure: lattice
bwplot(chgln20~LUlongish|element,ylim=c(-1,1),las=1,
       data=mrstkchgslim[mrstkchgslim$stand %in% 
                           c('BO.E','BO.P','Vg.E','Vg.N','Eu.E2','Eu.N',
                             'JP.E1','JP.N','JP.E2','JP.P','It.E1','It.N'),],
       as.table=T,col.line='black',box.ratio=0,#varwidth=T,
       ylab='Log change in stock over 12 years, 0-20 cm',
       layout=c(3,2),
       
       panel = function(x, y, ...){
         panel.bwplot(x, y, ...)
         superpose.symbol=list(fill=c('blue3','springgreen','darkgoldenrod1'))
         superpose.point=list(col=c('blue3','springgreen','darkgoldenrod1'))
         panel.axis(side=ifelse(panel.number()==3,'right','left'),
                    at=pct_to_L(c(-50,-25,25,75,125)),
                    labels=paste(c('-50','-25','+25','+75','+125'),
                                 '%',sep=''),outside = F,half=F,
                    # I want outside=T but it won't draw anything
                    draw.labels=ifelse(panel.number()%in%c(3,4),T,F),
                    ticks=ifelse(panel.number()%in%c(3,4),T,F))
         panel.abline(h=0,lty=3)
       })





simple20_2$LU=factor(simple20_2$LU,levels=c('E','N','P'))
simple20_2lim=simple20_2[simple20_2$element %in% c('C','N','K','P','Ca'),]

ggplot(data=simple20_2lim, aes(x=LU, y=stock20,
                     fill=LU,colour=year))+
  geom_boxplot(varwidth =T,size=.8,#show.legend = F,
               position=position_dodge2(.8,preserve='single'),
               outlier.shape = 1) +
  scale_colour_manual(values=c('red','blue'),
                      name='Year',labels=c('2004','2016'))+
  scale_fill_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                    name='Vegetation type',
                    labels=c('Eucalyptus','Native vegetation',
                             'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  labs(y='Stock, 0-20 cm, Mg/ha', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())

simple20_2lim=simple20_2[simple20_2$element %in% c('C','N','K','P','Ca'),]
# Make a new column to show which euc stands are paired with native vs pasture
simple20_2lim=mutate(simple20_2lim,LU3=as.character(LU))
simple20_2lim$LU3[simple20_2lim$stand=='JP.E1']='EP'
simple20_2lim$LU3[simple20_2lim$LU3=='E']='EN'
simple20_2lim$LU3=factor(simple20_2lim$LU3,levels=c('EN','N','EP','P'))
simple20_2lim$element =factor(simple20_2lim$element,
                              levels=c('C','N','K','P','Ca'))

ggplot(data=simple20_2lim, aes(x=LU3, y=stock20,
                               fill=LU,colour=year))+
  geom_boxplot(varwidth =T,size=.8,#show.legend = F,
               position=position_dodge2(.8,preserve='single'),
               outlier.shape = 1) +
  scale_colour_manual(values=c('red','blue'),
                      name='Year',labels=c('2004','2016'))+
  scale_fill_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                    name='Vegetation type',
                    labels=c('Eucalyptus','Native vegetation',
                             'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  labs(y='Stock, 0-20 cm, Mg/ha', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())

# Simpler version with less data?
simple20_2limd=group_by(simple20_2lim,element,year,LU3) %>%
  mutate(mn20=mean(stock20,na.rm=T),nobs=n(),
         se20=sd(stock20,na.rm=T)/sqrt(nobs),
         min20=min(stock20,na.rm=T),
         max20=max(stock20,na.rm=T))
simple20_2limd=distinct(simple20_2limd,element,year,LU,LU3,.keep_all=T)

ggplot(data=simple20_2limd, aes(x=LU3, y=mn20,colour=LU,shape=year))+
  geom_vline(xintercept=7.5)+
  geom_pointrange(aes(ymin=min20,ymax=max20),size=.8,#show.legend = F
  #geom_pointrange(aes(ymin=mn20-se20,ymax=mn20+se20),size=.8,#show.legend = F,
                                  position=position_dodge2(.8,preserve='single')) +
  scale_x_discrete()+
  scale_shape_manual(values=c(18,15),
                      name='Year',labels=c('2004','2016'))+
  scale_colour_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                    name='Vegetation type',
                    labels=c('Eucalyptus','Native vegetation',
                             'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  labs(y='Stock, 0-20 cm, Mg/ha', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())

# One point per stand?
simple20_2lime=group_by(simple20_2lim,element,year,stand,LU3) %>%
  mutate(mn20=mean(stock20,na.rm=T),nobs=n(),
         se20=sd(stock20,na.rm=T)/sqrt(nobs),
         min20=min(stock20,na.rm=T),
         max20=max(stock20,na.rm=T))
simple20_2lime=distinct(simple20_2lime,element,year,LU,LU3,.keep_all=T)


ggplot(data=simple20_2lime, aes(x=LU3, y=mn20,colour=LU,shape=year))+
  geom_point(size=3, position=position_dodge(.8))+
  scale_x_discrete()+
  scale_shape_manual(values=c(18,15),
                     name='Year',labels=c('2004','2016'))+
  scale_colour_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                      name='Vegetation type',
                      labels=c('Eucalyptus','Native vegetation',
                               'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  labs(y='Stock, 0-20 cm, Mg/ha', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())


simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]
simp100limd=group_by(simple20_2lim[simple20_2lim$site2!='JP2' & 
                                     simple20_2lim$site2!='It',],
                     element,year,LU3) %>%
  mutate(mn100=mean(stock100,na.rm=T),nobs=n(),
         se100=sd(stock100,na.rm=T)/sqrt(nobs),
         min100=min(stock100,na.rm=T),
         max100=max(stock100,na.rm=T))
simp100limd=distinct(simp100limd,element,year,LU,LU3,.keep_all=T)

ggplot(data=simp100limd, aes(x=LU3, y=mn100,colour=LU,shape=year))+
  geom_vline(xintercept=7.5)+
  geom_pointrange(aes(ymin=min100,ymax=max100),size=.8,#show.legend = F
                  #geom_pointrange(aes(ymin=mn20-se20,ymax=mn20+se20),size=.8,#show.legend = F,
                  position=position_dodge2(.8,preserve='single')) +
  scale_x_discrete()+
  scale_shape_manual(values=c(18,15),
                     name='Year',labels=c('2004','2016'))+
  scale_colour_manual(values=c('blue3','springgreen','darkgoldenrod1'),
                      name='Vegetation type',
                      labels=c('Eucalyptus','Native vegetation',
                               'Pasture'))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  labs(y='Stock, 0-100 cm, Mg/ha', x="") +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(1.2,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank())



Krat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='K',],
                        na.action = na.omit,family='binomial')
qqr(Krat100simp.pql) # ok? upper tail off
summary(Krat100simp.pql) # ratio starts bigger in noneuc (p=.053)
#   and decreases (p=.027)
# without It and JP2 (i.e. native Cerrado), no significant changes
# Native (just AF) almost higher with 3 separate LUs (p=.097)

Crat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='C',],
                        na.action = na.omit,family='binomial')
qqr(Crat100simp.pql) # ok
summary(Crat100simp.pql) # no change

Nrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='N',],
                        na.action = na.omit,family='binomial')
qqr(Nrat100simp.pql) # tails quite off; ok without Cerr nat
summary(Nrat100simp.pql) # with JP2 and It, decreases overall
# without, no change

Carat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='Ca',],
                        na.action = na.omit,family='binomial')
qqr(Carat100simp.pql) #good
summary(Carat100simp.pql) # marginal increase in euc, decrease in non
# starts higher in non
# without Cerr nat, same deal, but changes are signif only for pasture
# (marginally higher starting value p=.083 and decrease p=.067 in nat)
tapply(simp100$stockratio[simp100$element=='Ca'&simp100$year=='04'],
       simp100$LU[simp100$element=='Ca'&simp100$year=='04'],summary)
tapply(simp100$stockratio[simp100$element=='Ca'&simp100$year=='16'],
       simp100$LU[simp100$element=='Ca'&simp100$year=='16'],summary)

Prat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                         data=simp100[simp100$element=='P2',],
                         na.action = na.omit,family='binomial')
qqr(Prat100simp.pql) 
summary(Prat100simp.pql) # shallower under non-euc (if using P not P2), no change

Alrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                        data=simp100[simp100$element=='Al',],
                        na.action = na.omit,family='quasibinomial')
qqr(Alrat100simp.pql) 
summary(Alrat100simp.pql) # gets deeper in euc but not in native, shallower in pasture

Zrrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                         data=simp100[simp100$element=='Zr',],
                         na.action = na.omit,family='quasibinomial')
qqr(Zrrat100simp.pql) 
summary(Zrrat100simp.pql) # no significant anything, good
Znrat100simp.pql=glmmPQL(stockratio~year*LU,random=~1|site,
                         data=simp100[simp100$element=='Zn',],
                         na.action = na.omit,family='quasibinomial')
qqr(Znrat100simp.pql) # 3 outliers at upper tail
summary(Znrat100simp.pql) # also nothing
cor(widedats4$C,widedats4$Al,method='pear',use='pair') #.38
# why is this positive? because lots of C and Al in Vg and It?
# because C is on clay
AlC.lme=lme(log(C)~Al,random=~1|site/stand,na.action=na.omit,
            data=widedats4[widedats4$site!='TM'& widedats4$site!='Cr'&
                             widedats4$LU!='A',])
qqr(AlC.lme) # pretty good with log, some outliers
summary(AlC.lme) # yeah, negative, p < 10-4


yrdiffstockplot20_LU(tstock[tstock$element=='C'&tstock$stand %in%
                              unique(simple20_2$stand),])
points(stock20_16~stock20_04,pch=16,col='white',
       data=tstock[tstock$element=='C'&tstock$stand %in% c('JP.E2','JP.N'),])
