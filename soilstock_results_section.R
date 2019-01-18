# Soil stock analysis

#setwd('C:\\Users\\Devin\\Documents\\Soil data')
source('soil_data_reader.R')

# Figure 2
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


# Year trends between stocks in each biome (Figure 3)

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
#eucC20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
#                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
#summary(eucC20bm.lme2) # increase in Cerrado now barely signif (p=.046)
#qqr(eucC20bm.lme2) # slightly better?
#anova(eucC20bm.lme,eucC20bm.lme2) # not different

# N: doesn't converge with a random slope (does with optim)
eucN20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) # no change
qqr(eucN20bm.lme) # mostly ok

#eucN20bm.lme2=lme(stock20~year*biome,random=~1+year|site/stand,
#                  data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
#summary(eucN20bm.lme2) # no change
#qqr(eucN20bm.lme2) # tails off
#anova(eucN20bm.lme,eucN20bm.lme2)
# 2 has lower AIC but higher BIC; differ at p=.02 (no longer; 2 is worse)

eucP20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
qqr(eucP20bm.lme) # tails off without Eu, but much less bad
summary(eucP20bm.lme) # no change
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
#eucCa20bm.lme2=lme(log(stock20)~year*biome,random=~1+year|site/stand,
#                 data=euc2deps[euc2deps$element=='Ca2',],
#                 na.action = na.omit,control = lmeControl(opt='optim'))
#qqr(eucCa20bm.lme2)
#summary(eucCa20bm.lme2) 
#anova(eucCa20bm.lme,eucCa20bm.lme2)
# with a random slope, no signif effects
# but that seems wrong, clearly a lot of Ca was added in the Cerrado
# although not consistently across cerrado stands, hence loss of significance?
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


# 0-100 cm
eucCN100bm.lme=lme(conc100~year*biome,random=~1|site/stand,
                  data=euc2deps2[euc2deps2$element=='CN',],na.action = na.omit)
summary(eucCN100bm.lme) # decreases in AF, increases in Cer
qqr(eucCN100bm.lme) # one point off at each tail, upper far off
plot(eucCN100bm.lme)

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
summary(eucP100bm.lme) 
eucP100bm.lme2=lme(stock100~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2'&
                                   euc2deps2$site!='Eu',],
                  random=~1|site/stand,na.action=na.omit)
qqr(eucP100bm.lme2) # one outlier each end; almost signif in AF w/o Eu
summary(eucP100bm.lme2)

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

#Mean stocks
table(shorttstk$stockunit,shorttstk$element)
nutstk=group_by(dats2deps[dats2deps$element %in% c('C','N','K','P2','Ca2') &
                            dats2deps$stand !='It.E1',],
                element,LU,biome,year)
nutstk$stock100[nutstk$site=='Bp'&nutstk$element=='K']=NA
nutstk$stock20[nutstk$site=='Bp'&nutstk$element=='K']=NA
nutstk= summarise(nutstk,stk20down=mean(stock100-stock20,na.rm=T),
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
nutstkE=nutstk2[nutstk2$LU=='E',]
#plot(stock~element,data=nutstkE)
#library(ggplot2)
nutstkE <- with(nutstkE, nutstkE[order(element,year,depth),])
nutstkE=mutate(nutstkE,sig=rep(NA))
nutstkE$sig[nutstkE$depth=='20down'&nutstkE$element=='N'&nutstkE$biome=='AF'&
              nutstkE$year=='16']=1
nutstkE$sig[nutstkE$element=='Ca'& nutstkE$biome=='Cer'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig[nutstkE$element=='C'&nutstkE$biome=='Cer'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig[nutstkE$element=='K'&nutstkE$biome=='AF'&
              nutstkE$year=='16']=1 # both depths
nutstkE$sig=as.integer(nutstkE$sig)
labls <- c(AF = "Atlantic Forest", Cer = "Cerrado")

ggplot(data=nutstkE, aes(x=year, y=stock, fill=depth)) +
  geom_bar(stat="identity") + 
  facet_grid(element~biome,labeller = labeller(biome=labls)) + 
  coord_flip() +
  labs(y="Stock (10 Mg/ha for C; Mg/ha for other elements)",
       x="Year", fill="Depth") +
  theme(strip.text.y = element_text(angle = 0),
        legend.position=c(0.8,0.1),
        panel.background = element_rect(fill='white'),
        panel.grid.major.x = element_line(colour='grey80'),
        panel.grid.major.y = element_blank()) +
  geom_errorbar(aes(ymax=hibar,  ymin=lobar), width=0.15) +
  scale_fill_manual(values=c('slategray','lightblue'),
                      labels=c('20-100 cm','0-20 cm'),
                      guide = guide_legend(reverse=TRUE,title=NULL))+
  geom_point(mapping = aes(y = (hibar+1)*(sig>0)),
             shape=18,size=3,show.legend=F)
  #geom_text(mapping = aes(y = (hibar+1)*(sig>0)), label = '*',
  #          na.rm=T,show.legend=F,size=8,hjust=0)


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
       shortE$element[shorttstk$LU=='E'],
       function(x){mean(x,na.rm=T)})
tapply(shorttstk$stock20_04[shorttstk$LU=='E'],
       shortE$element[shorttstk$LU=='E'],
       function(x){mean(x,na.rm=T)})

tapply(shortE$rat_16,shortE$element,
       function(x){mean(x,na.rm=T)})
tapply(shortE$rat_16,shortE$element,
       function(x){sd(x,na.rm=T)})
mean(shortE$BD100_16[shortE$element=='C'])
sd(shortE$BD100_16[shortE$element=='C'])

tapply(shortE$stock100_16[shortE$element=='C'],
       shortE$biome[shortE$element=='C'],
       function(x){mean(x,na.rm=T)}) 

mean(shorttstk$stock20_04[shorttstk$LU=='E' & shorttstk$element=='C']) # 46.0
mean(shorttstk$stock20_16[shorttstk$LU=='E' & shorttstk$element=='C']) # 50.1


# just-euc ratios
# These are proportion data and don't come from normal distributions
# Try generalized mixed models
#library(MASS)
Krateuc2.pql=glmmPQL(stockratio~year*biome,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='K' &
                                      euc2deps2$site!='Bp',],
                     na.action = na.omit,family='quasibinomial')
summary(Krateuc2.pql) # same deal as lme, increases at p=.07
qqr(Krateuc2.pql)
# with Bp included, same result, somewhat smaller increase

Nrateuc2.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                     data=euc2deps2[euc2deps2$element=='N',],
                     na.action = na.omit,family='quasibinomial')
summary(Nrateuc2.pql) # again, same as lme: decreases, p=.026
qqr(Nrateuc2.pql) #ok

Crateuc.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                    data=euc2deps2[euc2deps2$element=='C',],
                    na.action = na.omit,family='quasibinomial')
summary(Crateuc.pql) # nada
qqr(Crateuc.pql) # not great
Crateuc.pql=glmmPQL(stockratio~year,random=~1|site/stand,
                    data=euc2deps2[euc2deps2$element=='C',],
                    na.action = na.omit,family='quasibinomial')
summary(Crateuc.pql) # nada
qqr(Crateuc.pql) # not great

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
             data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                         dats$element =='N',],na.action=na.omit) # log makes qqr worse
summary(Nelt.lme) # L lower again at p=.06
Neltbm.lme=lme(value~elt*biome, random=~1|stand,
               data=dats[dats$depth==5 & dats$year=='16' & dats$elt %in% c('E','L','T') &
                           dats$element =='N',],na.action=na.omit)
summary(Neltbm.lme) # L < E for Cerrado at p=.02; T < E at .06


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
# max N export = 478 kg ha-1 over the 12 years

summary(I((stkchgs$maxbudgconc-stkchgs$minbudgconc)/stkchgs$budget))
# -10.6 to +4.3 -- times, not %. Mean -5%, median + 45%

plot(Concentration~Egrandconc,data=budgets,pch=as.character(Nutrient),col=Nutrient)
abline(0,1)
# K and N vary a lot; one very low Ca value = BO.E (no bark measured)
# My measured N is almost always lower than Pagano estimate
# So if using their estimate, even larger N decreases would be predicted
# Plantar data has much less N, more K than my estimates
plot(Concentration~SantanaMG,data=budgets,pch=as.character(Nutrient),col=Nutrient)
# More N and Ca than Egrand; more of everything than my measurements
budgets$Stand[budgets$Concentration<.0008 & budgets$Nutrient=='Ca']
plot(budget~plconcbudg,data=budgets,pch=as.character(Nutrient),col=Nutrient)
abline(0,1)
# N and sometimes K influence the budget a lot.
plot(denserbudg~lessdensebudg,data=budgets,
     pch=as.character(Nutrient),col=Nutrient)
abline(0,1) # pretty close; only matters for C and N

# variation in percent bark could also be important
plot(bark20budg~bark5budg,data=budgets,
     pch=as.character(Nutrient),col=Nutrient)
abline(0,1)
abline(v=0,lty=3)
abline(h=0,lty=3)# but not that important? mostly matters for K, sometimes Ca

tapply(stkchgs$chgrt20,stkchgs$element,mean)
tapply(stkchgs$chgln20,stkchgs$element,mean)
tapply(stkchgs$budget,stkchgs$element,mean)
tapply(stkchgs$efs20,stkchgs$element,mean)
log(abs(tapply(stkchgs$budget,stkchgs$element,mean)))

t.test(stkchgs2$chg20,stkchgs2$budget,paired=T) 
t.test(stkchgs$chg20[stkchgs$element=='N'],
       stkchgs$budget[stkchgs$element=='N'],paired=T)
t.test(stkchgs$chg20,stkchgs$lessrotbudg,paired=T) 
t.test(stkchgs$chg20,stkchgs$woodonlybudg,paired=T) # p=.059
# 
t.test(stkchgs$bark20budg,stkchgs$woodonlybudg,paired=T) 
# those are different, good

palette(rainbow(9))
plot(chg100~budget,data=stkchgs,type='n', 
     xlab='Fertilizer - harvest, Mg ha-1',
     ylab='Observed change in stocks to 100 cm',las=1)
rect(xleft=-1, ybottom=-1, xright=2, ytop=2,border='gray50')
abline(0,1)
abline(h=0,lty=3)
abline(v=0,lty=3)
text(stkchgs$budget,stkchgs$chg100,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))


segments(x0=stkchgs$bark20budg,x1=stkchgs$woodonlybudg,y0=stkchgs$chg20,
         col=as.numeric(stkchgs$stand))


plot(chg20~budget,data=stkchgs,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 20 cm, Mg ha-1',
     #ylim=c(-1,1.8),
     las=1)
rect(xleft=-.2, ybottom=-.2, xright=.5, ytop=.5,border='gray50')
segments(x0=stkchgs$minbudgconc,x1=stkchgs$maxbudgconc,y0=stkchgs$chg20,
         col=as.numeric(stkchgs$stand))
# Much cleaner with conc
segments(x0=stkchgs$budget,y0=stkchgs$chg20-stkchgs$sdchg20,
         y1=stkchgs$chg20+stkchgs$sdchg20,
         col=as.numeric(stkchgs$stand))
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
#text(stkchgs$budget,stkchgs$chg20,labels=stkchgs$element,
text(stkchgs$plconcbudg,stkchgs$chg20,labels=stkchgs$element,
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
# Or maybe they left the bark onsite in JP so increases > expected? No
# Changing the concentrations to the literature values I cited for E. grandis
#   improves Ca in BO.E, but doesn't help much with other nutrients

# taking out the second rotation helps for Ca, K in BO.E and Bp.E1, not for JPs or N
stkchgs3=stkchgs[stkchgs$element!='Mg',]

plot(chg20~budget,data=stkchgs3,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     #xlim=c(-.2,.5),ylim=c(-.2,.5),
     ylab='Observed change in stocks to 20 cm, Mg ha-1',las=1)
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
segments(x0=stkchgs3$minbudgconc,x1=stkchgs3$maxbudgconc,y0=stkchgs3$chg20,
         col=as.factor(stkchgs3$stand))
segments(x0=stkchgs3$budget,y0=stkchgs3$chg20-stkchgs3$sdchg20,
         y1=stkchgs3$chg20+stkchgs3$sdchg20,
         col=as.factor(stkchgs3$stand))
text(stkchgs3$budget,stkchgs3$chg20,labels=stkchgs3$element,
#text(stkchgs3$lessrotbudg,stkchgs3$chg20,labels=stkchgs3$element,
          cex=stkchgs3$conc*2000,
     col=as.numeric(stkchgs3$stand))
# What is a realistic range of bark? How to present sensitivities?
# Table of ratios of budget to its variations?

JPchg=stkchgs3[stkchgs3$stand %in% c('JP.E1','JP.E2'),]
#text(stkchgs3$lessrotbudget,stkchgs3$chg20,labels=stkchgs3$element,
text(stkchgs3$moreharvbudg,stkchgs3$chg20,labels=stkchgs3$element,
     cex=stkchgs3$conc*2000,
     col=as.numeric(stkchgs3$stand))


stkchgs3=stkchgs[stkchgs$element!='Mg',]

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


plot(chg100~budget,data=stkchgs,type='n', 
     xlab='Net nutrient input (fertilizer - harvest), Mg ha-1',
     ylab='Observed change in stocks to 100 cm, Mg ha-1',las=1,
     ylim=c(-3,5.5))
     #xlim=c(-.2,.5),ylim=c(-.2,.5))
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(0,1)
#segments(x0=stkchgs$minbudg,x1=stkchgs$maxbudg,y0=stkchgs$chg20,
#         col=as.numeric(as.factor(stkchgs$biome))+2)
segments(x0=stkchgs$minbudgconc,x1=stkchgs$maxbudgconc,y0=stkchgs$chg20,
         col=as.numeric(as.factor(stkchgs$biome))+2)
segments(x0=stkchgs$budget,y0=stkchgs$chg100-stkchgs$sdchg100,
         y1=stkchgs$chg100+stkchgs$sdchg100,
         col=as.numeric(as.factor(stkchgs$biome))+2)
text(stkchgs$budget,stkchgs$chg100,labels=stkchgs$element,
     col=as.numeric(as.factor(stkchgs$biome))+2)
legend('bottomright',pch=15,col=c(3,4),
       legend=c('Atlantic Forest','Cerrado'),bty='n')

segments(x0=stkchgs$minbudg,x1=stkchgs$maxbudg,y0=stkchgs$chg100,
         col=as.factor(stkchgs$stand))
segments(x0=stkchgs$budget,y0=stkchgs$chg100-stkchgs$sdchg100,
         y1=stkchgs$chg100+stkchgs$sdchg100,
         col=as.factor(stkchgs$stand))
text(stkchgs$budget,stkchgs$chg100,labels=stkchgs$element,
     col=as.numeric(stkchgs$stand))




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
qqr(Ccaov) # horrible
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


table(simple20$LU2[simple20$element=='C'],simple20$year[simple20$element=='C']) 
# Not balanced--rep 5



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
summary(Zrsimp.lme4) #marginal year-native intrxn


CNsimp.lme=lme(conc20~year*LU,random=~1+year|site2,
  data=simple20_2[simple20_2$element=='CN',], na.action=na.omit)
summary(CNsimp.lme) # increases in euc at p=.046
# signif intrxn with native and pasture 
#   but the decreases in those veg types aren't significantly nonzero
qqr(CNsimp.lme) # pretty ok

simp100=simple20_2[simple20_2$site2!='JP2' & simple20_2$site2!='It',]

Ksimp100.lme=lme(log(stock100)~year*LU,random=~1|site,
                 data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme) # ok
summary(Ksimp100.lme) # with random intercept only, euc increases, others not different
Ksimp100.lme2=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='K',], na.action=na.omit)
qqr(Ksimp100.lme2) # mostly ok
summary(Ksimp100.lme2) # with random slopes, no change for euc, native decr
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
summary(mrstkchgs$chgln100[mrstkchgs$element=='K']) 
# median change is positive, change in medians is negative
simp100$stand[simp100$element=='K'&is.na(simp100$stock100)]
# several. What's up with that? Two BO.E, one JP.E1, one Vg.E, two Vg.N
# also missing values there for P, fewer for N

Nsimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='N',], na.action=na.omit)
qqr(Nsimp100.lme) # ok
summary(Nsimp100.lme)
#With random slope, nothing significant 
#(p for N-yr intrxn = .045); also for P intrxn with N as main

Csimp100.lme=lme(log(stock100)~year*LU,random=~1+year|site,
                 data=simp100[simp100$element=='C',], na.action=na.omit)
qqr(Csimp100.lme)
summary(Csimp100.lme) 
# with random ~1|site, increase in euc isn't signif, 
#   but decreases in native and pasture are! 
# with random ~1|year+site (more appropriate),
#   only signif thing is decrease in pasture (in native, p=.08)
# yes, when pasture or native set as default level, signif decrease
#   also signif intrxn with euc for both

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
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU!='E']='-'
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='N' & mr2$LU=='E']='+'
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='K' & mr2$LU=='N']='*'
mr2$sigyr[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='E']='+'
#mr2$sigveg[mr2$depth=='0-100'& mr2$element=='Ca' & mr2$LU=='N']='*' nope
mr2$sigveg[mr2$depth=='0-100'& mr2$element=='C' & mr2$LU=='P']='*'

mr2$depth=factor(mr2$depth,levels=c('0-20','0-100'))

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
