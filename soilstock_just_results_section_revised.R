# Soil stock analysis

#setwd('C:\\Users\\Devin\\Documents\\Soil data')
source('soil_data_reader.R')

# reviewer 2 requested table of concentrations
concsum=group_by(dats4[dats4$element %in% c('C','N','P2','K','Ca2'),],
                 element,depth) %>%
  summarise(minval=min(value,na.rm=T),maxval=max(value,na.rm=T),
            meanval=mean(value,na.rm=T),medval=median(value,na.rm=T))
#write.csv(concsum,'conc_by_depth.csv')

standnums=data.frame(stand=c('Bp.E1','Bp.E2','JP.E1','JP.P','JP.E2','JP.N',
                             'It.E1','It.N','It.E2','BO.E','BO.P','Vg.E',
                             'Vg.N','Eu.E1','Eu.E2','Eu.N'),
                     num=c(1,2,3,3,4,4,5,5,6,7,7,8,8,9,10,10))
dats=merge(dats,standnums,by='stand')

clay2 = data.frame(stand=c('Eu.E2','Eu.E1','Eu.N','BO.E','BO.P', 
                           'Vg.E','Vg.N','Bp.E1','Bp.E2',#not sure which Bp is which
                           'It.E1','It.E2','It.N',
                           'JP.E1','JP.E2','JP.N','JP.P'),
                   clay5=c(15,18,12,58,61,
                           74,57, 85,88,
                           77,65,74, 15,16,7,17),
                   clay15=c(16,21,14,61,61,
                            64,70, 78,88,
                            75,67,71, 16,18,6,17),
                   clay80=c(35,41,41,72,71,
                            74,58, 77,87,
                            79,74,76, 18,22,8,20))
clay2=mutate(clay2,clay20=(clay5+clay15)/2)
clay2$clay20

# Figure S3
# Highlighting the scale of the spatial heterogeneity in Bps
png('fig_s3_was_s2.png',res=100,height=5,width=6,units='in')
ggplot(data=dats[dats$element=='K' & dats$site=='Bp',],
       aes(x=value/1000,y=depth, color=as.factor(rep)))+
  geom_point(shape=16,size=2.5,show.legend=F)+ 
  scale_y_reverse()+
  facet_wrap(year~stand)+
  theme(panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank())+
  labs(y='Depth (cm)',x='Soil K (g/kg)')+
  scale_color_brewer(palette='Dark2')
dev.off()

# Figure S4
#########
tstock=merge(tstock,standnums,by='stand')
yrdiffstockplot20_LUnum=function(sub_ttests,label=T){
  #par(mar=c(5,5,2,2))
  par(mar=c(4,4,.5,.5))
  #palette(c('darkgoldenrod1','blue3','springgreen'))
  palette(c('#e79f00','#000000','#0072B2'))
  sub_ttests$LU=factor(sub_ttests$LU,levels=c('P','E','N'))
  sub_ttests$site=factor(sub_ttests$site,
                         levels=c('BO','Bp','It','JP','Eu','Vg'))
  sub_ttests$biome=factor(sub_ttests$biome,levels=c('Cer','AF'))
  plot(stock20_16~stock20_04,data=sub_ttests,type='n',las=1,
       xlab='2004',ylab='2016',#cex.lab=1.6,cex.axis=1.5,
       xlim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05),
       ylim=c(min(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*.9,
              max(c(sub_ttests$stock20_04,sub_ttests$stock20_16),na.rm=T)*1.05))
  abline(0,1)
  if(label==T){legend('topleft',bty='n',#cex=1.8,
                      legend=paste(unique(sub_ttests$element),'(Mg / ha)\n0-20 cm',
                                   sep=' '))}
  # Replace std devs by std errors
  segments(sub_ttests$stock20_04,sub_ttests$stock20_16-sub_ttests$sd20_16/sqrt(3),
           sub_ttests$stock20_04,sub_ttests$stock20_16+sub_ttests$sd20_16/sqrt(3),
           col='gray60',lwd=1.5)
  segments(sub_ttests$stock20_04-sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           sub_ttests$stock20_04+sub_ttests$sd20_04/sqrt(3),sub_ttests$stock20_16,
           col='gray60',lwd=1.5)
  text(sub_ttests$stock20_04,sub_ttests$stock20_16,
       col=as.numeric(sub_ttests$LU),labels=sub_ttests$num,cex=1.8)
  palette('default')
  par(mar=c(4,4,2,2))
}

#png('fig_s3_nums2.png',height=6,width=9,units='in',res=300)
par(mfrow=c(2,3))
yrdiffstockplot20_LUnum(tstock[tstock$element=='C',])
legend('bottomright',fill=c('#000000','#0072B2','#e79f00'),bty='n',
       border='white',legend=c('Eucalyptus','Native vegetation','Pasture'))
legend('bottomright',bty='n',legend='a',cex=1.5)
yrdiffstockplot20_LUnum(tstock[tstock$element=='N',])
legend('bottomright',bty='n',legend='b',cex=1.5)
yrdiffstockplot20_LUnum(tstock[tstock$element=='K' &
                                 tstock$site!='Bp',])
legend('bottomright',bty='n',legend='c',cex=1.5)

yrdiffstockplot20_LUnum(tstock[tstock$element=='P2',],label=F)
legend('topleft',bty='n',legend='P (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='d',cex=1.5)

yrdiffstockplot20_LUnum(tstock[tstock$element=='Ca2',],label=F)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='e',cex=1.5)

yrdiffstockplot20_LUnum(tstock[tstock$element=='Ca2' & tstock$stock20_16<1,],label=F)
legend('topleft',bty='n',legend='Ca (Mg / ha)\n0-20 cm')
legend('bottomright',bty='n',legend='f',cex=1.5)
par(mfrow=c(1,1))

dev.off()
#########

# how many samples had detection limit replaced?
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

# Year trends between stocks in eucalyptus stands
# without biome term (excluded from paper to avoid redundancy)
##############################
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
##############################

# Year trends between stocks in each biome (Figure 2)
# lmes are in the lme_tables_printing file
######################
eucC20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                    data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
summary(eucC20bm.lme) # increase in Cerrado only (also if using euc2deps2)
qqr(eucC20bm.lme) # again, log tranform helps a bit, increases significance

eucN20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='N',],na.action = na.omit)
summary(eucN20bm.lme) # no change
qqr(eucN20bm.lme) # mostly ok
intervals(eucN20bm.lme)

eucP20bm.lme=lme(stock20~year*biome,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='P2',],
                 na.action = na.omit)
qqr(eucP20bm.lme) # tails off without Eu, but much less bad
summary(eucP20bm.lme) # no change
intervals(eucP20bm.lme)

eucCa20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='Ca2',],
                  na.action = na.omit)#,control = lmeControl(opt='optim'))
qqr(eucCa20bm.lme)
summary(eucCa20bm.lme) # maybe increase in AF (p=.07), increase in Cer (.0005)
intervals(eucCa20bm.lme)

eucK20bm.lme=lme(log(stock20)~year*biome,random=~1|site/stand,
                  data=euc2deps[euc2deps$element=='K'&
                                  euc2deps$site!='Bp',],na.action = na.omit)
qqr(eucK20bm.lme)
summary(eucK20bm.lme) # increases in AF only
######################

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
# lme intercept estimate: equivalent to mean
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

eucP100bm.lme=lme(log(stock100)~year*biome,
                  data=euc2deps2[euc2deps2$element=='P2',],
                  random=~1|site/stand,na.action=na.omit)
qqr(eucP100bm.lme) # quite bad
summary(eucP100bm.lme) # without log, now significant! increases in AF
# qqr is bad either way; keep log for consistency

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

png('fig2.png',height=6,width=6,units='in',res=500)
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
dev.off()

# Presentation version
png('fig2_por_pres.png',height=6.45,width=6.04,units='in',res=500)
ggplot(data=nutstkE, aes(x=year, y=stock, fill=depth)) +
  geom_bar(stat="identity") + 
  facet_grid(element~biome,labeller = labeller(biome=labls_por)) + 
  coord_flip() +
  labs(y="Estoque (10 Mg/ha para C; Mg/ha para outros)",
       x="Ano", fill="Profundidade") +
  #facet_grid(element~biome,labeller = labeller(biome=labls)) + 
  #coord_flip() +
  #labs(y="Stock (10 Mg/ha for C; Mg/ha for other elements)",
  #     x="Year", fill="Depth") +
  theme(strip.text.y = element_text(angle = 0,size=16),
        strip.text.x = element_text(size=16),
        legend.position=c(0.75,0.14),
        panel.background = element_rect(fill='white'),
        panel.grid.major.x = element_line(colour='grey80'),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,margin = 
                                  margin(t = 20, b = 5,l=20,r=20)),
        legend.text=element_text(size=16)) +
  geom_errorbar(aes(ymax=hibar,  ymin=lobar), width=0.15) +
  scale_fill_manual(values=c('slategray','lightblue'),
                    labels=c('20-100 cm','0-20 cm'),
                    guide = guide_legend(reverse=TRUE,title=NULL))+
  geom_point(mapping = aes(y = (hibar+1)*(sig>0)),
             shape=18,size=3,show.legend=F)
dev.off()

# Portuguese
labls_por <- c(AF = "Mata Atlântica", Cer = "Cerrado")

png('fig2_por.png',height=6,width=6,units='in',res=500)
ggplot(data=nutstkE, aes(x=year, y=stock, fill=depth)) +
  geom_bar(stat="identity") + 
  facet_grid(element~biome,labeller = labeller(biome=labls_por)) + 
  coord_flip() +
  labs(y="Estoque (10 Mg/ha para C; Mg/ha para outros)",
       x="Ano", fill="Profundidade") +
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
             shape=18,size=3,show.legend=F)
dev.off()


shortE = shorttstk[shorttstk$LU=='E' &
                     shorttstk$stand!='It.E1',]

sefun=function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x))-1)}
summary(shorttstk$BD100_16[shorttstk$element=='C'])
summary(shorttstk$BD100_16[shorttstk$element=='C'&
                             shorttstk$LU=='E'])
summary(shorttstk$BD100_16[shorttstk$element=='C'&
                             shorttstk$LU=='N'])

summary(shorttstk$BD20_16[shorttstk$element=='C'&
                             shorttstk$LU=='E'])


sd(shortE$BD100_16[shortE$element=='C'])


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

# most recent rotation
summary(c(budgets$In_kgha_2[budgets$In_kgha_2!=0],
          budgets$In_kgha_1[budgets$In_kgha_2==0])) # by nutrient tho

# summaries by stand
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


# Figure 4
png(filename = 'fig4_newcen_bw.png',width=5,height=5, units='in',res=150)
png(filename = 'fig4_newcen_blue.png',width=5,height=5, units='in',res=150)
ggplot(obsdifm3,aes(x=element,y=disc,color=budgtype))+
  scale_color_manual(values=c("#9ad0f3", "#0072B2"),#c('grey40','black'),#c('darkblue','darkgreen'),
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

# Presentation size
ggplot(obsdifm3,aes(x=element,y=disc,color=budgtype))+
  scale_color_manual(values=c("#9ad0f3", "#0072B2"),#c('grey40','black'),#c('darkblue','darkgreen'),
                     name='Budget',
                     labels=c('Fertilizer - Harvest',
                              'Fertilizer - Harvest +\nInput from biomass change'))+
  guides(color=guide_legend())+
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_point(size=5, shape=1,position=position_dodge(.3),show.legend = F)+
  geom_pointrange(aes(x=element,y=median,colour=budgtype,ymax=max,
                      ymin=min),fatten=6,size=1,
                  data=distinct(obsdifm3,element,budgtype,.keep_all=T),
                  shape=18,show.legend = T,position=position_dodge2(.3))+
  #labs(y='(Observed change in soil stock - Budget) / Budget', x=NULL) +
  #labs(y='(Measured Δ soil stock - Budget) / Budget', x=NULL) +
  labs(y='Measured Δ soil stock / Budget', x=NULL) +
  theme(legend.position=c(0.6,0.15),#c(0.22,0.9),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_line(colour='grey80'),
        legend.key=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=16)) 


# Figure S4
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


# Figure 3
mr2=rbind(mru,mrd)
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
         chglnmin=min(chgln,na.rm=T), newnobs=n(),
         chglnse=sd(chgln,na.rm=T)/sqrt(newnobs-1)) #added 5-23-19
mr2=group_by(mr2,element) %>%
  mutate(maxchgln=max(chgln,na.rm=T),minchgln=min(chgln,na.rm=T),
         maxse=max(chglnmn+chglnse,minse=min(chglnmn-chglnse)))

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




# Katie Lotterhos http://dr-k-lo.blogspot.com/2013/07/a-color-blind-friendly-palette-for-r.html
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")
plot(x = 1:length(cbbPalette), y = rep(1, length(cbbPalette)), pch = 19, cex = 5, 
     col = cbbPalette, yaxt = "n", bty = "n", xaxt = "n", xlab = "", ylab = "")

# ggplot version
mr2=mr2[order(mr2$LU),]

# plot the points instead of making boxplots:
png('fig3_soil.png',res=500,height=6,width=6,units='in')
ggplot(data=mr2, aes(x=depth, y=chgln,#shape=depth,
                     shape=LUlongish,colour=LUlongish))+
  scale_x_discrete()+
  #scale_colour_manual(values=c('blue3','springgreen','darkgoldenrod1'),
  scale_colour_manual(values=c('#000000','#0072B2','#e79f00'),
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
dev.off()
# or just export the png from the screen with dimensions of 600x600?

# Simple version, with +/- 1 SE instead of points
#png('fig3_soil_new.png',res=150,height=6,width=6,units='in')
png('fig3_soil_shapes.png',res=150,height=6,width=6,units='in')
ggplot(data=distinct(mr2,LUlongish,element,depth,.keep_all=T),
       aes(x=depth, y=chgln,colour=LUlongish,shape=LUlongish))+
  scale_x_discrete()+
  scale_colour_manual(values=c('#000000','#0072B2','#e79f00'),
                      name='Vegetation type',
                      labels=c('Eucalyptus','Native vegetation',
                               'Pasture'))+
  scale_shape_manual(values=c(18,17,15),
                     name='Vegetation type',
                     labels=c('Eucalyptus','Native vegetation',
                              'Pasture'))+
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_pointrange(aes(x=depth,y=chglnmn,colour=LUlongish,ymax=chglnmn+chglnse,
                      ymin=chglnmn-chglnse,shape=LUlongish),#shape=18,
                  show.legend = T,position=position_dodge2(.8))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  scale_y_continuous(sec.axis = 
                       sec_axis(trans=~.,name='Percent change in stock',
                                breaks=pct_to_L(c(-75,-25,25,75,200, 1000)),
                                labels=paste(c('-75', '-25','+25','+75',
                                               '+200', '+1000'),'%',sep='')))+
  labs(y='ln(stock in 2016 / stock in 2004)', x="Depth (cm)") +
  guides(colour=guide_legend())+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill='grey80',size=.7),
        legend.position=c(0.9,0.2),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_line(colour='grey80'),
        legend.key=element_blank(),
        legend.title=element_text(size=10))+ 
  geom_text(mapping = aes( y=maxse*1.1*!is.na(sigveg),#x=LUlongish,
                           x=depth,colour=LUlongish,
                           label = sigveg),
            size=6,na.rm=T,show.legend=F,#colour='black',
            position=position_dodge2(.8))+
  geom_text(mapping = aes(y=(maxse)*1.1*!is.na(sigyr),#x=LUlongish, 
                          x=depth,colour=LUlongish,
                          label = sigyr),
            size=5,na.rm=T,show.legend=F, #colour='black',
            position=position_dodge2(.8))
dev.off()


# Big simple version
#mr2sub=mr2[mr2$element %in% c('C','N','Ca'),]
#png('fig3_big_CNCa.png',res=150,height=5,width=10,units='in')

mr2sub=mr2[mr2$element %in% c('C','Ca'),]
png('fig3_big_CCa.png',res=150,height=5,width=8,units='in')

ggplot(data=mr2sub, aes(x=depth, y=chgln,#shape=depth,
                     shape=LUlongish,colour=LUlongish))+
  scale_x_discrete()+
  scale_colour_manual(values=c('#000000','#0072B2','#e79f00'),
                      name=NULL,
                      labels=c('Eucalyptus','Native vegetation',
                               'Pasture'))+
  scale_shape_manual(values=c(1,5,0),
                     name=NULL,
                     labels=c('Eucalyptus','Native vegetation',
                              'Pasture'))+
    geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_point(size=3, position=position_dodge(.8))+
  geom_pointrange(aes(x=depth,y=chglnmn,colour=LUlongish,ymax=chglnmax,
                      ymin=chglnmin),fatten=8,
                  data=distinct(mr2sub,LUlongish,element,depth,.keep_all=T),
                  shape=16,show.legend = F,position=position_dodge2(.8))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  scale_y_continuous(sec.axis = 
                       sec_axis(trans=~.,name='Percent change in stock',
                                breaks=pct_to_L(c(-75,-25,25,75,200, 1000)),
                                labels=paste(c('-75', '-25','+25','+75',
                                               '+200', '+1000'),'%',sep='')))+
  labs(y='ln(Stock in 2016 / Stock in 2004)', x="Depth (cm)") +
  theme(strip.text.x = element_text(angle = 0,size=16),
        strip.background = element_rect(fill='grey80'),#,size=1.3),
#        legend.position=c(0.6,0.8),
        legend.position=c(0.85,0.8),
        panel.spacing = unit(2, "lines"),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank(),
        legend.title=element_text(size=16),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=12))+ 
  geom_text(mapping = aes( y=maxchgln*1.15*!is.na(sigveg),#x=LUlongish,
                           x=depth,#colour=LUlongish,
                           label = sigveg),
            size=8,na.rm=T,show.legend=F,colour='black',
            position=position_dodge(.8))+
  geom_text(mapping = aes(y=(maxchgln)*1.15*!is.na(sigyr),#x=LUlongish, 
                          x=depth,#colour=LUlongish,
                          label = sigyr),
            size=7,na.rm=T,show.legend=F, colour='black',
            position=position_dodge(.8))
dev.off()

# new big simple
png('fig3_solo_CCa.png',res=150,height=6,width=9,units='in')
ggplot(data=distinct(mr2sub,LUlongish,element,depth,.keep_all=T),
       aes(x=depth, y=chgln,colour=LUlongish,shape=LUlongish))+
  scale_x_discrete()+
  scale_colour_manual(values=c('#000000','#0072B2','#e79f00'),
                      name=NULL,
                      labels=c('Eucalipto','Nativa',
                               'Pastagem'))+
  scale_shape_manual(values=c(18,17,15),
                     name=NULL,
                     labels=c('Eucalipto','Nativa',
                              'Pastagem'))+
  geom_hline(yintercept=0,show.legend = F,colour='grey60') +
  geom_pointrange(aes(x=depth,y=chglnmn,colour=LUlongish,ymax=chglnmn+chglnse,
                      ymin=chglnmn-chglnse,shape=LUlongish),fatten=8,#shape=18,
                  show.legend = T,position=position_dodge2(.8))+
  facet_wrap(~element,ncol=3,scales='free_y') +
  scale_y_continuous(sec.axis = 
                       sec_axis(trans=~.,name='Δ Estoque (%)',
                                breaks=pct_to_L(c(-75,-25,25,75,200, 1000)),
                                labels=paste(c('-75', '-25','+25','+75',
                                               '+200', '+1000'),'%',sep='')))+
  labs(y='ln(Estoque 2016 / Estoque 2004)', x="Profundidade (cm)") +
  guides(colour=guide_legend())+
  theme(strip.text.x = element_text(angle = 0,size=16),
        strip.background = element_rect(fill='grey80'),#,size=1.3),
        #        legend.position=c(0.6,0.8),
        #legend.position=c(0.9,0.2),
        panel.spacing = unit(2, "lines"),
        legend.spacing.y = unit(.5,'lines'), 
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank(),
        legend.key=element_blank(),
        legend.title=element_text(size=16),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16))+ 
  geom_text(mapping = aes( y=maxse*1.1*!is.na(sigveg),#x=LUlongish,
                           x=depth,colour=LUlongish,
                           label = sigveg),
            size=8,na.rm=T,show.legend=F,#colour='black',
            position=position_dodge2(.8))+
  geom_text(mapping = aes(y=(maxse)*1.1*!is.na(sigyr),#x=LUlongish, 
                          x=depth,colour=LUlongish,
                          label = sigyr),
            size=7,na.rm=T,show.legend=F, #colour='black',
            position=position_dodge2(.8))
dev.off()


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

# Figure S1
######## Row position and other heterogeneity things
##################
eltdats=dats[dats$elt %in% c('E','L','T') & dats$LU=='E'&dats$year=='16',]
eltdats=mutate(eltdats,eltlong=ifelse(elt=='E','Inter-',
                                      ifelse(elt=='L','Current','Previous')))
eltdats$eltlong=factor(eltdats$eltlong,levels=c('Current','Previous',
                                                'Inter-'))
eltdats=group_by(eltdats,stand,depth,element,rep,elt)%>%
  mutate(eltmn=mean(value,na.rm=T))

eltstands=c(Bp.E1 = "Bom Despacho 1", Eu.E2 = "Eunápolis 2",
            It.E1 = "Itacambira 1", JP.E2 = "João Pinheiro 2")

png('figs2_was_s1.png',res=100,height=5,width=6,units='in')

ggplot(data=eltdats[eltdats$element=='C'&
                      eltdats$stand %in% c('Bp.E1','Eu.E2','It.E1','JP.E2')&
                      eltdats$depth==5,],
       aes(x=eltlong,y=eltmn, color=as.factor(rep)))+
  geom_point(shape=18,size=3,show.legend=F)+
  geom_line(aes(group=as.factor(rep)),show.legend = F)+
  facet_wrap(~stand,ncol=2,scales='free_y',
             labeller = labeller(stand=eltstands))+
  theme(panel.background = element_rect(fill='white'),
        panel.grid.major = element_blank())+
#  labs(y='Soil carbon content (g / 100 g), 0-10 cm',x='Row position')+
  labs(y='Soil carbon (%), 0-10 cm',x='Row position')+
  scale_color_brewer(palette='Dark2')
dev.off()

# Table S1
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

# Table S2
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
# CV within depth, stand, and year
depcvs=group_by(droplevels(dats4[dats4$site!='TM'&dats4$site!='Cr'&
                                   dats4$LU!='A',]),
                stand,LU,biome,depth,element,year) %>%
  summarise(CV=sd(repval,na.rm=T)/mean(repval,na.rm=T))
# Euc
tapply(droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K')&
                           depcvs$LU=='E',])$CV,
       droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K')&
                           depcvs$LU=='E',])$element,summary)
# Native
tapply(droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K')&
                           depcvs$LU=='N',])$CV,
       droplevels(depcvs[depcvs$element %in% c('C','N','P2','Ca2','K')&
                           depcvs$LU=='N',])$element,summary)

# CV of stock
stockcvs=group_by(dats2deps,stand,LU,biome,element,year) %>%
  summarise(CV20=sd(stock20,na.rm=T)/mean(stock20,na.rm=T),
            CV100=sd(stock100,na.rm=T)/mean(stock100,na.rm=T))
tapply(droplevels(stockcvs[stockcvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn')&
                             stockcvs$LU=='E',])$CV20,
       droplevels(stockcvs[stockcvs$element %in% 
                             c('C','N','P2','Ca2','K','S','Zn')&
                             stockcvs$LU=='E',])$element,summary)


