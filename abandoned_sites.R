source('soil_data_reader.R')

dats4a=dats4[is.element(dats4$site, c('TM','Cr','JP'))&
               !is.element(dats4$stand,c('JP.P'))&dats4$year==16,]
dats4a$LU=factor(dats4a$LU,levels=c('E','N','A'))

library(broom)
library(car)
#anova4a=group_by(droplevels(dats4a[dats4a$element %in% c('C','N','P','K','Ca'),]),
#                depth,element,unit,site) %>%
#  summarise(pval=anova(aov(repval~LU,data = .))$`Pr(>F)`[1]) #no
head(anova4a)

dats4asimp=droplevels(dats4a[dats4a$element %in% c('C','N','P2','K','Ca2'),])
anova4a=group_by(dats4asimp,depth,element,site) %>%
  do(tidy(Anova(aov(repval~LU,data = .),type="III"))) %>% filter(term=="LU")
anova4asig=anova4a[anova4a$p.value<0.05,]
View(anova4asig)
#anova4a=group_by(dats4asimp,depth,element,site) %>%
#  group_map(~anova(aov(repval~LU,data=.x))) 
# works, but doesn't keep the term names

dats4asig=merge(dats4asimp,anova4a,by=c('depth','element','site'))


library(dplyr)
Pdepaov=aov(repval~depth*LU,data=dats4a[dats4a$element=='P2',])
qqr(Pdepaov) # pretty ok!
summary(Pdepaov) # depth effect, no LU effect or interaction
# with P2, almost an LU effect, fewer missing values

Cdepaov_JP=aov(repval~depth*LU,data=dats4a[dats4a$element=='C' &
                                          dats4a$site=='JP',])
qqr(Cdepaov_JP) 
summary(Cdepaov_JP) 
Anova(Cdepaov_JP)
anova(Cdepaov_JP)
C5aov_JP=aov(repval~LU,data=dats4a[dats4a$element=='C' &dats4a$depth==5&
                                             dats4a$site=='JP',])
Anova(C5aov_JP,type='III')
anova(C5aov_JP) # same thing except Anova also includes info on the intercept
qqr(C5aov_JP)
TukeyHSD(C5aov_JP)

Pdepaov_TM=aov(log(repval)~depth*LU,data=dats4a[dats4a$element=='P' &
                                             dats4a$site=='TM',])
qqr(Pdepaov_TM) # pretty ok; log maybe not needed
TukeyHSD(Pdepaov_TM,which='LU') 
# all different; native vs aband most different
Cdepaov_TM=aov(log(repval)~depth*LU,data=dats4a[dats4a$element=='C' &
                                             dats4a$site=='TM',])
qqr(Cdepaov_TM) # upper tail is off; log is good
TukeyHSD(Cdepaov_TM,which='LU') # A has least, N=E

Pdeplme=lme(repval~depth*LU+I(depth^2),data=dats4a[dats4a$element=='P',],
            random=~depth|site,na.action=na.omit)
qqr(Pdeplme) # high outlier
summary(Pdeplme) # less P in A and N than E, maybe less superficial in N,
# definite decrease with depth
# driven mostly by single JP.E stand? but the other JP.E counteracts it
# adding a squared term: lots of significant things, except no depth*A intrxn

Cdeplme=lme(log(repval)~depth*LU,data=dats4a[dats4a$element=='C',],
            random=~depth|site,na.action=na.omit)
qqr(Cdeplme) # couple of possible high outliers but pretty ok with log
summary(Cdeplme)
# also decrease with depth, noneuc; lesser decrease with depth in N
 
xyplot(depth~mn|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='C' ,],
       ylim=c(90,0),xlab='C (g / 100 g)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# looks like euc has the most 0-20 cm C in all 3 sites

xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
        data=abdatsmn[abdatsmn$element=='K' ,],
        ylim=c(90,0),xlab='K (g / kg)',as.table=T,layout=c(3,1),
 #       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
  #                                                lwd = 2)),
        auto.key=list(space='top', columns=3,lines=T, points=F))
# why does Curvelo have so much K?
# wrt K, A1 matches E, A2 matches N (more K) there; in TM, A has less K
# slope position? At Cr, A1 and E are adjacent; A2 above, N below
# in TM, E above N, abuts A-- more K on oe side of road?
# slope pretty minimal in both places
# K distribution in JP.N weird, too 
plot(repval~depth,data=dats4a[dats4a$site=='JP'&dats4a$element=='K',],
     col=as.numeric(stand)+1,pch=18,ylab='K (mg kg-1)')
legend('topleft',col=as.numeric(unique(dats4a$stand[dats4a$site=='JP']))+1,
       legend=unique(dats4a$stand[dats4a$site=='JP']),pch=15)

xyplot(depth~mn*10|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='N' ,],
       ylim=c(90,0),xlab='N (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# Elevated N at 0-10 and 60-100 cm in TM.N and (to a lesser extent) Cr.N
#   affected by artifacts of N drift correction (fixed now?)
# re-ran JP.N 60-100 but N concs still higher there than 40-60 and other veg
xyplot(depth~repval*10|site,groups=stdonly,type='p',ylab='Depth (cm)',
       data=dats4a[dats4a$element=='N' ,],
       ylim=c(90,0),xlab='N (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.symbol = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)))

xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='Zr' ,],
       ylim=c(90,0),xlab='Zr (g / kg)',as.table=T,layout=c(3,1),
       #par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
      #                                           lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# Pretty variable; same pattern as K and Cu in TM and Cr,
#   but also variable in JP; less in active euc there
# More Zr = more weathered, but also more K
xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='P' ,],
       ylim=c(90,0),xlab='P (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='Ca' ,],
       ylim=c(90,0),xlab='Ca (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# added in large quantities in E in JP and TM, not in Cr or A
xyplot(depth~mn|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='Al' ,],
       ylim=c(90,0),xlab='Al (g / 100 g)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# acceptably similar I guess



xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='Cu' ,],
       ylim=c(90,0),xlab='Cu (g / kg)',as.table=T,layout=c(3,1),
      # par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
      #                                           lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# Cu follows same patterns as K
TMCrw=droplevels(widedats4[is.element(widedats4$site, c('TM','Cr'))&
                  widedats4$year==16,])
plot(Cu~K,data=TMCrw,col=stand,pch=18)#,
     #pch=as.numeric(as.factor(depth)))
legend('bottomright',col=as.numeric(unique(TMCrw$stand)),
       legend=unique(TMCrw$stand),pch=15,bty='n')
# K varies by rep, Cu more by depth within a stand
# Native veg: higher Cu:K ratio?
plot(K~Zr,data=TMCrw,col=stand,pch=18) 
# interesting? - relationship in Cr, blobbily + in TM
plot(K~Zr,data=TMCrw,col=depth,pch=18)# cluster by rep not depth


plot(C~P,data=TMCrw,col=stand,pch=18)
# nicely correlated, with diff slope in each site, 
#   maybe diff intercepts per stand
plot(N~C,data=TMCrw,col=stand,pch=18) # pretty tight
plot(N~C,data=TMCrw,col=depth,pch=18)
legend('bottomright',col=as.numeric(unique(TMCrw$depth)),
       legend=unique(TMCrw$depth),pch=15,bty='n')
# TM.N 60-100 cm and 0-10 cm have more N in reps 1 and 2
# So does Cr.N, weird
# ohh that's bad. Because those are the ones analyzed separately
# earlier than all the other samples
# Was N lost in the oven?
# No, but there are artifacts of serious N drift and correction that day
palette(rainbow(20))
plot(N~as.numeric(CNevalday),data=widedats,col=stand)
palette('default')
plot(value~as.numeric(CNevalday),col=stand,xaxt='n',
     data=droplevels(dats[dats$site %in% c('TM','Cr') & 
                            dats$element=='N',]))
axis(side=1,at=seq(length(unique(droplevels(
  dats[dats$site %in% c('TM','Cr') & dats$element=='N',])$CNevalday))),
  labels = unique(droplevels(
    dats[dats$site %in% c('TM','Cr') & dats$element=='N',])$CNevalday))

# Density
xyplot(depth~avgBD|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='C' ,],
       ylim=c(90,0),xlab='Bulk density (g cm-3)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
# Cr: generally denser in N, except at surface, 
#     where euc soil is densest and N and As least dense
# TM: little variation within or among stands

abdatsmn=droplevels(datsmean[is.element(datsmean$site, c('TM','Cr','JP'))&
                               !is.element(datsmean$stand,c('JP.P'))&
                               datsmean$year==16,])

abdatsmn$mn[abdatsmn$element %in% c('C','N')]=
  abdatsmn$mn[abdatsmn$element %in% c('C','N')]*10
abdatsmn$sd[abdatsmn$element %in% c('C','N')]=
  abdatsmn$sd[abdatsmn$element %in% c('C','N')]*10
abdatsmn$mn[!is.element(abdatsmn$element, c('C','N'))]=
  abdatsmn$mn[!is.element(abdatsmn$element, c('C','N'))]/1000
abdatsmn$sd[!is.element(abdatsmn$element, c('C','N'))]=
  abdatsmn$sd[!is.element(abdatsmn$element, c('C','N'))]/1000
abdatsmn$LU=factor(abdatsmn$LU,levels=c('E','N','A'))
# All in units of g/kg now

abdatsmnsig=merge(abdatsmn,anova4a,by=c('depth','element','site'))
abdatsmnsig=group_by(abdatsmnsig,element,site)%>%
  mutate(maxval=max(mn,na.rm=T))
abdatsmnsig$element[abdatsmnsig$element=='Ca2']='Ca'
abdatsmnsig$element[abdatsmnsig$element=='P2']='P'
abdatsmnsig$element=factor(abdatsmnsig$element,
                           levels=c('C','N','P','K','Ca'))

ggplot(aes(x=depth,y=mn,colour=LU),
       data=abdatsmnsig[abdatsmnsig$site=='JP',])+
       #data=abdatsmn[abdatsmn$element %in% c('C','N','K','P','Ca')&#,])+
        #               abdatsmn$site=='JP',])+
  geom_line(aes(group=stand),size=1.2)+
  coord_flip()+
  scale_x_reverse(breaks=c(100,60,40,20,10,0),minor_breaks=NULL)+
#  facet_grid(site~element,scales = 'free_x')+
  facet_wrap(~element,scales = 'free_x')+
  theme(legend.position=c(0.85,0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_rect(fill='white',colour='grey80'),
        panel.grid.major.x = element_line(colour='grey80'),
        panel.grid.major.y = element_line(colour='grey80'))+
  labs(y='Concentração (g elemento / kg solo)',
       x='Profundidade (cm)')+
  geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                width=0.2)+ 
  scale_colour_discrete(labels=c('Manejo ativo','Reserva nativa',
                                 'Abandonado'),
                        guide = guide_legend(reverse=F,title=NULL))+
  geom_text(aes(y=maxval*1.1,
    label=ifelse(p.value<0.05 &LU=='N','*','')),
    colour='black',size=6,nudge_x=-1)


# English version
ggplot(aes(x=depth,y=mn,colour=LU),data=abdatsmnsig)+
  geom_line(aes(group=stand),size=1.2)+
  coord_flip()+
  scale_x_reverse(breaks=c(100,60,40,20,10,0),minor_breaks=NULL)+
  facet_grid(site~element,scales = 'free_x')+
  #facet_wrap(~element,scales = 'free_x')+
  theme(legend.position=c(0.72,0.48),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_rect(fill='white',colour='black'),
        panel.grid.major.x = element_line(colour='grey80'),
        panel.grid.major.y = element_line(colour='grey80'))+
  labs(y='Concentration (g/kg)',
       x='Depth (cm)')+
  geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                width=0.2)+ 
  scale_colour_discrete(labels=c('Active stand','Native reserve',
                                 'Abandoned stand'),
                        guide = guide_legend(reverse=F,title=NULL))+
  geom_text(aes(y=maxval*1.1,
                label=ifelse(p.value<0.05 &LU=='N','*','')),
            colour='black',size=6,nudge_x=-1)
#ggsave('abdeps_all.png',device=png(height=8,width=10,units='in',res=150))

ab2deps=group_by(datsstk[is.element(datsstk$site,c('TM','Cr','JP'))&
                              datsstk$stand!='JP.P'&
                           datsstk$year=='16',],stand,rep,
                   element,site,LU,biome,unit,stockunit) %>% 
  summarise(ndeps20=length(unique(inc_to[inc_to<=20])),
            stock20=ifelse(ndeps20==2,sum(stock[inc_to<=20]),NA),
            stocklo20=ifelse(ndeps20==2,sum(stocklo[inc_to<=20]),NA),
            stockhi20=ifelse(ndeps20==2,sum(stockhi[inc_to<=20]),NA),
            oldBDstock20=ifelse(ndeps20==2,sum(oldBDstock[inc_to<=20]),NA),
            conc20=ifelse(ndeps20==2,
                          sum(repval[inc_to<=20]*inc[inc_to<=20]*avgBD[inc_to<=20])/
                            sum(inc[inc_to<=20]*avgBD[inc_to<=20]),NA),
            # weight concentrations by mass of soil in each layer
            # to get average density of the whole 1-m block of soil
            BD20=ifelse(ndeps20==2,
                        sum(avgBD[inc_to<=20]*inc[inc_to<=20])/
                          sum(inc[inc_to<=20]),NA),
            #BDsd20=ifelse(ndeps20==2,
            #           sum(BDsd[inc_to<=20]*inc[inc_to<=20])/
            #             sum(inc[inc_to<=20]),NA),
            BDsd20=ifelse(ndeps20==2,sqrt(sum(BDsd[inc_to<=20])^2),NA),
            stocksd20=ifelse(ndeps20==2,sqrt(sum(stockvar[inc_to<=20])),NA),
            # Variance of sum = sum of variances; check how best to do this
            ndeps100=length(unique(inc_to)),
            stock100=ifelse(ndeps100==5,sum(stock),NA),
            #ifelse(maxrep==5 & year=='16',
            #       sum(stock[inc_to<=40])+sum(stock),NA)),
            
            stocklo100=ifelse(ndeps100==5,sum(stocklo),NA),
            stockhi100=ifelse(ndeps100==5,sum(stockhi),NA),
            oldBDstock100=ifelse(ndeps100==5,sum(oldBDstock),NA),
            #conc100=ifelse(ndeps100==5,
            #               sum(repval*inc*avgBD)/sum(inc*avgBD),NA),
            conc100=ifelse(ndeps100==5,
                           sum(repval*inc*avgBD)/sum(inc*avgBD),NA),
            conc60to100=mean(repval[inc_to==100],na.rm=T),
            BD100=ifelse(ndeps100==5,sum(avgBD*inc)/sum(inc),NA),
            #BDsd100=ifelse(ndeps100==5,sqrt(sum(BDsd^2)),NA),
            # weight by inc
            suminc=sum(inc),
            BDsd100=ifelse(ndeps100==5,sqrt(sum((BDsd^2)*inc)/sum(inc)),NA),
            stocksd100=ifelse(ndeps100==5,sqrt(sum(stockvar*inc)/sum(inc)),NA),
            stockratio=stock20/stock100,concratio=conc20/conc100,
            concrat2=conc20/conc60to100)

abN20.lme=lme(stock20~LU,random=~1|site,
               data=ab2deps[ab2deps$element=='N',],na.action = na.omit)
qqr(abN20.lme) # not very good, with or without log
summary(abN20.lme) # native different (more N)

astksum =group_by(ab2deps[ab2deps$element %in% 
                            c('C','N','K','Ca2','P2','Al','Fe'),],
                  stand,element,site,LU) %>%
  summarise(stk20=mean(stock20,na.rm=T),
            nstocks20=sum(!is.na(stock20)),
            se20=sd(stock20,na.rm=T)/sqrt(nstocks20-1),
            stk100=mean(stock100,na.rm=T),
            nstocks100=sum(!is.na(stock100)),
            se100=sd(stock100,na.rm=T)/sqrt(nstocks100-1),
            BD20=mean(BD20),BD100=mean(BD100))
astksum$element[astksum$element=='Ca2']='Ca'
astksum$element[astksum$element=='P2']='P'
astksumr=mutate_if(astksum,is.numeric,round,digits=2)
astksumr=astksumr[,-which(names(astksumr) %in% 
                            c('nstocks100','nstocks20'))]
#write.csv(astksumr,'abd_stocksummaries.csv')

ab2deps$LU=factor(ab2deps$LU,levels=c('E','A','N'))
astksum2 =group_by(ab2deps[ab2deps$element %in% 
                            c('C','N','K','Ca2','P2','Al','Fe'),],
                  element,site) %>%
  summarise(stk20E=mean(stock20[LU=='E'],na.rm=T),
            stk20A=mean(stock20[LU=='A'],na.rm=T),
            stk20N=mean(stock20[LU=='N'],na.rm=T),
            nstocks20E=sum(!is.na(stock20[LU=='E'])),
            nstocks20A=sum(!is.na(stock20[LU=='A'])),
            nstocks20N=sum(!is.na(stock20[LU=='N'])),
            se20=sd(stock20,na.rm=T)/sqrt(nstocks20-1),
            stk100=mean(stock100,na.rm=T),
            nstocks100=sum(!is.na(stock100)),
            se100=sd(stock100,na.rm=T)/sqrt(nstocks100-1),
            BD20=mean(BD20),BD100=mean(BD100))
astksum$element[astksum$element=='Ca2']='Ca'
astksum$element[astksum$element=='P2']='P'
astksumr=mutate_if(astksum,is.numeric,round,digits=2)
astksumr=astksumr[,-which(names(astksumr) %in% 
                            c('nstocks100','nstocks20'))]
#write.csv(astksumr,'abd_stocksummaries.csv')


C20aov_TM=aov(stock20~LU,data=ab2deps[ab2deps$site=='TM' &
                                        ab2deps$element=='C',])
qqr(C20aov_TM) # not much different with log transform
anova(C20aov_TM) # no difference
anova20a=group_by(ab2deps[ab2deps$element %in% 
                            c('C','N','K','Ca2','P2','Al','Fe'),],
                  element,site) %>%
  do(tidy(Anova(aov(stock20~LU,data = .),type="III"))) %>% filter(term=="LU")
head(anova20a)
anova20asig=anova20a[anova20a$p.value<0.05,]
View(anova20asig)
# Ca and P vary in all 3 sites (ok)
# in TM, Fe and N also
N20aov_TM=aov(stock20~LU,data=ab2deps[ab2deps$site=='TM' &
                                        ab2deps$element=='N',])
qqr(N20aov_TM) 
TukeyHSD(N20aov_TM) # N marginally more than E and a lot more than A
Fe20aov_TM=aov(stock20~LU,data=ab2deps[ab2deps$site=='TM' &
                                        ab2deps$element=='Fe',])
qqr(Fe20aov_TM) #ugh (log worse)
TukeyHSD(Fe20aov_TM) # A has the least Fe...not enough data to tell?
