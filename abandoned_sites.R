source('soil_data_reader.R')

dats4a=dats4[is.element(dats4$site, c('TM','Cr','JP'))&
               !is.element(dats4$stand,c('JP.P'))&dats4$year==16,]

Pdepaov=aov(repval~depth*LU,data=dats4a[dats4a$element=='P',])
qqr(Pdepaov) # pretty ok!
summary(Pdepaov) # depth effect, no LU effect or interaction

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

ggplot(aes(x=depth,y=mn,colour=LU),
       data=abdatsmn[abdatsmn$element %in% c('C','N','K','P','Ca'),])+
  geom_line(aes(group=stand),size=1.2)+
  coord_flip()+
  scale_x_reverse(breaks=c(100,60,40,20,10,0),minor_breaks=NULL)+
  facet_grid(site~element,scales = 'free_x')+
  theme(legend.position=c(0.9,0.1),
        legend.background = element_blank(),
        legend.key = element_blank())+
  #labs(y=ifelse(rockder==T,'Concentração (mg elemento / kg solo)',
  #              'Concentração (g elemento / 100g solo)'),
  #     x='Profundidade (cm)')+
  geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                width=0.2) 
  scale_colour_discrete(labels=c('2004','2016'),
                        guide = guide_legend(reverse=F,title=NULL))+
  geom_text(aes(#y=mn+(I(sd/sqrt(ndepyr)))*1.1,
    label=ifelse(pval<0.05 &year=='16','*','')),
    colour='black',size=6,nudge_x=-1)