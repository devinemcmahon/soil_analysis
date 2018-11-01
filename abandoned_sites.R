source('soil_data_reader.R')

dats4a=dats4[is.element(dats4$site, c('TM','Cr','JP'))&
               !is.element(dats4$stand,c('JP.P'))&dats4$year==16,]

Pdepaov=aov(repval~depth*LU,data=dats4a[dats4a$element=='P',])
qqr(Pdepaov) # pretty ok!
summary(Pdepaov) # depth effect, no LU effect or interaction

Pdeplme=lme(repval~depth*LU,data=dats4a[dats4a$element=='P',],
            random=~depth|site,na.action=na.omit)
qqr(Pdeplme) # high outlier
summary(Pdeplme) # less P in A and N than E, maybe less superficial in N,
# definite decrease with depth
# driven mostly by single JP.E stand? but the other JP.E counteracts it

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
        par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                  lwd = 2)),
        auto.key=list(space='top', columns=3,lines=T, points=F))
# why does Curvelo have so much K?

xyplot(depth~mn*10|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='N' ,],
       ylim=c(90,0),xlab='N (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
xyplot(depth~mn/1000|site,groups=stdonly,type='l',ylab='Depth (cm)',
       data=abdatsmn[abdatsmn$element=='Zr' ,],
       ylim=c(90,0),xlab='Zr (g / kg)',as.table=T,layout=c(3,1),
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))
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
       par.settings = list(superpose.line = list(col = c(2,2,2,4,4,4,3),
                                                 lwd = 2)),
       auto.key=list(space='top', columns=3,lines=T, points=F))

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

