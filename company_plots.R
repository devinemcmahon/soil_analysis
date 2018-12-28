# Plots of changes in individual stands
source('soil_data_reader.R')

# Plots by year and depth in each stand

myplot=function(dats,stand,lmt){
  mysub=dats[dats$element==lmt & dats$stand==stand,]
  mysub=mysub[order(mysub$depth),]
  minx=min(c(dats$mn04[dats$element==lmt&dats$site %in% mysub$site],
             dats$mn16[dats$element==lmt&dats$site %in% mysub$site]),
           na.rm=T)
  maxx=max(c(dats$mn04[dats$element==lmt&dats$site %in% mysub$site],
             dats$mn16[dats$element==lmt&dats$site %in% mysub$site]),
           na.rm=T)
  plot(depth~mn04,type='l',data=mysub,
       xlim=c(minx*0.8,maxx*1.2),las=1,
       ylim=c(90,0),lwd=2,col=2,xlab=paste(lmt,unique(mysub$unit)),
       ylab='Profundidade (cm)'
       )
  lines(mysub$mn16,mysub$depth,col=4,lwd=2)
  segments(x0=mysub$mn04-I(mysub$sd04/sqrt(mysub$n04)),y0=mysub$depth-1,
    x1=mysub$mn04+I(mysub$sd04/sqrt(mysub$n04)),col=2)
  segments(x0=mysub$mn16-I(mysub$sd16/sqrt(mysub$n16)),y0=mysub$depth+1,
           x1=mysub$mn16+I(mysub$sd16/sqrt(mysub$n16)),col=4)
  #points(I(mysub$mn16+mysub$sd16),mysub$depth,pch=4,cex=.5,col=4)
  #points(I(mysub$mn16-mysub$sd16),mysub$depth,pch=4,cex=.5,col=4)
  text(rep(maxx*1.1),mysub$depth,label='*',cex=1.1,col=mysub$pval<0.05)
}
myplot(ttests,'BO.E','C')

pts=ttests[,c('stand','element','depth','pval','tstat')]
datsmnok2=merge(datsmnok,pts,by=c('stand','element','depth'),all.x=T)
myxy=function(rockder,mysite){
  lmts=c('C','N')
  if(rockder==T) lmts=c('K','P2','Ca2')
  mysub=datsmnok2[datsmnok2$site==mysite & datsmnok2$element %in% lmts,]
  mysub=mutate(mysub,element2=element)
  mysub$element2[mysub$element=='Ca2']='Ca'
  mysub$element2[mysub$element=='P2']='P'
  if(rockder==T) mysub$element2=factor(mysub$element2,
                                             levels=c('K','Ca','P'))
  mysub=mysub[order(mysub$depth),]
  #mypanel=panel.xyplot()
  if(rockder==T) xlims=list(c(min(mysub$mn[mysub$element2=='K']),
                              max(mysub$mn[mysub$element2=='K'])),
                            c(min(mysub$mn[mysub$element2=='Ca']),
                              max(mysub$mn[mysub$element2=='Ca'])),
                            c(min(mysub$mn[mysub$element2=='P']),
                              max(mysub$mn[mysub$element2=='P'])))
  if(rockder==F) xlims=list(c(min(mysub$mn[mysub$element2=='C']),
                             max(mysub$mn[mysub$element2=='C'])),
                           c(min(mysub$mn[mysub$element2=='N']),
                             max(mysub$mn[mysub$element2=='N'])))
  #ifelse(rockder==T,c('K','P','Ca'),c('C','N'))# lmts gets set to first element
  #ifelse(rockder==T,lmts=list('K','P','Ca'),lmts=list('C','N'))# doesn't work
  xyplot(depth~mn|element2+stand,groups=year,type='l',
         data=mysub,ylab='Profundidade (cm)',
         ylim=c(90,0),as.table=T,las=2,
         scales=list(relation="free",rot=c(0,90)),
         xlab=ifelse(rockder==T,'Concentração (mg elemento / kg solo)',
                     'Concentração (g elemento / 100g solo)'),
         xlim=xlims,
       par.settings = list(superpose.line = list(col = c(2,4),lwd = 2),
                           strip.background = list(col = c('grey80','steelblue2')),
                           par.strip.text=list(cex=.8)),
         auto.key=list(space='top', columns=2,lines=TRUE, points=FALSE))
}
myxy(F,'Vg')
myxy(T,'JP')

shortsum=shorttstk[shorttstk$element%in%
                     c('C','N','K','Ca','P','Al','Fe','Zr'),
          c('stand','LU','biome','site','element',
            'stock20_04','stock20_16','sd20_04','sd20_16',
            'stock100_04','stock100_16','sd100_04','sd100_16',
            'BD20_16','BD100_16','tstat20','pval20','tstat100','pval100')]
View(shortsum)
shortsumr=round(shortsum,2)
#write.csv(shortsum,'allcos_stocksummaries.csv')
