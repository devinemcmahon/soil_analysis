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
names(pts)[which(names(pts)=='element')]='element2'

datsmnok2=mutate(datsmnok,element2=element)
datsmnok2$element2[datsmnok2$element=='P2']='P'
datsmnok2=merge(datsmnok2,pts,by=c('stand','element2','depth'),all.x=T)
datsmnok2$element2[datsmnok2$element=='Ca2']='Ca'
datsmnok2=datsmnok2[order(datsmnok2$depth),]
# rearranged this part to attach significance from P to P2
standlong=datsmnok2$stand
standlong=gsub("Vg.E","Talhão 8 Cataquinho",standlong)
standlong=gsub("Vg.N","Reserva Cataquinho",standlong)
standlong=gsub("BO.E","Talhão 30 Lagoa Cristal",standlong)
standlong=gsub("BO.P","Pastagem Lagoa Cristal",standlong)
standlong=gsub("Eu.E1","Talhão 23 Sucupira",standlong)
standlong=gsub("Eu.E2","Talhão 19 Inhaíba",standlong)
standlong=gsub("Eu.N","Reserva Inhaíba",standlong)
standlong=gsub("Bp.E2","Talhão 683 Garça",standlong)
standlong=gsub("Bp.E1","Talhão 567 Extrema",standlong)
standlong=gsub("It.E1","Talhão 401 Itacambira",standlong)
standlong=gsub("It.E2","Talhão 123 Itacambira",standlong)
standlong=gsub("It.N","Reserva Pau Preto",standlong)
#standlong=gsub("JP.E1","Talhão 1423 Campo Alegre",standlong)
standlong=gsub("JP.E1","Talhão 1423",standlong)
standlong=gsub("JP.E2","Talhão 1671",standlong)
standlong=gsub("JP.P","Pastagem Faz. Pontes",standlong)
standlong=gsub("JP.N","Reserva Campo Alegre",standlong)
datsmnok2$standlong=standlong

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



mygg=function(rockder,mysite){
  lmts=c('C','N')
  if(rockder==T) lmts=c('K','P2','Ca2')
  mysub=datsmnok2[datsmnok2$site==mysite & datsmnok2$element %in% lmts,]
  mysub=mutate(mysub,element2=element)
  mysub$element2[mysub$element=='Ca2']='Ca'
  mysub$element2[mysub$element=='P2']='P'
  if(rockder==T) mysub$element2=factor(mysub$element2,
                                       levels=c('K','Ca','P'))
  mysub=mysub[order(mysub$depth),]
  ggplot(aes(x=depth,y=mn,colour=year),data=mysub)+
    geom_line(aes(group=year),size=1.2)+
    coord_flip()+
    scale_x_reverse(breaks=c(100,60,40,20,10,0),minor_breaks=NULL)+
    facet_grid(standlong~element2,scales = 'free_x')+
    theme(legend.position=c(0.9,0.1),
          legend.background = element_blank(),
          legend.key = element_blank())+
    labs(y=ifelse(rockder==T,'Concentração (mg elemento / kg solo)',
                  'Concentração (g elemento / 100g solo)'),
         x='Profundidade (cm)')+
    geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                  width=0.2) +
    scale_colour_discrete(labels=c('2004','2016'),
                          guide = guide_legend(reverse=F,title=NULL))+
    geom_text(aes(#y=mn+(I(sd/sqrt(ndepyr)))*1.1,
      label=ifelse(pval<0.05 &year=='16','*','')),
      colour='black',size=6,nudge_x=-1)
}
mygg(F,'BO')

ggplot(aes(x=depth,y=mn,colour=year),data=datsmnok2[datsmnok2$site=='Vg'&
                                             datsmnok2$element %in%
                                             c('K','Ca2','P2'),])+
  geom_line(aes(group=year),size=1.2)+
  coord_flip()+
  scale_x_reverse(limits = c(90, 0))+
  facet_grid(stand~element2,scales = 'free_x')+
  theme(legend.position=c(0.9,0.1),
        legend.background = element_blank(),
        legend.key = element_blank())+
  labs(y=ifelse(rockder==T,'Concentração (mg elemento / kg solo)',
                   'Concentração (g elemento / 100g solo)'),
       x='Profundidade (cm)')+
  geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                width=0.2) +
  scale_colour_discrete(labels=c('2004','2016'),
                    guide = guide_legend(reverse=F,title=NULL))+
  geom_text(aes(#y=mn+(I(sd/sqrt(ndepyr)))*1.1,
                label=ifelse(pval<0.05 &year=='16','*','')),
                colour='black',size=6,nudge_x=-1)

ggplot(aes(x=depth,y=mn,colour=year),data=datsmnok2[datsmnok2$element =='CN',])+
  geom_line(aes(group=year),size=1.2)+
  coord_flip()+
  scale_x_reverse(limits = c(90, 0))+
  facet_wrap(stand~.)+
  theme(legend.position=c(0.9,0.3),
        legend.background = element_blank(),
        legend.key = element_blank())+
  labs(y='C:N',
       x='Profundidade (cm)')+
  geom_errorbar(aes(ymax=mn+I(sd/sqrt(ndepyr)),  ymin=mn-I(sd/sqrt(ndepyr))),
                width=0.2) +
  scale_colour_discrete(labels=c('2004','2016'),
                        guide = guide_legend(reverse=F,title=NULL))+
  geom_text(aes(#y=mn+(I(sd/sqrt(ndepyr)))*1.1,
    label=ifelse(pval<0.05 &year=='16','*','')),
    colour='black',size=6,nudge_x=-1)


  labls <- c(AF = "Atlantic Forest", Cer = "Cerrado")

ggplot(data=nutstkE, aes(x=year, y=stock, fill=depth)) +
  geom_bar(stat="identity") + 
  facet_grid(element~biome,labeller = labeller(biome=labls))  
  

shortsum=shorttstk[shorttstk$element%in%
                     c('C','N','K','Ca','P','Al','Fe','Zr'),
          c('stand','LU','biome','site','element',
            'stock20_04','stock20_16','sd20_04','sd20_16',
            'stock100_04','stock100_16','sd100_04','sd100_16',
            'BD20_16','BD100_16','tstat20','pval20','tstat100','pval100')]
View(shortsum)
shortsumr=round(shortsum,2)
#write.csv(shortsum,'allcos_stocksummaries.csv')


stkchgs3=stkchgs[stkchgs$element!='Mg',]
stkchgs3$element=factor(stkchgs3$element,levels=c('N','P','K','Ca'))
cochgtypes=group_by(stkchgs3,stand,element) %>%
  summarise(fertilizante=(In_kgha_1+In_kgha_2),
            retirado=(Wood_m3_1+Wood_m3_2)*Concentration*511,
            balanco=budget*1000,minbal=minbudg*1000,maxbal=maxbudg*1000,
            mudanca20=chg20*1000,mudEP=sdchg20*1000,
            mudaerea=agbchg*-1000,balaerea=agbbudg*1000,
            minaerbal=minagbbudg*1000,maxaerbal=maxagbbudg*1000)%>%
  mutate_if(is.numeric,round,digits=0)
View(cochgtypes[cochgtypes$stand=='BO.E',])
write.csv(cochgtypes,'cochgsum.csv')

names(stkchgs3)
stkchgs$element[stkchgs$stand=='BO.E']
stkchgs$sdchg20[stkchgs$stand=='BO.E']
stkchgs$maxbudg[stkchgs$stand=='BO.E']
stkchgs$minbudg[stkchgs$stand=='BO.E']

