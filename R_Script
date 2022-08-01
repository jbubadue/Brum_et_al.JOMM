##############################################################################################
####################################Script for Brum et al.####################################
###Evolutionary rates, disparity, and ecomorphology of the mandible of American marsupials####
################################Journal of Mammalian Evolution################################
##############################################################################################

##Load packages
library(geomorph)#pacotes
library(ape)
library(ggplot2)
library(Morpho)
library(phytools)
library(ggtree)
library(ggplot2)
library(RRPP)
library(ggridges)
setwd("~/")####select directory


mandible<-readland.nts("mandible.nts") #opening nts file
datamandible<-read.table("datamandible.txt",h=T,sep=",")#opening additional data file
tree<-read.tree("AMtreeUpham.txt")#tree

#Plotting tree
plotTree(tree, fsize=0.9, lwd=2, ftype="i")

#sliders<-define.sliders(all,nsliders = 34)###manually defining sliders -> or use the ready to go file "sliders.txt"
sliders<-read.table("sliders.txt",h=F) ###ready to go sliders 
gpa<-gpagen(mandible,curves=sliders)### GPA with sliders
plot(gpa) #check overall variation cloud and the consensus shape


##trait gram -> Figure 3
size<-datamandible$Log_Centroid_Size###size
names(size)<-datamandible$species
plot<-phenogram(tree,size,spread.labels=TRUE,spread.cost=c(2,0))

####create geomorph.data.frame
gdf<-geomorph.data.frame(gpa, phy=tree, speciesID=datamandible$label, clade=datamandible$Clade,
                         families=datamandible$Family, lc=datamandible$Locomotor_Category,
                         diet1=datamandible$Diet_1,diet2=datamandible$Diet_2,
                         diet3=datamandible$Diet_3,diet4=datamandible$Diet_4,size=datamandible$LOGCS)

###plotting preparation
clades<-as.factor(datamandible$Clade)
### Make clade-level color vector
cols <- c("black", "yellow", "red", "pink", "dark green", "light green", 
          "#9966FF") 
# specific colors for clades; can be hex codes or real names
clades<-as.factor(datamandible$Clade)
names(cols)=levels(as.factor(datamandible$Clade))
col.gp <- cols[match(clades, names(cols))] 


####PCA
PCA.w.phylo<-gm.prcomp(gpa$coords,tree)


#PCA setup
# preparing axis labels. 
PCA_summary <- summary(PCA.w.phylo)

pc1lab <-paste("PC 1 ","(",
               round(100*PCA_summary$PC.summary[2,1],1),"%)",sep="")
pc2lab <-paste("PC 2 ","(",
               round(100*PCA_summary$PC.summary[2,2],1),"%)",sep="")
pc3lab <-paste("PC 3 ","(",
               round(100*PCA_summary$PC.summary[2,3],1),"%)",sep="")
pc4lab <-paste("PC 4 ","(",
               round(100*PCA_summary$PC.summary[2,4],1),"%)",sep="")

#linkmandible<-define.links(gpa$consensus, ptsize = 2, links = NULL) #define links manually, or use the ready to go file.
linkmandible<-read.table("linksmandible.txt",h=F)#ready to go link file

#####PCA plot

#Plot PC1 vs PC2
par(mar=c(4, 5, 0, 1))  # sets the margins c(bottom, left, top, right)

plot(PCA.w.phylo, phylo = TRUE, pch = 21, cex = 1,  cex.axis=1.3, cex.lab=1.3, col = col.gp, bg = col.gp, phylo.par = list(edge.color="black", node.cex=NULL,  tip.labels = F, node.labels = FALSE, edge.color = "gray10"),bty="n",asp=T) #if tip.labels = T, species names are plotted
clades <- levels(as.factor(clades))

plotRefToTarget(PCA.w.phylo$shapes$shapes.comp1$max,PCA.w.phylo$shapes$shapes.comp1$min,links=linkmandible,method="TPS") #PC1 min deformation
plotRefToTarget(PCA.w.phylo$shapes$shapes.comp1$min,PCA.w.phylo$shapes$shapes.comp1$max,links=linkmandible,method="TPS") #PC1 max deformation


plotRefToTarget(PCA.w.phylo$shapes$shapes.comp2$max,PCA.w.phylo$shapes$shapes.comp2$min,links=linkmandible,method="TPS") #PC2 min deformation
plotRefToTarget(PCA.w.phylo$shapes$shapes.comp2$min,PCA.w.phylo$shapes$shapes.comp2$max,links=linkmandible,method="TPS") #PC2 max deformation

#Plot PC3 vs PC4
plot(PCA.w.phylo, phylo = TRUE, axis1=3, axis2=4, pch = 21, cex = 1,  cex.axis=1.3, cex.lab=1.3, col = col.gp, bg = col.gp, phylo.par = list(edge.color="black", node.cex=NULL,  tip.labels = F, node.labels = FALSE, edge.color = "gray10"),bty="n",asp=T ) #if tip.labels = T, species names are plotted

plotRefToTarget(PCA.w.phylo$shapes$shapes.comp3$max,PCA.w.phylo$shapes$shapes.comp3$min,links=linkmandible,method="TPS")#PC3 min deformation
plotRefToTarget(PCA.w.phylo$shapes$shapes.comp3$min,PCA.w.phylo$shapes$shapes.comp3$max,links=linkmandible,method="TPS")#PC3 max deformation


plotRefToTarget(PCA.w.phylo$shapes$shapes.comp4$max,PCA.w.phylo$shapes$shapes.comp4$min,links=linkmandible,method="TPS")#PC4 min deformation
plotRefToTarget(PCA.w.phylo$shapes$shapes.comp4$min,PCA.w.phylo$shapes$shapes.comp4$max,links=linkmandible,method="TPS")#PC4 max deformation


#####phylogenetic signal
physignal(size,tree,iter=9999) #size
physignal(gpa$coords,tree,iter=9999) #shape

####pgls models for size
fit1<-procD.pgls(size~as.factor(diet1),tree,iter=99999,data=gdf)
summary(fit1)
fit2<-procD.pgls(size~as.factor(diet2),tree,iter=99999,data=gdf)
summary(fit2)
fit3<-procD.pgls(size~as.factor(diet3),tree,iter=99999,data=gdf)
summary(fit3)
fit4<-procD.pgls(size~as.factor(diet4),tree,iter=99999,data=gdf)
summary(fit4)
fit5<-procD.pgls(size~as.factor(lc),tree,iter=99999,data=gdf)
summary(fit5)

modcomp<-model.comparison(fit1,fit2,fit3,fit4,fit5,type="logLik")#model comparison
modcomp #best model diet4

####plot size density plots
ggplot(datamandible, aes(x = Log_Centroid_Size, y = Diet_4,fill=stat(x))) +
  geom_density_ridges_gradient(jittered_points = F) +
  scale_fill_viridis_c(name = "Size", option = "C") +
  coord_cartesian(clip = "off") + # To avoid cut off
  theme_minimal()+
  stat_summary(
    mapping = aes(x = Log_Centroid_Size, y=as.factor(Diet_4)),
    fun.min = min,
    fun.max = max,
    fun = mean
  )

####pgls models for shape

fit1<-procD.pgls(coords~size,phy,iter=9999,data=gdf)#allometry
summary(fit)

fit2<-procD.pgls(coords~diet1,phy,iter=9999,data=gdf)#diet1 no effect
summary(fit)

fit3<-procD.pgls(coords~diet2,phy,iter=9999,data=gdf)#diet2 no effect
summary(fit3)

fit4<-procD.pgls(coords~diet3,phy,iter=9999,data=gdf)#diet3 no effect
summary(fi4)

fit5<-procD.pgls(coords~diet4,phy,iter=9999,data=gdf)#diet4 no effect
summary(fit5)

fit6<-procD.pgls(coords~lc,phy,iter=9999,data=gdf)#lc no effect
summary(fit6)

modcomp<-model.comparison(fit1,fit2,fit3,fit4,fit5,fit6,type="logLik")#model comparison
modcomp #best model size

###plot allometry
fit<-procD.lm(coords~size,iter=9999,data=gdf)#alometry 
summary(fit)

alo<-plotAllometry(fit, size = gdf$size, logsz = F, method = "RegScore", pch = 21, cex=1,col = col.gp, bg = col.gp, xlab="(Log)Centroid Size") #if method = "PredLine" to plot prediction line
pred.score<-alo$RegScore #prediction scores for each species
preds <- shape.predictor(gpa$coords, x= gdf$size, Intercept = TRUE, predmin = min(gdf$size), predmax = max(gdf$size)) #estimating shape configurations based on allometric predictions

plotRefToTarget(preds$predmax, preds$predmin, mag=1, links = linkmandible) #shape related to small sizes

plotRefToTarget(preds$predmin, preds$predmax, mag=1, links = linkmandible)#shape related to large sizes

####compare evolutionary rates and morphological disparity
clades<-as.factor(datamandible$Clade) ##factor 1
names(clades)<-datamandible$species  

ecomorphs<-as.factor(datamandible$Diet_4)  ##factor 2
names(ecomorphs)<-datamandible$species

###SIZE##########
###  1 Caenolestinae; 2 Caluromyinae; 3 Didelphini; 4 Marmosini; 5 Thylamyini
#compare evolutionary rates between clades ->mandible size
mandsize<-compare.evol.rates(A=size, phy=tree, gp=clades,  iter = 9999)
ratesmandsize<-as.matrix(mandsh$sigma.d.gp)#create sigma obj

#####plot evol rates size  without Metachirini and Microbiotheriidae - insufficient data
plot(ratesmandsize[-(5:6),1], xlab="Clade",  
     ylab=expression(paste(sigma^2)), pch=21,cex=1.5)
abline(h = mandsh$sigma.d.all, lty=2) # plots the overal rate

#compare disparity between clades ->mandible size
dispmandsize<-morphol.disparity(size~1, group=clades,  iter = 9999,data=gdf)#pergroup
dispsize<-morphol.disparity(size~1, iter = 9999,data=gdf)#all disp is 0.214326
dispsizeclades<-as.matrix(dispmandsize$Procrustes.var)
plot(dispsizeclades[-(5:6),1], xlab="Clade",  
     ylab="Disparity", pch=21,cex=1.5, pt.bg = "black")
abline(h = 0.214326, lty=2) # plots the overal disp

#compare evolutionary rates between ecomorphs ->mandible size
###  1 Carnivores; 2 Frugivores; 3 Insectivores; 4 Omnivores
ecoevomandsize<-compare.evol.rates(A=size, phy=tree, gp=ecomorphs,  iter = 9999)
ratesmandsizeeco<-as.matrix(ecoevomandsize$sigma.d.gp)#create sigma obj
plot(ratesmandsizeeco[-4,1], xlab="Ecomorphs",  
     ylab="expression(paste(sigma^2))", pch=21,cex=1.5) #####plot evol rates size  without IC - insufficient data
abline(h = ecoevomandsize$sigma.d.all, lty=2) # plots the overal rate

#compare disparity between ecomorphs ->mandible size
ecodispmandsize<-morphol.disparity(size~1, group=ecomorphs,  iter = 9999,data=gdf)
dispsizeeco<-as.matrix(ecodispmandsize$Procrustes.var)
plot(dispsizeeco[-4,1], xlab="Ecomorph",  
     ylab="Disparity", pch=21,cex=1.5, pt.bg = "black")
abline(h = 0.214326, lty=2) # plots the overal disp


#comparing evolutionary rates between clades ->mandible shape
mandsh<-compare.evol.rates(A=gpa$coords, phy=tree, gp=clades,  iter = 9999)
ratesmandsh<-as.matrix(mandsh$sigma.d.gp)#create sigma obj

#####plot evol rates shape without Metachirini and Microbiotheriidae
###  1 Caenolestinae; 2 Caluromyinae; 3 Didelphini; 4 Marmosini; 5 Thylamyini
plot(ratesmandsh[-(5:6),1], xlab="Clade",  
     ylab=expression(paste(sigma^2)), pch=21,cex=1.5)
abline(h = mandsh$sigma.d.all, lty=2) # plots the overal rate

#compare disparity between clades ->mandible shape
dispmandshape<-morphol.disparity(coords~1, group=clades,  iter = 9999,data=gdf)
dispoverallshape<-morphol.disparity(coords~1, iter = 9999,data=gdf) #overall shape disp is 0.004957546
dispshclades<-as.matrix(dispmandshape$Procrustes.var)
#####plot shape disparity without Metachirini and Microbiotheriidae
plot(dispshclades[-(5:6),1], xlab="Clade",  
     ylab="Disparity", pch=21,cex=1.5, pt.bg = "black")
abline(h = 0.004957546, lty=2) # plots the overal rate
