library(Morpho) # Version 2.8
library(geomorph) # Version 4.0.7
library(phytools) # Version 2.1-1
library(castor)
library(RRphylo) # Version 2.8
library(viridisLite) # Version 0.4.0
library(landvR) # Version 0.3
library(dispRity) # Version 1.8
library(mvMORPH)
library(ggplot2) # Version 3.3.3
library(ggpubr) # Version 0.4.0
library(ggrepel)
library(ggdist)
library(ggthemes)
library(ggtree)
library(PupillometryR)
library(reshape2)
library(ggphylomorpho)
library(plotrix)

setwd("C:/Users/favie/Documents/Geomorph/NoHypno/Fintastic - Home/Images RAD")

########################################################################
###### GPA
########################################################################
## Read landmarks file

RJL <- readland.tps("BatSel_200524_FAna.TPS",specID = "imageID",
                    warnmsg = TRUE)

### Only Batoids, skip lines 35 and 36
Bat.set <- RJL[,,-c(127,128,143,144,156:158,223:235,
                    309:313,334,337,345,347,353,354,355)]

### Only Exant Batoids. Do not run if analysis is done with the previous line
Bat.set <- RJL[,,-c(127,128,143,144,156:158,223:235,
                    286:331,334,337,345,347,353,354,355,362)]

## Read curves
Jcurves<-as.matrix(read.csv("Sliders_BatSel_110224_FAna.CSV",
                            header = T, sep = ","))

### GPA
Bat.gpa<-gpagen(Bat.set,curves = Jcurves, ProcD = F)
plotOutliers(Bat.gpa$coords)

#Coracoid bar
CB.LM<-c(1,2,5,6,7,10:23,48:58,79:86)
#Propterygium
Pro.LM<-c(2:5,24:47)
#Metapterygium
Met.LM<-c(7:10,59:78)

## subset the landmarks
CB.coords<-Bat.set[CB.LM, , ]
PR.coords<-Bat.set[Pro.LM, , ]
MT.coords<-Bat.set[Met.LM, , ]

## read curves for each configuration
CBcurv<-as.matrix(read.csv("PG_Slide.csv",
                           header = T, sep = ","))
Procurv<-as.matrix(read.csv("Pro_Slide.csv",
                            header = T, sep = ","))
Metcurv<-as.matrix(read.csv("Met_Slide.csv",
                            header = T, sep = ","))

### GPA for each structure
CB.gpa<-gpagen(CB.coords,curves = CBcurv, ProcD = F, verbose = TRUE)
plotOutliers(CB.gpa$coords)

PR.gpa<-gpagen(PR.coords,curves = Procurv, ProcD = F, verbose = TRUE)
plotOutliers(PR.gpa$coords)

MT.gpa<-gpagen(MT.coords,curves = Metcurv, ProcD = F, verbose = TRUE)
plotOutliers(MT.gpa$coords)

########################################################################
###### Read Classifiers table
########################################################################
## Subset the classifiers. 
class.Raw <- read.csv("BatSel_Class_280524.csv",header = TRUE, sep = ",",
                      stringsAsFactors = TRUE)

## Dont run if all batoids are included
class.Raw<-class.Raw[class.Raw$EF=="Extant", ]

class.Bat<-class.Raw[class.Raw$SO=="Bat", ]
class.Bat <- droplevels(class.Bat)

## colour vector with fossils
p10 <- c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F","#C5C69D",
         "#93DA97","#38383B","#C24D67","#EF793D", "#A8AA6D")
names(p10) <- levels(class.Bat$Order)
p10.gp <- p10[match(class.Bat$Order,names(p10))]

## color vector only for extant
p4 <- c("#38383B","#C24D67","#EF793D","#A8AA6D")
names(p4) <- levels(class.Bat$Order)
p4.gp <- p4[match(class.Bat$Order,names(p4))]

########################################################################
###### Principal component analysis
########################################################################
## substitute the object if it's done with only extant
pc.bat.al <-gm.prcomp(Bat.gpa$coords)
pc.bat.cb <-gm.prcomp(CB.gpa$coords)
pc.bat.pr <-gm.prcomp(PR.gpa$coords)
pc.bat.mt <-gm.prcomp(MT.gpa$coords)

summary(pc.bat.al)

## simple plot
plot(pc.bat.al, axis1 = 1, axis2 = 2,pch = c(21,24)[as.numeric(class.Bat$EF)],
     cex = 1.5, bg = p10.gp)

## only extant
plot(pc.bat.al, axis1 = 1, axis2 = 2,pch = 21,
     cex = 1.5, bg = p4.gp)

### dataframe to use with ggplot
df.bat.al <- data.frame(Order = as.factor(class.Bat$Order),
                        SL = class.Bat$SpeLab,
                        FE = as.factor(class.Bat$EF),
                        PC1 = pc.bat.al$x[,1],
                        PC2 = pc.bat.al$x[,2])

df.bat.cb <- data.frame(Order = as.factor(class.Bat$Order),
                        SL = class.Bat$SpeLab,
                        FE = as.factor(class.Bat$EF),
                        PC1 = pc.bat.cb$x[,1],
                        PC2 = pc.bat.cb$x[,2])

df.bat.pr <- data.frame(Order = as.factor(class.Bat$Order),
                        SL = class.Bat$SpeLab,
                        FE = as.factor(class.Bat$EF),
                        PC1 = pc.bat.pr$x[,1],
                        PC2 = pc.bat.pr$x[,2])

df.bat.mt <- data.frame(Order = as.factor(class.Bat$Order),
                        SL = class.Bat$SpeLab,
                        FE = as.factor(class.Bat$EF),
                        PC1 = pc.bat.mt$x[,1],
                        PC2 = pc.bat.mt$x[,2])


sea3<- ggplot(df.bat.al, aes(PC1,PC2)) +
  geom_point(aes(colour = Order, fill = Order, shape = factor(FE)), size = 4,
             colour="black", stroke = 0.7) +
  scale_shape_manual(values=c(21,24))+
  geom_text_repel(mapping =  aes(label = SL, fontface = "italic"),
                  colour = "black", size = 1.5)+
  scale_fill_manual(values=c("#E0A6B3","#F39B6D","#31688EFF",
                             "#8D8D8F","#C5C69D",
                             "#93DA97","#38383B","#C24D67",
                             "#EF793D", "#A8AA6D"))+
  xlab("PC1: 54.97%") + ylab("PC2: 28.96%")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 0.5))

sea3

## same as above but with extant only
sea3<- ggplot(df.bat.al, aes(PC1,PC2)) +
  geom_point(aes(colour = Order, fill = Order, shape = factor(FE)), size = 4,
             colour="black", stroke = 0.7) +
  scale_shape_manual(values=c(21,24))+
  geom_text_repel(mapping =  aes(label = SL, fontface = "italic"),
                  colour = "black", size = 1.5)+
  scale_fill_manual(values=c("#38383B","#C24D67",
                             "#EF793D", "#A8AA6D"))+
  xlab("PC1: 54.97%") + ylab("PC2: 28.96%")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 0.5))

sea3

########################################################################
### Disparity by groups
########################################################################
### inc fossils
ArrayBT<-two.d.array(Bat.gpa$coords[,,c(1:330)])
ArrayCB<-two.d.array(CB.gpa$coords[,,c(1:330)])
ArrayMT<-two.d.array(MT.gpa$coords[,,c(1:330)])
ArrayPR<-two.d.array(PR.gpa$coords[,,c(1:330)])

### extant batoids only
ArrayBT<-two.d.array(Bat.gpa$coords[,,c(1:288)])
ArrayPG<-two.d.array(CB.gpa$coords[,,c(1:288)])
ArrayMT<-two.d.array(MT.gpa$coords[,,c(1:288)])
ArrayPR<-two.d.array(PR.gpa$coords[,,c(1:288)])

## substitute accordingly
DispAll <- dispRity.per.group(ArrayBT,
                              list(Rhinop=c(16,70,71,130:134,152:156,191:205,
                                            209:213,264,265,308,321,323),
                                   Myliob=c(1:7,46:60,65,66,74:93,102:109,
                                            123,127:129,140:151,160:162,206,
                                            217,225:263,309:311,314,317:319,
                                            322,325:327,329),
                                   Torped=c(28,29,61,62,94,110:122,135:139,
                                            218:224,312,315,316,324,328),
                                   Rajifo=c(8:15,17:27,30:45,63,64,67:69,72,73,
                                            95:101,124:126,157:159,163:190,
                                            207,208,214:216,307,313,320),
                                   Jr_Bato=c(266:272,277,278,303,304),
                                   Cr_Rhin=c(290:299,305,330),
                                   Cr_Scle=c(286,287,289,300:302),
                                   Cr_Raji=c(279,280,281),
                                   Eo_Myli=c(273:276,282:285,288)),
                              metric=c(sum, variances))

summary(DispAll)
test.dispRity(DispAll, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(DispAll)

## plot with ggplot
Rhinopristiformes<-as.matrix(t(DispAll$disparity$Rhinop[[2]]))
Myliobatiformes<-as.matrix(t(DispAll$disparity$Myliob[[2]]))
Torpediniformes<-as.matrix(t(DispAll$disparity$Torped[[2]]))
Rajiformes<-as.matrix(t(DispAll$disparity$Rajifo[[2]]))
JrBatoids<-as.matrix(t(DispAll$disparity$Jr_Bat[[2]]))
Cr_Raj<-as.matrix(t(DispAll$disparity$Cr_Raj[[2]]))
Cr_Scl<-as.matrix(t(DispAll$disparity$Cr_Scl[[2]]))
Cr_Rhi<-as.matrix(t(DispAll$disparity$Cr_Rhi[[2]]))
Eo_Myl<-as.matrix(t(DispAll$disparity$Eo_Myl[[2]]))

DispBGP <- as.data.frame(cbind(Rhinopristiformes,
                               Myliobatiformes,
                               Torpediniformes,
                               Rajiformes,
                               JrBatoids,
                               Cr_Scl,
                               Cr_Raj, 
                               Cr_Rhi,
                               Eo_Myl))

colnames(DispBGP) <- c("Rhino","Mylio",
                       "Torpe","Rajif",
                       "JrBat",
                       "Cr_Scl","Cr_Raj","Cr_Rhi","Eo_Myl")

ggplot(stack(DispBGP), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1.9, alpha=0.8,
              aes(colour = factor(ind))) + 
  scale_colour_manual(values = c("#EF793D","#38383B","#A8AA6D","#C24D67",
                                 "#93DA97","#31688EFF","#E0A6B3","#F39B6D",
                                 "#8D8D8F")) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  theme_minimal()+
  ylab("Sum Variances")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 0.5))+
  theme(legend.position = "none", axis.title.x = element_blank())

## Disparity for habitat groups 
## Only for extant batoids
DispEcBT <- dispRity.per.group(ArrayBT,
                               list(Deep=c(8:15,20:45,63,64,67:69,72,73,
                                           88,95:101,124,159,163:165,167,169:173,
                                           175:179,181:190,207,208,214,215,
                                           246,247,266),
                                    FreshW=c(65,66,87,92,140:151,225,236,277),
                                    Reef=c(4,48,49,85,86,93,94,115,123,
                                           127:134,162,191,209,210,217:220,
                                           230,231,234,237,238,251,264,265,
                                           268,271,275,276,280,283,287,288),
                                    Shelf=c(1:3,5:7,16:19,46,47,50:62,70,71,
                                            74:84,89:91,102:114,116:122,125,
                                            126,135:139,152:158,160,161,
                                            166,168,174,180,192:206,211:213,
                                            216,221:224,226:229,232,233,235,
                                            239:245,248:250,252:263,267,
                                            269,270,272:274,278,279,281,282,
                                            284:286)),
                               metric=c(sum, variances))
summary(DispEcBT)
test.dispRity(DispEcBT, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(DispEcBT)


## Swimming type
DispMovBT <- dispRity.per.group(ArrayBT,
                                list(Undul=c(8:15,17:27,30:60,63:69,
                                             72,73,87:93,95:101,123:151,
                                             157:159,163:190,207,208,
                                             214:217,225:266,268,269,272,
                                             276,277,279,284:286,288),
                                     Oscil=c(1:7,74:86,102:109,
                                             160:162,206,270,273,278,281),
                                     AxUnd=c(16,28,29,61,62,70,71,94,
                                             110:122,152:156,191:205,
                                             209:213,218:224,267,271,274,275,
                                             280,282,283,287)),
                                metric=c(sum, variances))

summary(DispMovBT)

plot(DispMovBT)
test.dispRity(DispMovBT, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")

### repeat with the different landmark subsets and merge the objects
### to create the boxplots as in the taxonomic groups plot

###################################################################
## Phylogentic analysis
###################################################################
########################################################################
###### Read Tree(s) and classifiers for the phylogeny
########################################################################
class.SP <- read.csv("Class_SpBatSel_Both_300524.csv",header = TRUE, sep = ",",
                     stringsAsFactors = TRUE)

class.BTre<-class.SP[class.SP$SO=="Bat",]

### do not run if including all batoids
class.BTre<-class.BTre[class.BTre$EF=="Extant",]

class.BTre <- droplevels(class.BTre)

### read calibrated tree and prune it for only batoids 
Phy.cal<-read.tree("PhyCalBatSel.tre")
plotTree(Phy.cal, fsize = 0.5)

sp.Pec <- as.vector(unique(class.BTre$TreeName2))
Phy.bat <- keep.tip(Phy.cal,tip = sp.Pec)
Phy.bat$root.time<-185
plotTree(Phy.bat, fsize = 0.5)

## coordinates as a 2D array
Bt.res<-two.d.array(Bat.gpa$coords)
CB.res<-two.d.array(CB.gpa$coords)
PR.res<-two.d.array(PR.gpa$coords)
MT.res<-two.d.array(MT.gpa$coords)

### average coordinates by species
means<-(aggregate(Bt.res~class.Bat$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Bat$TreeName2)
BatG<-arrayspecs(means,86,2)
means<-(aggregate(CB.res~class.Bat$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Bat$TreeName2)
Bt.sp.CB<-arrayspecs(means,38,2)
means<-(aggregate(PR.res~class.Bat$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Bat$TreeName2)
Bt.sp.PR<-arrayspecs(means,28,2)
means<-(aggregate(MT.res~class.Bat$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Bat$TreeName2)
Bt.sp.MT<-arrayspecs(means,24,2)

## For centroid size and explore possible adjustment with residuals
CsizeB<-(aggregate(Bat.gpa$Csize~class.Bat$TreeName2, FUN = mean))[,-1]
names(CsizeB)<-levels(class.Bat$TreeName2)

GDF.Phy <- geomorph.data.frame(shape1 = BatG,
                               size = log(CsizeB),
                               phy = Phy.bat) 

All.pDlm<- procD.lm(shape1~size, data = GDF.Phy, 
                    print.progress = TRUE)
summary(All.pDlm)

## including phylogeny
All.pgls <- procD.pgls(shape1~size, phy = phy, 
                       data = GDF.Phy, iter = 999) 
summary(All.pgls)
## size not significant
########################################################################
### Phylogenetically aligned PCA
########################################################################

### substitute the landmark data and tree deopending if all specimens
### or only extanant 
Bat.al.paca <- gm.prcomp(BatG, phy = Phy.bat,
                         align.to.phy = TRUE)
Bat.cb.paca <- gm.prcomp(Bt.sp.CB, phy = Phy.bat,
                         align.to.phy = TRUE)
Bat.pr.paca <- gm.prcomp(Bt.sp.PR, phy = Phy.bat,
                         align.to.phy = TRUE)
Bat.mt.paca <- gm.prcomp(Bt.sp.MT, phy = Phy.bat,
                         align.to.phy = TRUE)

summary(Bat.al.paca)

phy10 <- c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F","#C5C69D",
           "#93DA97","#38383B","#C24D67","#EF793D", "#A8AA6D")
class.BTre<-class.BTre[order(class.BTre$TreeName2), ]
names(phy10) <- levels(class.BTre$Order)
phy10.gp <- phy10[match(class.BTre$Order,names(phy10))]

phy4 <- c("#38383B","#C24D67","#EF793D","#A8AA6D")
class.BTre<-class.BTre[order(class.BTre$TreeName2), ]
names(phy4) <- levels(class.BTre$Order)
phy4.gp <- phy4[match(class.BTre$Order,names(phy4))]

## change xlab and ylab depending on the results from the summary
plot(Bat.al.paca, axis1 = 1, axis2 = 2,  pch=c(21,24)[as.numeric(class.BTre$EF)],
     cex = 1.5,
     xlab = "PC1: 80.99%", ylab = "PC2: 18.5 %", bg=phy10.gp,
     time.plot = FALSE, 
     phylo = TRUE,
     phylo.par = list(edge.color="grey50", anc.states =TRUE,
                      node.cex=0, node.labels = FALSE,
                      tip.labels = FALSE,
                      tip.txt.cex = 0.5))
## extant
plot(Bat.al.paca, axis1 = 1, axis2 = 2,  pch=c(21,24)[as.numeric(class.BTre$EF)],
     cex = 1.5,
     xlab = "PC1: 80.99%", ylab = "PC2: 18.5 %", bg=phy4.gp,
     time.plot = FALSE, 
     phylo = TRUE,
     phylo.par = list(edge.color="grey50", anc.states =TRUE,
                      node.cex=0, node.labels = FALSE,
                      tip.labels = FALSE,
                      tip.txt.cex = 0.5))

## phylogenetic signal, substitute for each data set 
PS.BT <- physignal.z(BatG, Phy.bat, iter = 999,
                     lambda = "all", PAC.no = 3)
summary(PS.BT)

## ggphylomorphospace data. Run only at once and substitute or make a
## data frame for each PACA
PCs<-cbind(Bat.al.paca$x[,1],Bat.al.paca$x[,2])

PCs<-cbind(Bat.pg.paca$x[,1],Bat.pg.paca$x[,2])
PCs<-cbind(Bat.pr.paca$x[,1],Bat.pr.paca$x[,2])
PCs<-cbind(Bat.mt.paca$x[,1],Bat.mt.paca$x[,2])

df <- data.frame(Order = as.factor(class.BTre$Order),
                 SL = class.BTre$SpeLab,
                 tip = class.BTre$TreeName2,
                 FE = as.factor(class.BTre$EF),
                 PC1 = PCs[,1],
                 PC2 = PCs[,2])

### I change line 29 of the function for "size = 0.1"
trace(ggphylomorpho, edit = T)

PMSpace <- ggphylomorpho(tree= Phy.bat, 
                         tipinfo = df,
                         labelvar = df$tip,
                         factorvar= df$Order, 
                         tree.alpha = 1, 
                         repel = TRUE)+
  geom_point(aes(x = df$PC1,y = df$PC2,
                 color = df$Order, fill = df$Order, 
                 shape = factor(df$FE)),
             size = 3.5,
             colour="black", stroke = 1) +
  
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F","#C5C69D",
                             "#93DA97","#38383B","#C24D67","#EF793D", "#A8AA6D"))+
  #scale_fill_manual(values=clade_colours)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="#FBF9FB"),
        legend.position="none",
        panel.grid.minor=element_blank())

PMSpace  + labs(x=paste('PAC1 = 80.99%'),
                y=paste('PAC2 = 18.5 %'))

## only extant ones
PMSpace <- ggphylomorpho(tree= Phy.bat, 
                         tipinfo = df,
                         labelvar = df$tip,
                         factorvar= df$Order, 
                         tree.alpha = 1, 
                         repel = TRUE)+
  geom_point(aes(x = df$PC1,y = df$PC2,
                 color = df$Order, fill = df$Order, 
                 shape = factor(df$FE)),
             size = 3.5,
             colour="black", stroke = 1) +
  
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c("#38383B","#C24D67","#EF793D", "#A8AA6D"))+
  #scale_fill_manual(values=clade_colours)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="#FBF9FB"),
        legend.position="none",
        panel.grid.minor=element_blank())

PMSpace  + labs(x=paste('PAC1 = 80.99%'),
                y=paste('PAC2 = 18.5 %'))

########################################################################
##### Disparity through time
########################################################################

## Only for data with fossils
MPCs<-Bat.al.paca$x[,c(1:3)]
NPCs<-Bat.al.paca$anc.x[,c(1:3)]

MPCs<-Bat.cb.paca$x[,c(1:3)]
NPCs<-Bat.cb.paca$anc.x[,c(1:3)]

MPCs<-Bat.pr.paca$x[,c(1:3)]
NPCs<-Bat.pr.paca$anc.x[,c(1:3)]

MPCs<-Bat.mt.paca$x[,c(1:3)]
NPCs<-Bat.mt.paca$anc.x[,c(1:3)]


DisDat<-rbind(MPCs,NPCs)
nodes<-paste("n",(1:193),sep="")
row.names(DisDat)[195:387]<-nodes
Phy.bat$node.label<-nodes

time_slices <- chrono.subsets(data   = DisDat,
                              tree   = Phy.bat,
                              method = "continuous",
                              model  = "proximity",
                              inc.nodes = T,
                              time   = 10)

boot_time_slices<-boot.matrix(time_slices, bootstraps = 1000)
disparity_time_slices<-dispRity(boot_time_slices,
                                metric=c(sum, variances))
plot.dispRity(disparity_time_slices, type = "continuous", 
              ylab="Sum of Variances",xlab ="Time (Mya)")


########################################################################
###### Modularity
########################################################################

## make partitions of landmarks for each module tested
## use the set which includes fossils

## three modules
mods<-c("A","A","B","B","A","A","A","C","C","A",
        "A","A","A","A","A","A","A","A","A","A","A","A","A",
        "B","B","B","B","B","B","B","B","B","B","B",
        "B","B","B","B","B","B","B","B","B","B","B","B","B",
        "A","A","A","A","A","A","A","A","A","A","A",
        "C","C","C","C","C","C","C","C","C",
        "C","C","C","C","C","C","C","C","C","C","C",
        "A","A","A","A","A","A","A","A")

## two modules PG_PR
modPG_Pr<-c("A","A","A","A","A","A","A","B","B","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A")

## two modules PG_Mt
modPG_Mt<-c("A","A","B","B","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A")

## Two modules Pr_Mt
modPr_Mt<-c("A","A","B","B","A","A","A","B","B","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A")

### Divide the data by taxonomic group
Raji<-BatG[,,c(11,12,17:24,27:30,32:35,39,40,53:60,65,66,
               87,88,90:92,94,95,109:111,130:132,136:146,
               166,167,169,170,172)]
Myli<-BatG[,,c(6:10,41:50,67:74,76:85,96:100,108,113,114,
               117:127,133,134,158,159,173,180:192)]
Torp<-BatG[,,c(26,51,52,86,101:107,115,116,176:179)]
Rhin<-BatG[,,c(1:4,13,14,61:64,128,129,147:151,153,154,
               162,163,164,193,194)]

### Divide the data by extant or extint batoids
BatoidsF <- BatG[,,c(31,152,155:157,161,171,174,175,36:38,
                     93,135,165,5,15,25,168,16,75,89)]
BatoidsE <- BatG[,,c(1:4,6:14,17:24,26:30,32:35,39:74,76:88,
                     90:92,94:134,136:151,153,154,158,
                     159,160,162:164,166,167,169,170,172,
                     173,176:194)]

### first the modularity hypothesis test
MT.AL <- phylo.modularity(BatG, mods,     Phy.bat, CI = TRUE, iter = 999)
MT.GP <- phylo.modularity(BatG, modPG_Pr, Phy.bat, CI = TRUE, iter = 999)
MT.GM <- phylo.modularity(BatG, modPG_Mt, Phy.bat, CI = TRUE, iter = 999)
MT.PM <- phylo.modularity(BatG, modPr_Mt, Phy.bat, CI = TRUE, iter = 999)

## compare effect sizes between different hypotheses
model.Z <- compare.CR(MT.AL, MT.GP, MT.GM, MT.PM, CR.null = TRUE)
model.Z

### Now each taxonomic subset. Trim a tree for each object. Example with Torpediniforms
Names <- dimnames(Torp)
species<-unique(Names[[3]])
Phy.Torp<-drop.tip(Phy.cal,Phy.cal$tip.label
                  [-match(species, Phy.cal$tip.label)])
plotTree(Phy.Torp, fsize = 0.5)

### extant orders selecting the supported hypothesis
MT.Raj <- phylo.modularity(Raji,mods,Phy.Raji,CI = TRUE,iter = 999)
MT.Myl <- phylo.modularity(Myli,mods,Phy.Myli,CI = TRUE,iter = 999)
MT.Rhi <- phylo.modularity(Rhin,mods,Phy.Rhin,CI = TRUE,iter = 999)
MT.Tor <- phylo.modularity(Torp,mods,Phy.Torp,CI = TRUE, iter = 999)

## extant vs extinct 
MT.BEx <- phylo.modularity(BatoidsE,mods,Phy.BEx,CI = TRUE, iter = 999)
MT.BFo <- phylo.modularity(BatoidsF,mods,Phy.BFo,CI = TRUE, iter = 999)


model.BO <- compare.CR(MT.Raj, MT.Myl, MT.Rhi, MT.Tor,CR.null = FALSE)
model.BO

model.BEF <- compare.CR(MT.BEx, MT.BFo, CR.null = FALSE)
model.BEF

########################################################################
###### Integration 
########################################################################

IT.Tor <- phylo.integration(Torp,partition.gp = mods,phy = Phy.Torp)
IT.Myl <- phylo.integration(Myli,partition.gp = mods,phy = Phy.Myli)
IT.Raj <- phylo.integration(Raji,partition.gp = mods,phy = Phy.Raji)
IT.Rhi <- phylo.integration(Rhin,partition.gp = mods,phy = Phy.Rhin)

IT.BEx <- phylo.integration(BatoidsE,partition.gp = mods,phy = Phy.BEx)
IT.BFo <- phylo.integration(BatoidsF,partition.gp = mods,phy = Phy.BFo)


summary(IT.Tor)

CINT.BO <- compare.pls(IT.Raj,IT.Myl,IT.Tor,IT.Rhi)
summary(CINT.BO)

CINT.BEF <- compare.pls(IT.BEx,IT.BFo)
summary(CINT.BEF)

globalIntegration(Torp)

globalIntegration(BatoidsF)


########################################################################
###### Phylogenetic comparative analyses
########################################################################

### Read all 100 calibrated trees
### The selection and calibration is found in the Script_SelBat_280524 file
### lines 275-298

Phy.mult<-read.tree(file="PhyCalBatSel1605.tre")

## trim the trees
sp.Pec <- as.vector(unique(class.BTre$TreeName2))

Bat.trees <- keep.tip.multiPhylo(Phy.mult, tip = sp.Pec)
plotTree(Bat.trees[[10]],fsize = 0.5)

### so far  trimming the tree to examine only batoid extant species


### create three matrices like this one for each model
Hab_L_ER <- matrix(NA, ncol= 2, nrow=length(Bat.trees))

## example with Equal Rates (ER) model

for(i in 1:length(Bat.trees)) {
  tree <- Bat.trees[[i]]
  ## Order the character to match the one from the input tree
  character <- class.BTre[match(Bat.trees[[i]]$tip.label,
                                class.BTre$TreeName2),]
  ## Select the trait from the data frame column (Col-10 habitat)
  ## (Col-11 Swimming type)
  tokens <- as.numeric(character[,10]) 
  ### change the rate model and Nstates accordingly
  fitER  <- fit_mk(trees = tree, Nstates = 4,
                    tip_states = tokens, 
                    rate_model="ER", Ntrials=2, 
                    Nthreads = 4, 
                    Nbootstraps = 100,
                    verbose = TRUE)
  ### Each new table change the object here "Hab_L_ARD[i,1]"
  ### and the parameter to extract
  Hab_L_ER[i,1] <- fitER$AIC
  Hab_L_ER[i,2] <- fitER$loglikelihood
}

### do the same with the remaining models and assign column names
colnames(Hab_L_ER)  <- c("AIC_ER_Hab" , "logLik_ER_Hab")
colnames(Hab_L_SYM) <- c("AIC_SYM_Hab", "logLik_SYM_Habb")
colnames(Hab_L_ARD) <- c("AIC_ARD_Hab", "logLik_ARD_Amb")

colnames(Mov_L_ER)  <- c("AIC_ER_Mov" , "logLik_ER_Mov")
colnames(Mov_L_SYM) <- c("AIC_SYM_Mov", "logLik_SYM_Mov")
colnames(Mov_L_ARD) <- c("AIC_ARD_Mov", "logLik_ARD_Mov")

Hab_Mod <- cbind(Hab_L_SYM, Hab_L_ARD, Hab_L_ER)
apply(Hab_Mod, 2, mean)
boxplot(Hab_Mod[,c(1,3,5)])
boxplot(Hab_Mod[,c(2,4,6)])

Mov_Mod <- cbind(Mov_L_SYM, Mov_L_ARD, Mov_L_ER)
apply(Mov_Mod, 2, mean)
boxplot(Mov_Mod[,c(1,3,5)])
boxplot(Mov_Mod[,c(2,4,6)])
### (lowest AIC is the selected model)

## map the trait to the many phylogenies, select the column
Amb<-as.factor(class.BTre$DepthCat3)
names(Amb)<-class.BTre$TreeName2

Mov<-as.factor(class.BTre$Mov)
names(Mov)<-class.BTre$TreeName2

### simmap with the selected model for each trait
map.er.mov <- make.simmap(Bat.trees,Mov,model="ER",
                          nsim = 1)
map.ard.amb <- make.simmap(Bat.trees,Amb,model="ARD",
                           nsim = 1)

## save the trees
save(map.er.mov,file="./Mov_map_ER.R")
save(map.ard.amb,file="./Amb_map_ARD.R")

load(file="./Hab_map_ARD.R")
load(file="./Mov_map_ER.R")

## simple plot to visualize
pd<-summary(map.ard.hab, plot=FALSE)

plot(pd,ftype="i", cex = 0.4, fsize = 0.5)

################################################
### from the previous analyses only with extant batoids
PCS <- Bat.al.paca$x

PCS <- Bat.cb.paca$x
PCS <- Bat.pr.paca$x
PCS <- Bat.mt.paca$x
### apparently 3 pcs help to describe better the data, up to 99% of the variation
MPCs<- PCS[,1:3]

Y <- MPCs
data <- list(Y=MPCs)

#######
### list for Sel landmark subset
data.sel=list(shape=MPCs,Amb=Amb, Mov=Mov)

#######################################
### evolutionary rates (extract sigma^2)  only with model "BMM"
## habitat mapping
simHabBat <- sapply(1:length(map.ard.hab), function(x) {
  fitBM <- mvgls(Y~1, data=data,
                 tree = map.ard.hab[[x]],
                 model = "BMM", method = "H&L",
                 error = TRUE, nbcores = 4L)
  fitBM$param
})

rowMeans(simHabBat)
boxplot(t(simHabBat))

## swimming mapping
simMovBat <- sapply(1:length(map.er.mov), function(x) {
  fitBM <- mvgls(Y~1, data=data,
                 tree = map.er.mov[[x]], 
                 model = "BMM", method = "H&L",
                 error = TRUE, nbcores = 4L)
  fitBM$param
})

rowMeans(simMovBat)
boxplot(t(simMovBat))

## save if necessary
save(simHabBat,file="./EVRates_BatAmb.R")
save(simMovBat,file="./EVRates_BatMov.R")

### repeat with the different landmark configurations

load(file="./EVRates_BatAmb.R")
load(file="./EVRates_BatMov.R")

## make a data frame, can be combined with the results from the other configurations
ER_all <- as.data.frame(cbind(t(simHabBat)))

ggplot(stack(ER_all), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1.9, alpha=0.8,
              aes(colour = factor(ind))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  scale_colour_manual(values = c("#0B0405FF","#357BA2FF","#3E356BFF","#ACBEB2")) + 
  theme_minimal()+
  ylab("Evolutonary Rates")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 0.5))+
  theme(legend.position = "none", axis.title.x = element_blank())+
  theme(legend.position = "none", axis.title.x = element_blank())

## same as above
ER_all <- as.data.frame(cbind(t(simMovBat)))
ggplot(stack(ER_all), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1.9, alpha=0.8,
              aes(colour = factor(ind))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  scale_colour_manual(values = c("#4686FBFF","#FABA39FF","#7A0403FF")) + 
  theme_minimal()+
  ylab("Evolutonary Rates")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 0.5))+
  theme(legend.position = "none", axis.title.x = element_blank())+
  theme(legend.position = "none", axis.title.x = element_blank())


#### mvgls, example with habitat type
fit.Amb<-mvgls(shape ~ Amb, data=data.sel,
               tree = Phy.bat, model="lambda", method="H&L",
               penalty = "RidgeArch",
               verbose = TRUE,
               nbcores = 4L)

summary(fit.Amb)

mult_Amb <- manova.gls(fit.Amb, nperm=999, type="II",
                       test="Pillai", verbose = TRUE,
                       nbcores = 4L)
mult_Amb
## repeat for the other configurations


fit.Amb<-mvgls(shape ~ Mov, data=data.sel,
               tree = Phy.bat, model="lambda", method="H&L",
               penalty = "RidgeArch",
               verbose = TRUE,
               nbcores = 4L)

summary(fit.Amb)

mult_Amb <- manova.gls(fit.Amb, nperm=999, type="II",
                       test="Pillai", verbose = TRUE,
                       nbcores = 4L)
mult_Amb
###########################
fit1  <- mvgls(Y~1, tree = Phy.bat, model="BM", penalty="RidgeArch",
               method = "PL-LOOCV", nbcores = 4L)
fit2  <- mvgls(Y~1, tree = Phy.bat, model="OU", penalty="RidgeArch",
               method = "PL-LOOCV", nbcores = 4L)
fit3  <- mvgls(Y~1, tree = Phy.bat, model="EB", penalty="RidgeArch",
               method = "PL-LOOCV" , nbcores = 4L)
fit4  <- mvgls(Y~1, tree = map.sym.tax[[1]], model="BMM", penalty="RidgeArch",
               method = "PL-LOOCV" , nbcores = 4L)

GIC(fit1); GIC(fit2); GIC(fit4); GIC(fit4)
GIC(fitBM); GIC(fitEB); GIC(fitOU)


simLA_Coords <- matrix(NA, ncol= 3, nrow=length(Bat.trees))

simLA <- for(i in 1:length(Bat.trees)) {
  Y <- Y[Bat.trees[[i]]$tip.label,]
  
  fitLA <- mvgls(Y~1,
                 tree = Bat.trees[[i]], penalty = "RidgeArch",
                 model = "lambda", method = "H&L",
                 error = TRUE, nbcores = 4L)
  simLAgic <- GIC(fitLA)
  simLA_Coords[i,1] <- simLAgic$LogLikelihood
  simLA_Coords[i,2] <- simLAgic$GIC
}

colnames(simBM_Coords) <- c("logLik_BM","GIC_BM","bias_BM")
colnames(simEB_Coords) <- c("logLik_EB","GIC_EB","bias_EB")
colnames(simOU_Coords) <- c("logLik_OU","GIC_OU","bias_OU")
colnames(simLA_Coords) <- c("logLik_LA","GIC_LA","bias_LA")

ECoord_Mod <- cbind(simBM_Coords, simEB_Coords, simOU_Coords,simLA_Coords)

apply(ECoord_Mod, 2, mean)
