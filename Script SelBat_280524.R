library(geomorph) # Version 4.0.7
library(phytools) # Version 2.1-1
library(RRphylo) # Version 2.8
library(viridisLite) # Version 0.4.0
library(landvR) # Version 0.3
library(dispRity) # Version 1.8
library(ggplot2) # Version 3.3.3
library(ggpubr) # Version 0.4.0
library(ggrepel)
library(ggdist)
library(ggthemes)
library(ggtree)
library(PupillometryR)
library(reshape2)
library(ggphylomorpho)


setwd("C:/Users/favie/Documents/Geomorph/NoHypno/Fintastic - Home/Images RAD")

########################################################################
###### Read landmarks
########################################################################
RJL <- readland.tps("BatSel_200524_FAna.TPS",specID = "imageID",
                    warnmsg = TRUE)

## read curves
Jcurves<-as.matrix(read.csv("Sliders_BatSel_110224_FAna.CSV",
                            header = T, sep = ","))

### GPA
BS.gpa<-gpagen(RJL,curves = Jcurves, ProcD = F)
BS.gpa<-gpagen(RJL,curves = Jcurves, ProcD = F, verbose = TRUE)
plotOutliers(BS.gpa$coords)
plotAllSpecimens(BS.gpa$coords, mean = TRUE)

#Pectoral girdle subset
PG.LM<-c(1,2,5,6,7,10:23,48:58,79:86)
#Propterygium subset
Pro.LM<-c(2:5,24:47)
#Metapterygium subset
Met.LM<-c(7:10,59:78)

PG.coords<-RJL[PG.LM, , ]
Pro.coords<-RJL[Pro.LM, , ]
Met.coords<-RJL[Met.LM, , ]

# read curves for each subset
PGcurv<-as.matrix(read.csv("PG_Slide.csv", header = T, sep = ","))
Procurv<-as.matrix(read.csv("Pro_Slide.csv",header = T, sep = ","))
Metcurv<-as.matrix(read.csv("Met_Slide.csv",header = T, sep = ","))

### GPA separated elements
PG.gpa<-gpagen(PG.coords,curves = PGcurv, ProcD = F, verbose = TRUE)
plotOutliers(PG.gpa$coords)

Pro.gpa<-gpagen(Pro.coords,curves = Procurv, ProcD = T, verbose = TRUE)
plotOutliers(Pro.gpa$coords)

Met.gpa<-gpagen(Met.coords,curves = Metcurv, ProcD = F, verbose = TRUE)
plotOutliers(Met.gpa$coords)


########################################################################
###### Read Classifiers table
########################################################################
class.Raw <- read.csv("BatSel_Class_280524.csv",header = TRUE, sep = ",",
                      stringsAsFactors = TRUE)

## colours vector
p14 <- c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F","#C5C69D","#8D76A6","#FDD57F","#93DA97",
         "#38383B","#492570","#87357D","#C24D67","#EF793D","#FCBD32","#A8AA6D")
names(p14) <- levels(class.Raw$Order)
p14.gp <- p14[match(class.Raw$Order,names(p14))]

########################################################################
###### PCA 
########################################################################
pc.all <-gm.prcomp(BS.gpa$coords)
pc.PG  <-gm.prcomp(PG.gpa$coords)
pc.Pro <-gm.prcomp(Pro.gpa$coords)
pc.Met <-gm.prcomp(Met.gpa$coords)

summary(pc.all)

plot(pc.all, axis1 = 1, axis2 = 2, pch = c(21,24)[as.numeric(class.Raw$EF)],
     cex = 1.5, bg = p14.gp)
legend("topright", legend= unique(class.Raw$Order), pch= c(19,17)[as.numeric(class.Raw$EF)],
       cex = 0.7, col = unique(p14.gp))


### data frame for plot with lables 
df <- data.frame(Order = as.factor(class.Raw$Order),
                 SL = class.Raw$SpeLab,
                 FE = as.factor(class.Raw$EF),
                 PC1 = pc.all$x[,1],
                 PC2 = pc.all$x[,2])
df.PG <- data.frame(Order = as.factor(class.Raw$Order),
                    SL = class.Raw$SpeLab,
                    FE = as.factor(class.Raw$EF),
                    PC1 = pc.PG$x[,1],
                    PC2 = pc.PG$x[,2])
df.Pr <- data.frame(Order = as.factor(class.Raw$Order),
                    SL = class.Raw$SpeLab,
                    FE = as.factor(class.Raw$EF),
                    PC1 = pc.Pro$x[,1],
                    PC2 = pc.Pro$x[,2])
df.Mt <- data.frame(Order = as.factor(class.Raw$Order),
                    SL = class.Raw$SpeLab,
                    FE = as.factor(class.Raw$EF),
                    PC1 = pc.Met$x[,1],
                    PC2 = pc.Met$x[,2])


sea3<- ggplot(df, aes(PC1,PC2)) +
  geom_point(aes(colour = Order, fill = Order, shape = factor(FE)), size = 4,
             colour="black", stroke = 0.7) +
  scale_shape_manual(values=c(21,24))+
  geom_text_repel(mapping =  aes(label = SL, fontface = "italic"),
                  colour = "black", size = 1.5)+
  scale_fill_manual(values=c("#E0A6B3","#F39B6D","#31688EFF",
                             "#8D8D8F","#C5C69D","#8D76A6",
                             "#FDD57F","#93DA97","#38383B",
                             "#492570","#87357D","#C24D67",
                             "#EF793D","#FCBD32","#A8AA6D"))+
  xlab("PC1: 59.51%") + ylab("PC2: 22.3%")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 0.5))

sea3


########################################################################
##### Procrustes variance plots 
########################################################################

# select each configuration
procrustes_2d_array <- geomorph::two.d.array(BS.gpa$coords)
ordination <- stats::prcomp(procrustes_2d_array)

PC1_variation <- variation.range(BS.gpa, axis = 1,
                                 ordination = ordination, type = "spherical")
PC2_variation <- variation.range(BS.gpa, axis = 2,
                                 ordination = ordination, type = "spherical")

## nrow depends on the number of landmarks for each configuration
hyp_P1min <- as.matrix(pc.all$shapes$shapes.comp1$min[1:86,],
                       nrow = 86, ncol = 2)
hyp_P1max <- as.matrix(pc.all$shapes$shapes.comp1$max[1:86,],
                       nrow = 86, ncol = 2)
hyp_P2min <- as.matrix(pc.all$shapes$shapes.comp2$min[1:86,],
                       nrow = 86, ncol = 2)
hyp_P2max <- as.matrix(pc.all$shapes$shapes.comp2$max[1:86,],
                       nrow = 86, ncol = 2)

procrustes.var.plot(hyp_P1min, hyp_P1min, col = inferno, 
                    magnitude = 1, pt.size = 2.5,
                    col.val = PC1_variation[, "radius"], labels = FALSE)


########################################################################
### Disparity by groups
########################################################################
### inc fossils
Array2D<-two.d.array(BS.gpa$coords[,,c(1:362)])
ArrayPG<-two.d.array(PG.gpa$coords[,,c(1:362)])
ArrayMT<-two.d.array(Met.gpa$coords[,,c(1:362)])
ArrayPR<-two.d.array(Pro.gpa$coords[,,c(1:362)])

DispAll <- dispRity.per.group(Array2D,
                              list(Rhinop=c(16,70,71,132:136,159:163,198:212,
                                            216:220,284,285,333,350,352),
                                   Myliob=c(1:7,46:60,65,66:93,102:109,123,129:131,
                                            142,145:155,167:169,213,237,245:283,335,
                                            336,338,341,344,346,348,351,357:359,361),
                                   Torped=c(28,29,61,62,94,110:122,137:141,238:244,
                                            339,342,343,356,360),
                                   Rajifo=c(8:15,17:27,30:45,63,64,67:69,72,73,95:101,
                                            124:126,164:166,170:197,214,215,221,222,
                                            236,332,340,349),
                                   Orecto=c(127,128,334,337,345),
                                   Pristi=c(143,144,156:158,347),
                                   Squati=c(223:235,353:355),
                                   Cr_Raj=c(299:301),
                                   Cr_Scl=c(306,307,314,325:327),
                                   Cr_Rhi=c(315:324,330,362),
                                   Eo_Myl=c(293:296,302:305,308),
                                   Jr_Squ=c(311:313),
                                   Jr_Bat=c(286:292,297,298,
                                            328,329)),
                              metric=c(sum, variances))


summary(DispAll)
test.dispRity(DispAll, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(DispAll)

### extract results for each group
Rhinopristiformes<-as.matrix(t(DispAll$disparity$Rhinop[[2]]))
Myliobatiformes<-as.matrix(t(DispAll$disparity$Myliob[[2]]))
Torpediniformes<-as.matrix(t(DispAll$disparity$Torped[[2]]))
Rajiformes<-as.matrix(t(DispAll$disparity$Rajifo[[2]]))
Orectolobiformes<-as.matrix(t(DispAll$disparity$Orecto[[2]]))
Pristiophoriformes<-as.matrix(t(DispAll$disparity$Pristi[[2]]))
Squatiniformes<-as.matrix(t(DispAll$disparity$Squati[[2]]))
JrBatoids<-as.matrix(t(DispAll$disparity$Jr_Bat[[2]]))
Cr_Raj<-as.matrix(t(DispAll$disparity$Cr_Raj[[2]]))
Cr_Scl<-as.matrix(t(DispAll$disparity$Cr_Scl[[2]]))
Cr_Rhi<-as.matrix(t(DispAll$disparity$Cr_Rhi[[2]]))
Eo_Myl<-as.matrix(t(DispAll$disparity$Eo_Myl[[2]]))
Jr_Squ<-as.matrix(t(DispAll$disparity$Jr_Squ[[2]]))

## combine in a data frame
DispGP <- as.data.frame(cbind(Rhinopristiformes,
                              Myliobatiformes,
                              Torpediniformes,
                              Rajiformes,
                              Orectolobiformes,
                              Pristiophoriformes,
                              Squatiniformes,
                              JrBatoids,
                              Jr_Squ,
                              Cr_Scl,
                              Cr_Raj, 
                              Cr_Rhi,
                              Eo_Myl))

colnames(DispGP) <- c("Rhino","Mylio",
                      "Torpe","Rajif",
                      "Orect","Prist",
                      "Squat","JrBat","Jr_Squ",
                      "Cr_Scl","Cr_Raj","Cr_Rhi","Eo_Myl")

ggplot(stack(DispGP), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1.9, alpha=0.8,
              aes(colour = factor(ind))) + 
  scale_colour_manual(values = c("#EF793D","#38383B","#A8AA6D","#C24D67",
                                 "#492570","#87357D","#FCBD32","#93DA97",
                                 "#FDD57F","#31688EFF","#E0A6B3","#F39B6D","#8D8D8F")) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  theme_minimal()+
  ylab("Sum Variances")+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 0.5))+
  theme(legend.position = "none", axis.title.x = element_blank())


###################################################################
## Phylogentic analysis
###################################################################

########################################################################
###### Read Tree(s) and classifiers for the phylogeny
########################################################################
class.SP <- read.csv("Class_SpBatSel_Both_300524.csv",header = TRUE, sep = ",",
                     stringsAsFactors = TRUE)
set.seed(42)

## 420 trees
Phy.mult<-read.tree(file="BatTreeBRL150424.tre")
## trim the trees
sp.Pec <- as.vector(unique(class.SP$TreeName2))
Chond3D.trees <- keep.tip.multiPhylo(Phy.mult, tip = sp.Pec)

## ages for tips and nodes
age.tips<-read.table("FADLAD.txt",header = FALSE, sep = ",",
                     row.names = 1)
age.nodes<-read.table("NodeAges.txt",header = FALSE, sep = ",",
                      row.names = 1)

## subset 100 trees
random.trees<-sample(Chond3D.trees,size=100)

## calibrate
treeS3 <- lapply(random.trees, scaleTree,tip.ages = age.tips, 
                 node.ages = age.nodes)

## save the 100 calibrated random trees
write.tree(treeS3,file = "PhyCalBatSel1605.tre")
Phy.mult<-read.tree(file="PhyCalBatSel1605.tre")

## single best score tree
Phy.cal<-read.tree("PhyCalBatSel.tre")
Phy.cal$root.time<-256


##### average landmarks by species
Al.res<-two.d.array(BS.gpa$coords)
PG.res<-two.d.array(PG.gpa$coords)
PR.res<-two.d.array(Pro.gpa$coords)
MT.res<-two.d.array(Met.gpa$coords)


means<-(aggregate(Al.res~class.Raw$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Raw$TreeName2)
ElasmoG<-arrayspecs(means,86,2)
means<-(aggregate(PG.res~class.Raw$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Raw$TreeName2)
PG.sp<-arrayspecs(means,38,2)
means<-(aggregate(PR.res~class.Raw$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Raw$TreeName2)
PR.sp<-arrayspecs(means,28,2)
means<-(aggregate(MT.res~class.Raw$TreeName2, FUN = mean))[,-1]
rownames(means)<-levels(class.Raw$TreeName2)
MT.sp<-arrayspecs(means,24,2)

## this one

CsizeB<-(aggregate(BS.gpa$Csize~class.Raw$TreeName2, FUN = mean))[,-1]
names(CsizeB)<-levels(class.Raw$TreeName2)

GDF.Phy <- geomorph.data.frame(shape1 = ElasmoG,
                               size = log(CsizeB),
                               phy = Phy.cal) 

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

All.paca <- gm.prcomp(ElasmoG, phy = Phy.cal,
                      align.to.phy = TRUE)
PG.paca <- gm.prcomp(PG.sp, phy = Phy.cal,
                     align.to.phy = TRUE)
PR.paca <- gm.prcomp(PR.sp, phy = Phy.cal,
                     align.to.phy = TRUE)
MT.paca <- gm.prcomp(MT.sp, phy = Phy.cal,
                     align.to.phy = TRUE)

summary(All.paca)
summary(PG.paca)
summary(PR.paca)
summary(MT.paca)


phy14 <- c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F",
           "#C5C69D","#93DA97","#8D76A6","#FDD57F",
           "#38383B","#492570","#87357D","#C24D67",
           "#EF793D","#FCBD32","#A8AA6D")
class.SP<-class.SP[order(class.SP$TreeName2), ]
names(phy14) <- levels(class.SP$Order)
phy14.gp <- phy14[match(class.SP$Order,names(phy14))]

### substitute arguments accordingly
plot(All.paca, axis1 = 1, axis2 = 2,  pch=c(21,24)[as.numeric(class.SP$EF)],
     cex = 1.5,
     xlab = "PC1: 80.99%", ylab = "PC2: 18.5 %", bg=phy14.gp,
     time.plot = FALSE, 
     phylo = TRUE,
     phylo.par = list(edge.color="grey50", anc.states =TRUE,
                      node.cex=0, node.labels = FALSE,
                      tip.labels = FALSE,
                      tip.txt.cex = 0.5))

### phylogenetic signal, change the landmark configuration 
PS.All <- physignal.z(ElasmoG, Phy.cal, iter = 999,
                     lambda = "all", PAC.no = 3)

summary(PS.All)


### phylomorphospace with ggplot
### I change line 29 of the function for "size = 0.1"
trace(ggphylomorpho, edit = T)

PCs<-cbind(All.paca$x[,1],All.paca$x[,2])

PCs<-cbind(PG.paca$x[,1],PG.paca$x[,2])
PCs<-cbind(PR.paca$x[,1],PR.paca$x[,2])
PCs<-cbind(MT.paca$x[,1],MT.paca$x[,2])

df <- data.frame(Order = as.factor(class.SP$Order),
                 SL = class.SP$SpeLab,
                 tip = class.SP$TreeName2,
                 FE = as.factor(class.SP$EF),
                 PC1 = PCs[,1],
                 PC2 = PCs[,2])

PMSpace <- ggphylomorpho(tree= Phy.cal, 
                         tipinfo = df,
                         labelvar = df$tip,
                         factorvar= df$Order, 
                         tree.alpha = 1, 
                         repel = TRUE)+
  geom_point(aes(x = df$PC1,y = df$PC2,
                 color = df$Order, fill = df$Order, 
                 shape = factor(df$FE)),
             size = 3.5,
             colour="black", stroke = 0.7) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F",
                             "#C5C69D","#93DA97","#8D76A6","#FDD57F",
                             "#38383B","#492570","#87357D","#C24D67",
                             "#EF793D","#FCBD32","#A8AA6D"))+
  #scale_fill_manual(values=clade_colours)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="#FBF9FB"),
        legend.position="none",
        panel.grid.minor=element_blank())

### set the values depenting on the results from the summary
PMSpace  + labs(x=paste('PAC1 = 94.00 %'),
                y=paste('PAC2 = 5.59 %'))

###############################################################
##### RRphylo
###############################################################
cc<-2/parallel::detectCores()

## select the PACs necessary for each configuration
PCS <- All.paca$x

PCS <- PG.paca$x
PCS <- PR.paca$x
PCS <- MT.paca$x

MPCs<- PCS[,1:3]

## ridge arch regression
RRphylo(tree=Phy.cal,y=MPCs)->RR

plotRR(RR,y=MPCs,multivariate="rates")->pRRmulti

pRRmulti$plotRRrates(variable=1,
                     tree.args=list(no.margin=TRUE,type="fan",
                                    edge.width=4,show.tip.label = FALSE),
                     color.pal=inferno,
                     colorbar.args = list(x="bottomleft",
                                          width = 0.3,
                                          labs.adj=0.4,
                                          xpd=TRUE,
                                          tck.pos="out"))


nodes<-c(406,403,401,216,384,322,249)
labels<-c("Selachii","Jurassic Batoids","Sclerorhyncoidei","Rhinopristiformes",
          "Torpediniformes","Rajiformes","Myliobatiformes")

for(i in 1:length(nodes))
  arc.cladelabels(tree=Phy.cal,labels[i],nodes[i],
                  orientation=if(labels[i]%in%c("Jurassic Batoids",
                                                "Sclerorhyncoidei"))
                    "horizontal" else "curved",mark.node=FALSE,
                  lwd=5, cex = 1)

### 406=Sel; 403=Jbat; 401=Scle; 216 Rhin; 
### 384 Torp; 322 Raji; 249 Mylio; 213 Bats
search.shift(RR=RR,status.type="clade",
             node=c(406,213,403,401,
                    216,384,322,249))->SSauto

plotShift(RR=RR,SS=SSauto)->plotSS
addShift(SS=SSauto)
plotSS$plotClades()

### phylogenetic uncertainty, test the shifts with the 100 random trees
overfitRR(RR=RR,y=MPCs, phylo.list = random.trees,
          shift.args = list(node=rownames(SSauto$single.clades)),
          nsim=10,clus=cc)->orr.ss
orr.ss$shift.results$clade

########### similarity/convergence

## make a distance matrix for the Procrustes ditances
Fins.PD<-as.matrix(dist(two.d.array(gpagen(ElasmoG, curves = Jcurves,
                                           ProcD = F)$coords),
                        diag = FALSE,upper = FALSE))

Fins.PD <- as.dist(Fins.PD)
jaw.clust <- hclust(Fins.PD, method = "average")

plot(jaw.clust)
## make the distance dendrogram as phylo tree
hcd <- as.phylo(jaw.clust)
plotTree(hcd, fsize = 0.5)

## tanglegram with the phylogeny and distance tree
obj<-cophylo(Phy.cal,hcd,  rotate = TRUE)
plot(obj,  fsize=0.2, link.lwd = 2)

class.phyOrd <- class.SP[match(Phy.cal$tip.label, 
                               class.SP$TreeName2),]

phy14 <- c("#E0A6B3","#F39B6D","#31688EFF","#8D8D8F",
           "#C5C69D","#93DA97","#8D76A6","#FDD57F",
           "#38383B","#492570","#87357D","#C24D67",
           "#EF793D","#FCBD32","#A8AA6D")

names(phy14) <- levels(class.phyOrd$Order)
phy14.gp <- phy14[match(class.phyOrd$Order,names(phy14))]
plot(obj, link.col = phy14.gp, fsize=0.3, link.lwd = 2)

### Convergence by clades-nodes
### 246 Zapt; 400 Plat; 412 squ;  403 Jbat; 216 Rhin
### 401 Scl; 237 Pris; 406=Sel; 213 Bats; 329 Cyclobatis
### 281 AsteroHelio 410=Hypno

## repeat with the selected nodes of interest
search.conv(RR=RR, y=MPCs, 
            nodes=c(329,281),clus=cc)->sc_clade
plotConv(SC=sc_clade, y=MPCs, variable=1, RR = RR)->plotSC
par(mfrow=c(1,2))
plotSC$plotPChull()
plotSC$plotTraitgram()

dev.off()

## convergence by custer groups from the tanglegram
traits_vector <- rep("nostate", 211)
names(traits_vector)<-rownames(MPCs)

JBat<-c("Spathobatus","Aellopobatis_bavarica",
        "Asterodermus_platypterus","Belemnobatis_sismondae")
FGuit<-c("Rhinobatos_whitfieldy", "Rhinobatus_latus", 
         "Rhinobatos_hakelensis", "Stahlraja_sertanensis",
         "Tlalocbatos")
EGuit<-c("Rhinobatos_glaucostigma", "Glaucostegus_thouin",
         "Platyrhinoidis_triseriata", "Zapteryx_exasperata",
         "Zapteryx_xyster")

traits_vector[which(names(traits_vector)%in%JBat)]<-"Jbatoid"
traits_vector[which(names(traits_vector)%in%FGuit)]<-"F_Guit"
traits_vector[which(names(traits_vector)%in%EGuit)]<-"E_Guit"

conv.test <- search.conv(tree=Phy.cal,y= MPCs, declust = TRUE,
                         state = traits_vector, clus = cc)
conv.test
plotConv(SC=conv.test, y=MPCs, variable=1, state=traits_vector)->plotSC_state

# plots
par(mfrow=c(1,2))
plotSC_state$plotPChull()
plotSC_state$plotPolar()


####################################################################
## Disparity through time
####################################################################

MPCs<-All.paca$x[,c(1:3)]
NPCs<-All.paca$anc.x[,c(1:3)]

DisDat<-rbind(MPCs,NPCs)
nodes<-paste("n",(1:210),sep="")
row.names(DisDat)[212:421]<-nodes
Phy.cal$node.label<-nodes

time_slices <- chrono.subsets(data   = DisDat,
                              tree   = Phy.cal,
                              method = "continuous",
                              model  = "proximity",
                              inc.nodes = T,
                              time   = 10)

boot_time_slices<-boot.matrix(time_slices, bootstraps = 1000)
disparity_time_slices<-dispRity(boot_time_slices,
                                metric=c(sum, variances))
plot.dispRity(disparity_time_slices, type = "continuous", 
              ylab="Sum of Variances",xlab ="Time (Mya)")

summary(disparity_time_slices)
###############################################################
### Modularity and Integration 
## make partitions

mods<-c("A","A","B","B","A","A","A","C","C","A",
        "A","A","A","A","A","A","A","A","A","A","A","A","A",
        "B","B","B","B","B","B","B","B","B","B","B",
        "B","B","B","B","B","B","B","B","B","B","B","B","B",
        "A","A","A","A","A","A","A","A","A","A","A",
        "C","C","C","C","C","C","C","C","C",
        "C","C","C","C","C","C","C","C","C","C","C",
        "A","A","A","A","A","A","A","A")

modPG_Pr<-c("A","A","A","A","A","A","A","B","B","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A")

modPG_Mt<-c("A","A","B","B","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A")

modPr_Mt<-c("A","A","B","B","A","A","A","B","B","A",
            "A","A","A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B",
            "B","B","B","B","B","B","B","B","B","B","B",
            "A","A","A","A","A","A","A","A")

Sharks1 <- ElasmoG[,,c(28,58,114,118,122,133:135,138,142,143,
                        182:187)]
Batoids1 <- ElasmoG[,,c(1:27,29:57,59:113,115:117,119:121,123:132,
                        136,137,139:141,144:181,188:211)]

M<-mshape(jaws3d.gpa$coords)
n<-2
p3 <- viridis(n)
names(p3) <- levels(as.factor(modPr_Mt))
col.land <- p3[match(modPr_Mt, names(p3))]
plotRefToTarget(M,M,  method = "point",
                gridPars = gridPar(pt.bg = c(col.land),
                                   pt.size = 2.5,
                                   tar.pt.bg = c(col.land),
                                   tar.pt.size = 3))


### test which hypothesis is better supported
MT.AL <- phylo.modularity(ElasmoG, mods,     Phy.all, CI = TRUE, iter = 999)
MT.GP <- phylo.modularity(ElasmoG, modPG_Pr, Phy.all, CI = TRUE, iter = 999)
MT.GM <- phylo.modularity(ElasmoG, modPG_Mt, Phy.all, CI = TRUE, iter = 999)
MT.PM <- phylo.modularity(ElasmoG, modPr_Mt, Phy.all, CI = TRUE, iter = 999)

## compare effect sizes
model.Z <- compare.CR(MT.AL, MT.GP, MT.GM, MT.PM, CR.null = TRUE)
model.Z

## in this scrip only test for the difference between batoids and sharks
MT.Bat <- phylo.modularity(Batoids1,mods,Phy.BT,CI = TRUE, iter = 999)
MT.Sha <- phylo.modularity(Sharks1,mods,Phy.Sha,CI = TRUE, iter = 999)

model.BS <- compare.CR(MT.Bat, MT.Sha, CR.null = FALSE)
model.BS


###############################################################
### Integration 

IT.Bat <- phylo.integration(Batoids1,partition.gp = mods,phy = Phy.BT)
IT.Sha <- phylo.integration(Sharks1,partition.gp = mods,phy = Phy.Sha)

CINT.BS <- compare.pls(IT.Bat,IT.Sha)
summary(CINT.BS)

globalIntegration(Shareks1)
globalIntegration(Batoids1)

