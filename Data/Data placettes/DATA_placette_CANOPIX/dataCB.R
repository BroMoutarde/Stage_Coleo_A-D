
# GLMM
#selection de la loi de proba la mieux ajustée à la distribution des données
library(fitdistrplus)
# sélection de la loi d'ajustement : poisson, gaussien, negative binomiale ou lognormal

descdist(Rscolytexoticfemer$rs.scol.non.ambrosia,discrete=TRUE,boot=1001)
plot(fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"norm"))
fitnb<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"nbinom")
fitp<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"pois")
fitn<-fitdist(Rscolytexoticflemer$rs.scol.non.ambrosia,"norm")
fitlnorm<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"lnorm")
fitexp<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"exp")
fitgeom<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"geom")
fitbeta<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"beta")
fitunif<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"unif")
fitg<-fitdist(Rscolytexoticfemer$rs.scol.non.ambrosia,"gamma",method="mme")
gofstat(fitnb)$chisqpvalue
gofstat(fitp)$chisqpvalue
gofstat(fitn)$chisqpvalue
gofstat(fitlnorm)$chisqpvalue
gofstat(fitexp)$chisqpvalue
gofstat(fitgeom)$chisqpvalue
gofstat(fitunif)$chisqpvalue
gofstat(fibeta)$chisqpvalue
gofstat(fitg)$chisqpvalue
#après log+1 transformation de la variable
descdist(log1p(Rscolytexoticfemer$rs.scol.non.ambrosia),discrete=TRUE,boot=1001)
plot(fitdist(log1p(Rscolytexoticfemer$rs.scol.non.ambrosia),"norm"))
fitn<-fitdist(log1p(Rscolytexoticfemer$rs.scol.non.ambrosia),"norm")
gofstat(fitn)$chisqpvalue
# la loi la mieux ajustée correspond à la p-value la plus élevée

# voir aussi : tweedie

# main predictor (stratum) used as either a linear term or a quadratic term to be better fit for our data set. 
# add forest and stand as nested random effects on the intercept in the mixed models to account for repeated measurements and the spatial configuration of the sampling design. 
# Since some traps were not continuously functional,for example, because they fell from the tree, add an offset of the log-number of effective traps across the entrapment season.

#addition de l'autocorrélation spatiale
#Avec ACS
pos<-numFactor((as.numeric(Rscolytexoticfcircul$longitude)),(as.numeric(Rscolytexoticfcircul$latitude)))#va faire des longitudes /latitudes un couple
ID <- factor(rep(1, 3887))#une placette = un identifiant
t3<-Sys.time()
model3<-glmmTMB(rs.scol.non.ambrosia~occur.germanus+(1|cluster)+(1|idtrapyear)+exp(pos+0|ID),family="nbinom1",data=Rscolytexoticfcircul,na.action=na.fail)

# #magnitude
#récupérer les valeurs de Estimate et Std. Error pour la variable d’intérêt X

Estimate_mod<-rnorm(10000,mean=0.3200,sd=0.1480)
DX<-exp(Estimate_mod*10)-1 #pour un delta de X de 10 (%)
mean(DX)
quantile(DX, c(0.01, 0.99))

#diagnostic glmm : ajustement, surdispersion, outliers
library(DHARMa)
sim<-simulateResiduals(model3)
testUniformity(sim)#le qqplot est issu de cette commande 
#KS Test p-value # Dispersion test  # Outliers 
testOutliers(sim)

#use the differences in AICc scores (function AICc, AICCmodavg R-package) to compare the fit among models.

# account for the increased type I error risk due to multiple testing (n parameters) with the Bonferroni correction




#############################################################################
#############################################################################
#############################################################################
# Analyses multivariées des variations de composition de communautés

################### ORDINATIONS CAP

distinfra<-vegdist(buprbuchetotmat[,6:15], method="bray",na.rm = TRUE)
bupr.cap<-capscale(distbupr~as.numeric(depefeu2)+as.factor(depetree3)+as.factor(annee),buprbuchetotmat)

#Pour conna?tre la contribution ? l'inertie de chaque variable X du mod?le et tester leur significativit?, construire tous les mod?les ?l?mentaires ? pour la contribution marginale , utiliser la fonction anova (by=margin)
library(vegan)
#toto est la matrice esp?ces-relev?s
#tabtrap est la matrice relev?s-environnement
#chaque variable a une contribution intrins?que (dite marginale) ? l'inertie du nuage et une contribution conjointe avec d'autres variables ; la contribution totale de chaque variable est la somme de (intrins?que + conjointe)
artbmmh2.cap<-capscale(toto2~altitude+site+forest.type+mixture+LSDW+LLDW+LT+dcav+dsap+dfun+dperdw+open,tabtrap2,dist="jaccard")
anova(artbmmh2.cap, by="margin")

#La colonne Var pour chaque ligne indique la contribution intrins?que de chaque variable ? l'inertie totale.
#Pour la ligne Residual, la colonne Var indique l'inertie non contrainte du nuage.
#Pour calculer l'inertie contrainte dans la nuage par l'ensemble des variables = inertie totale - inertie non contrainte
#L'inertie totale et la contribution totale d'une variable ? l'inertie (incluant les contributions partag?es avec les autres variables), sont d?duites de chaque mod?le ?l?mentaire :
#Test de la significativit? de l'effet de cette variable sur la composition de l'assemblage sp?cifique
artbmmh2.cap<-capscale(toto2~altitude,tabtrap2,dist="jaccard")
anova.cca(artbmmh2.cap)
artbmmh2.cap<-capscale(toto2~site,tabtrap2,dist="jaccard")
anova.cca(artbmmh2.cap)
artbmmh2.cap<-capscale(toto2~forest.type,tabtrap2,dist="jaccard")
anova.cca(artbmmh2.cap)
#...
#L'inertie totale du nuage est l'inertie du mod?le + l'inertie r?siduelle (ici = 17.265+274.405)
#Le % de l'inertie contrainte expliqu?e par cette variable : =17,26*100/155 = ?? %
#Il suffit ensuite de reconstituer un tableau avec l'inertie totale expliqu?e par chaque variable, son %, et sa significativit?	
#La colonne Var pour chaque ligne indique la contribution intrins?que de chaque variable ? l'inertie totale.
#Pour la ligne Residual, la colonne Var indique l'inertie non contrainte du nuage.
#Pour calculer l'inertie contrainte dans la nuage par l'ensemble des variables = inertie totale - inertie non contrainte
#L'inertie totale et la contribution totale d'une variable ? l'inertie (incluant les contributions partag?es avec les autres variables), sont d?duites de chaque mod?le ?l?mentaire :
#Test de la significativit? de l'effet de cette variable sur la composition de l'assemblage sp?cifique

############################ ordination NMDS inter-strates
dist.best<-vegdist(traitTREM[,4:42], method="bray",na.rm = TRUE)
nmds<-metaMDS(traitTREM[,4:42], distance = "bray")
ordiplot (nmds, display = 'sites', type = 'n')
orditorp (nmds, display = 'sites',label=traitTREM$id.TreM.type,cex=0.4)
ordiellipse(nmds, traitTREM$cluster17, display = "sites", kind = "sd", label = T,cex=0.8)
ordispider(nmds, traitTREM$cluster17,label=F)

#valeurs de dissmilarites par strate ?
dist.saproclim<-vegdist(speresPB[,5:30],distance = "bray")
mean(dist.saproclim)

############################################### dissimilarité intra-groupes de TreMs
TREM<-vegdist(traitTREM[,4:42],distance = "bray")
piege.ano<-anosim(TREM,traitTREM$id.TreM.group)
summary(piege.ano)
dissim_matrix <- as.matrix(piege.ano$dis.rank)
dissim_matrix
piege.ano$class.vec
write.table(piege.ano$dis.rank,file="",sep=";")
write.table(piege.ano$class.vec,file="",sep=";")

# permanova
#npmanova contraint par cldia
adonis(speres[,5:30]~speres$ess,speres, permutations = 999, method = "bray",strata = speres$cldia)


#ANOSIM
dist.saproclim<-vegdist(spefeu[,5:36],distance = "bray")
piege.ano<-anosim(dist.saproclim,spefeu$ess)
summary(piege.ano)

# identification des especes qui contribuent a la dissimilarité inter-categories
############################################### Simper
library(vegan)
traitTREM3<-read.table("P:/Mes documents/bois mort/Larrieu/ADD typo/table_traits_treMs3.txt",header=T)
simperTREM<-simper(traitTREM3[,4:42], traitTREM3$cluster, permutations = 100)
toto<-summary(simperTREM)

# recherche d'espèces caractéristiques de catégories (e.g. strates)
#Indval
library(indicspecies)

# Restreindre analyse aux esp?ces pr?sentes dans plus de 10% des sites
#sp.n.sites <- apply(coleoabond[,13:319], 2, function(x) sum(x > 0))
#sp.common <- sp.n.sites > 0.10*nrow(coleoabond[,13:319])
#coleo_common <- coleoabond[,13:319][, sp.common]

# Restreindre analyse aux esp?ces repr?sent?es par + de 10 individus
#sp.n.ind <- apply(coleo_common, 2, function(x) sum(x))
#sp.abund <- sp.n.ind > 10
#coleo_abund <- coleo_common[, sp.abund]

indvalfeu <- multipatt(spefeu[,5:36], spefeu$ess, func="IndVal.g", duleg=T, control=how(nperm=9999))
summary(indvalfeu, indvalcomp=TRUE)

############################################### diagramme de Venn inter-strates

library(ggvenn)

ggvenn(
  frisbee[,2:5], 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggvenn(frisbee, c("Cav", "Fon", "Lig", "Suc"), show_percentage = TRUE)

#partition de la variance ou de l'inertie expliquée par chaque prédicteur ou groupe de prédicteurs

#partition de variance sur RS par piege
# 4 groupes de variables
#1.	Micro-habitat density = dcav+dfun+dsap+dperdw
#2.	Large DW density = LSDW+LLDW
#3.	Stand conditions = open
#4.	Stand type = mixture+forest.type
library(vegan)
partrs<-varpart(tabtrap2$rs,~dcav+dfun+dsap+dperdw,~LSDW+LLDW,~open,~mixture+forest.type,data=tabtrap2)
partrs

#test de l'effet propre+conjoint de chaque groupe de variable
rda.result<-rda(tabtrap2D$rs~dcav+dfun+dsap+dperdw,data=tabtrap2D)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(subtrapD2$rs~dcav+dfun+dsap+dperdw+Condition(LSDW)+Condition(LLDW)+Condition(open)+Condition(forest.type),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(subtrapD2$rs~LSDW+LLDW,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(subtrapD2$rs~LSDW+LLDW+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(open)+Condition(forest.type),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(subtrapD2$rs~open,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(subtrapD2$rs~open+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(LSDW)+Condition(LLDW)+Condition(forest.type),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(subtrapD2$rs~forest.type,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(subtrapD2$rs~forest.type+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(LSDW)+Condition(LLDW)+Condition(open),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

########################partition d'inertie sur composition
partcompo<-varpart(toto3D,~dcav+dfun+dsap+dperdw,~LSDW+LLDW,~open,~forest.type,data=subtrapD2)
partcompo

#test de l'effet propre+conjoint de chaque groupe de variable
rda.result<-rda(toto3D~dcav+dfun+dsap+dperdw,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(toto3D~dcav+dfun+dsap+dperdw+Condition(LSDW)+Condition(LLDW)+Condition(open)+Condition(forest.type),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(toto3D~LSDW+LLDW,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(toto3D~LSDW+LLDW+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(open)+ Condition(forest.type),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(toto3D~open,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(toto3D~open+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(LSDW)+Condition(LLDW)+ Condition(forest.type))
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

rda.result<-rda(toto3D~forest.type,data=subtrapD2)
anova(rda.result, step=200, perm.max=200)
rda.result<-rda(toto3D~forest.type+Condition(dsap)+Condition(dperdw)+Condition(dcav)+Condition(dfun)+Condition(LSDW)+Condition(LLDW)+Condition(open),data=subtrapD2)
RsquareAdj(rda.result)
anova(rda.result, step=200, perm.max=200)

# additive diversity partitioning
#evaluate the contribution of α- and β-diversity to the γ-diversity of canopy-dwelling Hymenoptera over the entire sampled area (function adipart, index = ‘richness’, 1000 permutations, vegan R-package).
# At the plot scale, we used three levels of differentiation:  plot, plot decline category and the whole sampled area. α plot corresponds to the average species diversity per plot, β plot to the diversity
# among plots and β cat_plot to the diversity among plot decline categories.γ diversity is the sum of α plot, β plot and β cat_plot (γ = α plot + β plot + β cat_plot). Consequently, four levels of differentiation
# were included: plot, stand, stand decline category and the whole sampled
#area. At the stand scale, four levels of differentiation were
#included: plot, stand, stand decline category and the whole sampled
#area. As an additional measure of β-diversity, we calculated the dissimilarity
#between plots, stands, decline categories and pairs of
#decline categories to separate β-diversity (Sorensen dissimilarity index) into β-turnover (Simpson dissimilarity index) and β-nestedness
# (nestedness-resultant fraction of Sorensen dissimilarity) (function beta.multi, Sorensen family, betapart R-package;. 
# β-turnover indicates the replacement of some species by others, whereas β-nestedness indicates that species assemblages are subsets of species occurring at larger spatial scales.

# test de comparaison multiple entre >2 catégories en cas d'effet significatif d'une variable catégorielle
library(multcomp)
RES<-glht(glm1,linfct=mcp(forest.type2="Tukey"))
summary(RES)

# identifier le meilleur modèle multi-variables expliquant une variable réponse
library(arm)
library(MuMIn)

model1=lmer(ab.myco~LSDW+LLDW+LT+dcav+dsap+dfun+dperdw+open+(1|site),family=poisson,data=tabtrap2T)
stdz.model1=standardize(model1,standardize.y=F)
dd1<-dredge(stdz.model1,m.max=2)
model2<-model.avg(get.models(dd1,subset=delta<100))
model2
summary(model2)

# test de comparaison multiple entre >2 catégories en cas d'effet significatif d'une variable catégorielle
library(multcomp)
RES<-glht(glm1,linfct=mcp(forest.type2="Tukey"))
summary(RES)

# identifier le meilleur modèle multi-variables expliquant une variable réponse
library(arm)
library(MuMIn)

model1=lmer(ab.myco~LSDW+LLDW+LT+dcav+dsap+dfun+dperdw+open+(1|site),family=poisson,data=tabtrap2T)
stdz.model1=standardize(model1,standardize.y=F)
dd1<-dredge(stdz.model1,m.max=2)
model2<-model.avg(get.models(dd1,subset=delta<100))
model2
summary(model2)