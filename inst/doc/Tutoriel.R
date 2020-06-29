## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
library(knitr)
library(rmdformats)

oldpar <- par() 
oldoptions <- options()

## Global options
options(max.print="75")
opts_chunk$set(echo=TRUE,
	             cache=FALSE,
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE)
opts_knit$set(width=75)

## ---- message=FALSE-----------------------------------------------------------
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(WeightedCluster)
library(FactoMineR)
library(ade4)
library(RColorBrewer)
library(questionr)
library(GDAtools)
library(dplyr)
library(purrr)
library(ggplot2)
library(seqhandbook)

## -----------------------------------------------------------------------------
# chargement des trajectoires
data(trajact)
str(trajact)

## -----------------------------------------------------------------------------
# définition du corpus de séquences
labs <- c("études","temps plein","temps partiel","petits boulots","inactivité","serv. militaire")
palette <- brewer.pal(length(labs), 'Set2')
seqact <- seqdef(trajact, labels=labs, cpal=palette)

## -----------------------------------------------------------------------------
# nombre de séquences distinctes
seqtab(seqact, idx=0) %>% nrow

## ----fig1, fig.align="center", out.width="80%"--------------------------------
# chronogramme
seqdplot(seqact, xtlab=14:50, cex.legend=0.7)

## -----------------------------------------------------------------------------
# chargement des variables socio-démographiques
data(socdem)
str(socdem)

## -----------------------------------------------------------------------------
indics <- seqinepi(seqact)
head(indics)

## -----------------------------------------------------------------------------
# matrice de distance à partir des indicateurs
dissim <- dist(indics, method='euclidean') %>% as.matrix

# matrice de distance à partir des résultats d'une ACP
acp_coords <- PCA(indics, scale.unit=FALSE, ncp=5, graph=FALSE)$ind$coord
dissim <- dist(acp_coords, method='euclidean') %>% as.matrix

## -----------------------------------------------------------------------------
# codage disjonctif complet
disjo <- dichotom(seqact)
disjo <- disjo[,colSums(disjo)>0]

# distance euclidienne
dissim <- dist(disjo, method='euclidean') %>% as.matrix

# distance du chi2
dissim <- map_df(disjo, as.factor) %>%
          dudi.acm(scannf=FALSE, nf=ncol(disjo)) %>%
          dist.dudi() %>%
          as.matrix

# après une ACP
acp_coords <- PCA(disjo, scale.unit=FALSE, ncp=5, graph=FALSE)$ind$coord
dissim <- dist(acp_coords, method='euclidean') %>% as.matrix

# après une ACM
acm_res <- purrr::map_df(disjo, as.factor) %>%
           MCA(ncp=5, graph=FALSE)
dissim <- dist(acm_res$ind$coord, method='euclidean') %>% as.matrix

## -----------------------------------------------------------------------------
# codage AHQ
ahq <- seq2qha(seqact, c(1,3,7,10,15,20,28))
ahq <- ahq[,colSums(ahq)>0]

# distance du chi2
dissim <- dudi.coa(ahq, scannf=FALSE, nf=ncol(ahq)) %>%
          dist.dudi() %>%
          as.matrix

# après une AFC
afc_coord <- CA(ahq, ncp=5, graph=FALSE)$row$coord
dissim <- dist(afc_coord, method='euclidean') %>% as.matrix

## ---- message=FALSE-----------------------------------------------------------
# construction de la matrice de distance
couts <- seqsubm(seqact, method="CONSTANT", cval=2)
dissim <- seqdist(seqact, method="OM", sm=couts, indel=1.5)

## ---- eval=FALSE--------------------------------------------------------------
#  # sequencing
#  dissim <- seqdist(seqact, method="OMstran", otto=0.1, sm=couts, indel=1)
#  dissim <- seqdist(seqact, method="OMspell", expcost=0, sm=couts, indel=1)
#  dissim <- seqdist(seqact, method="SVRspell", tpow=0)
#  
#  # timing
#  dissim <- seqdist(seqact, method="HAM", sm=couts)
#  dissim <- seqdist(seqact, method="CHI2", step=1)
#  
#  # duration (distribution aver the entire period)
#  dissim <- seqdist(seqact, method="EUCLID", step=37)
#  
#  # duration (spell lengths)
#  dissim <- seqdist(seqact, method="OMspell", expcost=1, sm=couts, indel=1)
#  dissim <- seqdist(seqact, method="LCS")

## -----------------------------------------------------------------------------
# classification ascendante hiérarchique
agnes <- as.dist(dissim) %>% agnes(method="ward", keep.diss=FALSE)

## ---- fig.align="center", out.width="80%"-------------------------------------
# dendrogramme
as.dendrogram(agnes) %>% plot(leaflab="none")

## ---- fig.align="center", out.width="80%"-------------------------------------
# heatmap (dendrogramme + index plot)
seq_heatmap(seqact, agnes)

## ---- fig.align="center", out.width="80%"-------------------------------------
# graphique des sauts d'inertie
plot(sort(agnes$height, decreasing=TRUE)[1:20], type="s", xlab="nombre de classes", ylab="inertie")

## -----------------------------------------------------------------------------
# indicateurs de qualité
wardRange <- as.clustrange(agnes, diss=dissim)
summary(wardRange, max.rank=2)

## ---- fig.align="center", out.width="80%"-------------------------------------
plot(wardRange, stat=c('ASW','R2','CH'), norm="zscore")

## -----------------------------------------------------------------------------
# choix de la partition en 5 classes
nbcl <- 5
part <- cutree(agnes, nbcl)

## -----------------------------------------------------------------------------
# consolidation de la partition
newpart <- wcKMedoids(dissim, k=nbcl, initialclust=part, cluster.only=TRUE)
table(part, newpart)
wcClusterQuality(dissim, part)$stats['R2sq'] %>% round(3)
wcClusterQuality(dissim, newpart)$stats['R2sq'] %>% round(3)

## -----------------------------------------------------------------------------
part <- as.numeric(as.factor(newpart))

## -----------------------------------------------------------------------------
# classification floue (fuzzy clustering)
fanny <- as.dist(dissim) %>% fanny(k=5, metric='euclidean', memb.exp=1.2)
fanny$membership %>% round(2) %>% .[1:3,]

## ---- fig.align="center", out.width="80%"-------------------------------------
# chronogrammes de la typologie
seqdplot(seqact, group=part, xtlab=14:50, border=NA, cex.legend=0.8)

## ---- fig.align="center", out.width="80%"-------------------------------------
# index plots de la typologie
seqIplot(seqact, group=part, xtlab=14:50, yaxis=FALSE, cex.legend=0.8)

## ---- fig.align="center", out.width="80%"-------------------------------------
# index plots de la typologie, triés par multidimensional scaling
mds.order <- cmdscale(dissim,k=1)
seqIplot(seqact, sortv=mds.order, group=part, xtlab=14:50, yaxis=FALSE, cex.legend=0.8)

## -----------------------------------------------------------------------------
smoothed <- seqsmooth(seqact, dissim, k=30)$seqdata
seqIplot(smoothed, sortv=mds.order, group=part, xtlab=14:50, yaxis=FALSE, cex.legend=0.8)

## ---- fig.align="center", out.width="80%", message=FALSE----------------------
# relative frequency sequence plots
seqplot.rf(seqact, diss=dissim, group=part, xtlab=14:50)

## ---- fig.align="center", out.width="80%"-------------------------------------
# frequency plots
seqfplot(seqact, group=part, ylab="", xtlab=14:50, cex.legend=0.8)

## ---- fig.align="center", out.width="80%"-------------------------------------
# modal state plots
seqmsplot(seqact, group=part, xtlab=14:50, cex.legend=0.8)

## ---- fig.align="center", out.width="80%", message=FALSE, include=FALSE-------
# representative sequence plots
seqrplot(seqact, group=part, diss=dissim, nrep=10, xtlab=14:50)

## -----------------------------------------------------------------------------
# effectifs
table(part)

## -----------------------------------------------------------------------------
# pourcentages
100*(table(part))/length(part)

## -----------------------------------------------------------------------------
# distances intra-classes
Dintra <- integer(length=nbcl)
for(i in 1:nbcl) Dintra[i] <- round(mean(dissim[part==i,part==i]),1)
Dintra

## -----------------------------------------------------------------------------
# distances moyennes au centre de la classe
dissassoc(dissim, part)$groups

## -----------------------------------------------------------------------------
# entropie transversale moyenne par classe
entropie <- vector()
for(i in 1:nbcl) entropie[i] <- round(mean(seqstatd(seqact[part==i,])$Entropy),2)
entropie

## ---- message=FALSE-----------------------------------------------------------
# durées dans les états selon la classe
dur <- seqistatd(seqact)
durees <- round(aggregate(dur, by=list(part), FUN=mean), 1)
rownames(durees) <- NULL
durees

## ---- message=FALSE-----------------------------------------------------------
# au moins un épisode dans les états, selon la classe
epi <- seqi1epi(seqact)
episodes <- round(aggregate(epi, by=list(part), FUN=mean), 2)
rownames(episodes) <- NULL
episodes

## -----------------------------------------------------------------------------
assoc.twocat(factor(part), socdem$sexe)

## -----------------------------------------------------------------------------
catdesc(factor(part), socdem, min.phi=0.1)

## -----------------------------------------------------------------------------
# parangon de chaque classe (numéros de ligne dans le fichier de données)
medoids(dissim, part)

## ---- message=FALSE-----------------------------------------------------------
ref <- seqdef(as.matrix("(1,4)-(2,33)"), informat="SPS", alphabet=alphabet(seqact))
distref <- seqdist(seqact, refseq = ref, method="OM", sm=couts, indel=1.5)

## ---- message=FALSE-----------------------------------------------------------
socdem %>% select(sexe,nbenf) %>%
           tibble(distref=distref) %>%
           ggplot(aes(x=nbenf, y=distref)) + 
             geom_boxplot(aes(fill=sexe), notch=T) +
             xlab("nombre d'enfants") +
             ylab("distance à la référence") +
             theme_bw()

## -----------------------------------------------------------------------------
# distances intra-classes selon le sexe
sapply(levels(socdem$sexe), function(x) round(mean(dissim[socdem$sexe==x,socdem$sexe==x]),1))

## ---- fig.align="center", out.width="80%"-------------------------------------
mds <- cmdscale(dissim, k=2)
plot(mds, type='n', xlab="axe 1", ylab="axe 2")
abline(h=0, v=0, lty=2, col='lightgray')
points(mds, pch=20, col=part)
legend('topleft', paste('classe',1:nbcl), pch=20, col=1:nbcl, cex=0.8)
text(aggregate(mds, list(socdem$sexe), mean)[,-1], levels(socdem$sexe), col='orange', cex=1, font=2)

## ---- message=FALSE-----------------------------------------------------------
# durées dans les états selon le sexe
dur <- seqistatd(seqact)
durees_sexe <- aggregate(dur, by=list(socdem$sexe), function(x) round(mean(x),1))
rownames(durees_sexe) <- NULL
colnames(durees_sexe) <- c("classe",labs)
durees_sexe

## -----------------------------------------------------------------------------
# turbulence
turbu <- aggregate(seqST(seqact), list(socdem$annais), mean)
plot(turbu, type='l', ylim=c(0,10), xlab='Année de naissance')

## -----------------------------------------------------------------------------
# analyse de variance avec le sexe comme facteur
dissassoc(dissim, socdem$sexe)

## ---- fig.align="center", out.width="80%"-------------------------------------
# analyse de variance selon la position dans le temps
diff <- seqdiff(seqact, group=socdem$sexe)
rownames(diff$stat) <- rownames(diff$discrepancy) <- 14:49
plot(diff, stat="Pseudo R2")

## ---- fig.align="center", out.width="80%"-------------------------------------
pal <- brewer.pal(ncol(diff$discrepancy), "Set2")
plot(diff, stat="discrepancy", legend.pos=NA, col=pal, lwd=1.5)
legend('topright', fill=pal, legend=colnames(diff$discrepancy), cex=0.7)

## -----------------------------------------------------------------------------
# analyse de la variance avec facteurs multiples
dissmfacw(dissim ~ annais+nbenf+sexe+diplome, data=socdem)

## -----------------------------------------------------------------------------
# arbre d'induction
arbre <- seqtree(seqact ~ annais+nbenf+sexe+diplome, data=socdem, diss=dissim, min.size=0.1, max.depth=3)

## -----------------------------------------------------------------------------
# résultats de l'arbre sous forme textuelle
print(arbre)

## ---- eval=FALSE--------------------------------------------------------------
#  # résultats de l'arbre sous forme graphique
#  seqtreedisplay(arbre,type="d",border=NA,show.depth=TRUE)

## ---- fig.align="center", out.width = '100%', echo=FALSE----------------------
knitr::include_graphics("http://nicolas.robette.free.fr/Docs/seqtree.png")

## ---- fig.align="center", out.width="80%"-------------------------------------
# statistiques implicatives
implic <- seqimplic(seqact, group=socdem$sexe)
# par(mar=c(2,2,2,2))
plot(implic, xtlab=14:50, lwd=2, conf.level=c(0.95, 0.99), cex.legend=0.7)

## -----------------------------------------------------------------------------
data(seqmsa)
trajmat <- seqmsa %>% select(starts_with('smat'))
str(trajmat)

## -----------------------------------------------------------------------------
trajenf <- seqmsa %>% select(starts_with('nenf'))
str(trajenf)

## -----------------------------------------------------------------------------
trajlog <- seqmsa %>% select(starts_with('slog'))
str(trajlog)

## -----------------------------------------------------------------------------
# définition de la trajectoire matrimoniale
labs <- c("jamais","union libre","marié","separé")
palette <- brewer.pal(length(labs), 'Set2')
seqmat <- seqdef(trajmat, labels=labs, cpal=palette)

## -----------------------------------------------------------------------------
# définition de la trajectoire parentale
labs <- c("0","1","2","3+")
palette <- brewer.pal(length(labs), 'YlOrRd')
seqenf <- seqdef(trajenf, labels=labs, cpal=palette)

## -----------------------------------------------------------------------------
# définition de la trajectoire d'indépendance résidentielle
labs <- c("non indépendant","indépendant")
palette <- brewer.pal(3, 'Set1')[1:2]
seqlog <- seqdef(trajlog, labels=labs, cpal=palette)

## ---- message=FALSE-----------------------------------------------------------
# matrices de distances des différentes dimensions
dmat <- seqdist(seqmat, method="OM", indel=1.5, sm=seqsubm(seqmat,"CONSTANT",cval=2))
denf <- seqdist(seqenf, method="OM", indel=1.5, sm=seqsubm(seqenf,"CONSTANT",cval=2))
dlog <- seqdist(seqlog, method="OM", indel=1.5, sm=seqsubm(seqlog,"CONSTANT",cval=2))

## ---- message=FALSE-----------------------------------------------------------
# matrice de distances entre séquences multidimensionnelles
dissim.MSA <- seqdistmc(list(seqmat,seqenf,seqlog), method="OM", indel=1.5, sm=as.list(rep("CONSTANT",3)), cval=2)

## -----------------------------------------------------------------------------
asso <- assoc.domains(list(dmat,denf,dlog), c('mat','enf','log'), dissim.MSA)
asso

## ---- fig.align="center", out.width="80%"-------------------------------------
# ACP à partir des corrélations entre dimensions
matcor <- asso$correlations$pearson[1:3,1:3]
PCA <- PCA(matcor, scale.unit=F, graph=F)
plot.PCA(PCA, choix='varcor')

## ---- message=FALSE-----------------------------------------------------------
dnolog <- seqdistmc(list(seqmat,seqenf), method="OM", indel=1.5, sm=as.list(rep("CONSTANT",2)), cval=2)
dnomat <- seqdistmc(list(seqenf,seqlog), method="OM", indel=1.5, sm=as.list(rep("CONSTANT",2)), cval=2)
dnoenf <- seqdistmc(list(seqmat,seqlog), method="OM", indel=1.5, sm=as.list(rep("CONSTANT",2)), cval=2)
dnolog <- as.numeric(as.dist(dnolog))
dnomat <- as.numeric(as.dist(dnomat))
dnoenf <- as.numeric(as.dist(dnoenf))

## -----------------------------------------------------------------------------
dmat %>% as.dist %>% as.numeric %>% cor(dnomat) %>% round(3)
denf %>% as.dist %>% as.numeric %>% cor(dnoenf) %>% round(3)
dlog %>% as.dist %>% as.numeric %>% cor(dnolog) %>% round(3)

## ---- fig.align="center", out.width="80%"-------------------------------------
mds.msa <- cmdscale(dmat,k=1)
par(mfrow=c(3,2), mar=c(2.1,2.1,2.1,2.1))
seqIplot(seqmat, sortv=mds.msa, xtlab=14:35, with.legend=FALSE, yaxis=FALSE, ylab="")
seqlegend(seqmat, cex=0.7)
seqIplot(seqenf, sortv=mds.msa, xtlab=14:35, with.legend=FALSE, yaxis=FALSE, ylab="")
seqlegend(seqenf, cex=0.7)
seqIplot(seqlog, sortv=mds.msa, xtlab=14:35, with.legend=FALSE, yaxis=FALSE, ylab="")
seqlegend(seqlog, cex=0.7)

## ---- fig.align="center", out.width="80%"-------------------------------------
# classification ascendante hiérarchique
agnes.MSA <- agnes(as.dist(dissim.MSA), method="ward", keep.diss=FALSE)
plot(as.dendrogram(agnes.MSA), leaflab="none")

## -----------------------------------------------------------------------------
# choix d'une typologie en 5 classes
nbcl.MSA <- 5
part.MSA <- cutree(agnes.MSA, nbcl.MSA) %>% factor

## ---- fig.align="center", out.width = '100%', eval=FALSE----------------------
#  # chronogrammes de la typologie
#  # ngroups <- nlevels(part.MSA)
#  par(mfrow=c(3,nbcl.MSA+1), mar=c(2.5, 2.1, 2.1, 2.1))
#  for(i in 1:nbcl.MSA) seqdplot(seqmat[part.MSA==i,], xtlab=14:35, border=NA, with.legend=FALSE, main=paste('classe',i))
#  seqlegend(seqmat, cex=0.5)
#  for(i in 1:nbcl.MSA) seqdplot(seqenf[part.MSA==i,], xtlab=14:35, border=NA, with.legend=FALSE)
#  seqlegend(seqenf, cex=0.5)
#  for(i in 1:nbcl.MSA) seqdplot(seqlog[part.MSA==i,], xtlab=14:35, border=NA, with.legend=FALSE)
#  seqlegend(seqlog, cex=0.5)

## -----------------------------------------------------------------------------
# chargement des données
data(seqgimsa)
trajfilles <- seqgimsa %>% select(starts_with('f'))
str(trajfilles)

## -----------------------------------------------------------------------------
trajmeres <- seqgimsa %>% select(starts_with('m'))
str(trajmeres)

## -----------------------------------------------------------------------------
# définition des séquences
lab.meres <- c("indép","moyen/sup","popu","inactivité","études")
pal.meres <- brewer.pal(5, "Set1")
seqmeres <- seqdef(trajmeres,lab=lab.meres, cpal=pal.meres)

## -----------------------------------------------------------------------------
lab.filles <- c("études","inactivité","temps partiel","temps plein")
pal.meres <- brewer.pal(4, "Set1")
seqfilles <- seqdef(trajfilles,lab=lab.filles, cpla=pal.filles)

## ---- message=FALSE-----------------------------------------------------------
# étape 1 : mesure de dissimilarité
dmeres <- seqdist(seqmeres,method="LCS")
cout.filles <- seqsubm(seqfilles, method="CONSTANT", cval=2)
dfilles <- seqdist(seqfilles, method="HAM", sm=cout.filles)

## -----------------------------------------------------------------------------
# étape 2 : multidimensional scaling
mds.meres <- cmdscale(dmeres, k=20, eig=TRUE)
mds.filles <- cmdscale(dfilles, k=20, eig=TRUE)

## ---- fig.align="center", out.width = '100%'----------------------------------
# choix du nombre de dimensions à retenir pour les mères
par(mfrow=c(1,2))
# mesure de stress
seqmds.stress(dmeres, mds.meres) %>% plot(type='l', xlab='nombre de facteurs', ylab='stress')
# part de variance expliquée
(mds.meres$eig[1:10]/mds.meres$eig[1]) %>% plot(type='s', xlab='nombre de facteurs', ylab='part de variance expliquée')

## ---- fig.align="center", out.width = '100%'----------------------------------
# choix du nombre de dimensions à retenir pour les filles
par(mfrow=c(1,2))
# mesure de stress
seqmds.stress(dfilles, mds.filles) %>% plot(type='l', xlab='nombre de facteurs', ylab='stress')
# part de variance expliquée
(mds.filles$eig[1:10]/mds.filles$eig[1]) %>% plot(type='s', xlab='nombre de facteurs', ylab='part de variance expliquée')

## -----------------------------------------------------------------------------
# étape 3 : PLS symétrique
a <- mds.meres$points[,1:5]
b <- mds.filles$points[,1:4]
pls <- symPLS(a,b)

## ----  eval=FALSE-------------------------------------------------------------
#  # étape 4 : distance et classification
#  
#  # pas de pondération
#  F <- pls$F
#  G <- pls$G
#  
#  # pondération par la variance des composantes de la PLS (w1)
#  F <- apply(pls$F,2,scale,center=FALSE)
#  G <- apply(pls$G,2,scale,center=FALSE)
#  
#  # pondération par le nombre de séquences distinctes (w2)
#  F <- pls$F/nrow(seqtab(seqmeres,tlim=0))
#  G <- pls$G/nrow(seqtab(seqfilles,tlim=0))
#  
#  # pondération par la 1ère valeur propre du MDS (w3)
#  F <- pls$F/mds.meres$eig[1]
#  G <- pls$G/mds.filles$eig[1]

## -----------------------------------------------------------------------------
# pondération par la variance des composantes de la PLS (w1)
F <- apply(pls$F,2,scale,center=FALSE)
G <- apply(pls$G,2,scale,center=FALSE)

## -----------------------------------------------------------------------------
# calcul de distance
diff2 <- function(X) return(as.matrix(dist(X,upper=T,diag=T)^2,nrow=nrow(X)))
D <- (diff2(F)+diff2(G))^0.5

## ---- fig.align="center", out.width = '80%'-----------------------------------
# classification
dist.GIMSA <- as.dist(D)
agnes.GIMSA <- agnes(dist.GIMSA, method="ward", keep.diss=FALSE)
plot(as.dendrogram(agnes.GIMSA), leaflab="none")

## -----------------------------------------------------------------------------
# partition en 5 classes
nbcl.GIMSA <- 5
part.GIMSA <- cutree(agnes.GIMSA, nbcl.GIMSA) %>% factor

## ---- fig.align="center", out.width = '100%'----------------------------------
par(mfrow=c(3,nbcl.GIMSA), mar=c(2.5, 2.1, 2.1, 2.1))
for(i in 1:nbcl.GIMSA) seqdplot(seqmeres[part.GIMSA==i,], xtlab=14:60, border=NA, with.legend=FALSE, main=paste('classe',i))
for(i in 1:nbcl.GIMSA) seqdplot(seqfilles[part.GIMSA==i,], xtlab=1:15, border=NA, with.legend=FALSE)
seqlegend(seqmeres, cex=0.6)
seqlegend(seqfilles, cex=0.6)

## ---- echo=FALSE--------------------------------------------------------------
par(oldpar)
options(oldoptions)

