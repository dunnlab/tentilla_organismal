#R code used for the analyses and figures in Damian-Serrano et al. Zoologica Scripta
library(tidyverse)
library(stringr)
library(reshape2)
## Biological
library(Rphylip)
library(arbutus)
library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(picante)
library(geiger)
library(phylobase)
library(fields)
library(phylosignal)
library(geomorph)
library(phylopath)
library(phylolm)
library(convevol)
library(surface)
## Graphics
library(FactoMineR)
library(factoextra)
library(corrplot)
library(BAMMtools)
library(gridExtra)
library(xtable)
library(colorRamps)
library(ggfortify)
library(ggpubr)

# Set paths to input data
setwd("~/tentilla_morph/Organismal_paper/")

#Load raw data
read.csv("~/tentilla_morph/Supplementary_materials/Dryad/raw_morphology_data.csv") -> numbers
numbers$Species = as.character(numbers$Species)
categorical <- read.csv("~/tentilla_morph/Organismal_paper/IOB-PostReview/raw_figures/CategoricalCharacters.csv")[,-2]
rownames(categorical) = categorical$Species

#Correct species spellings
numbers$Species[which(numbers$Species == "Agalma okeni")] <- "Agalma okenii"
numbers$Species[which(numbers$Species == "Forskalia edwardsi")] <- "Forskalia edwardsii"
numbers$Species[which(numbers$Species == "Rhizophysa eysenhardti")] <- "Rhizophysa eysenhardtii"

#Merge Nanomias
numbers$Species[which(numbers$Species == "Nanomia cara")] <- "Nanomia sp"
numbers$Species[which(numbers$Species == "Nanomia bijuga")] <- "Nanomia sp"

#Polish fields
castnumbers = numbers
castnumbers[castnumbers=="needDIC"] <- NA
castnumbers[castnumbers=="needConfocal"] <- NA
castnumbers[castnumbers=="needMature"] <- NA
castnumbers[castnumbers=="needReslide"] <- NA
castnumbers[castnumbers==-1] <- NA
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.character)
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.numeric)

#Define functions to deal with NAs and decimals
mean.na <- function(a){mean(a, na.rm = TRUE)}
var.na <- function(a){var(a, na.rm = TRUE)}
round3 <- function(a){round(a, 3)}

#Compute morphometric ratios
coiledness = castnumbers$Cnidoband.free.length..um./castnumbers$Cnidoband.length..um.
heteroneme_elongation = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.width..um.
haploneme_elongation = castnumbers$Haploneme.free.length..um./castnumbers$Haploneme.width..um.
desmoneme_elongation = castnumbers$Desmoneme.length..um./castnumbers$Desmoneme.width..um.
rhopaloneme_elongation = castnumbers$Rhopaloneme.length..um./castnumbers$Rhopaloneme.width..um.
heteroneme_shaft_extension = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.shaft.free.length..um.
Heteroneme_to_CB = castnumbers$Heteroneme.free.length..um./(castnumbers$Cnidoband.free.length..um.+0.001)
total_heteroneme_volume = castnumbers$Heteroneme.volume..um3.*castnumbers$Heteroneme.number
total_haploneme_volume = (castnumbers$Cnidoband.free.length..um./castnumbers$Haploneme.width..um.)*((4*pi/3)*0.5*castnumbers$Haploneme.free.length..um.*((0.5*castnumbers$Haploneme.width..um.)^2))
cnidomic_index = log(total_heteroneme_volume + total_haploneme_volume)
cnidomic_index[which(is.na(cnidomic_index))] <- log(total_haploneme_volume[which(is.na(cnidomic_index))])
cnidomic_index[which(is.na(cnidomic_index))] <- log(castnumbers[which(is.na(cnidomic_index)), "Heteroneme.volume..um3."])
cnidomic_index[cnidomic_index==-Inf]<-NA

#Define surface and volume functions
surface_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  S <- 4*pi*((((a*b)^1.6)+((a*c)^1.6)+((c*b)^1.6))/3)^(1/1.6)
  return(S)
}
volume_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  V <- (4/3)*pi*a*b*c
  return(V)
}

SAV_haploneme = surface_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)/volume_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)
names(SAV_haploneme) = castnumbers$Species
SAV_heteroneme = surface_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)/volume_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)
names(SAV_heteroneme) = castnumbers$Species

#Combine all morphometric variables
morphometrics = data.frame(castnumbers$slide_id, castnumbers$Species,coiledness, heteroneme_elongation, haploneme_elongation, desmoneme_elongation, rhopaloneme_elongation,heteroneme_shaft_extension, Heteroneme_to_CB, total_heteroneme_volume, total_haploneme_volume, cnidomic_index, SAV_haploneme)
names(morphometrics)[1:2] = c("slide_id"," Species")

#Combine morphometrics with regular characters
castnumbers = data.frame(castnumbers,morphometrics[,c(-1,-2)])

#Logtransform variables which are not normal
castlogs = castnumbers
non_normal = apply(castnumbers[,c(-1,-2,-3)], 2, shapiro.test) %>% lapply(function(x){x$p.value}) %>% unlist() 
non_normal = which(as.vector(non_normal) < 0.05)+3
castlogs[,non_normal] <- sapply(castnumbers[,non_normal], log)
castlogs[castlogs==-Inf]<-NA

#Means, varinces, standard errors
castmeans <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], mean.na, na.action = na.pass)
castmean_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], mean.na, na.action = na.pass)
castvariances <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], FUN = var.na, na.action = na.pass) #[,-2]
castvariance_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], FUN = var.na, na.action = na.pass)
castvariances[is.na(castvariances)] <- 0
castvariance_logs[is.na(castvariance_logs)] <- 0
castN <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], length, na.action = na.pass) 
castN = data.frame(castN$Species, rowMeans(castN[,-1])) #[,-2]
names(castN) = c("Species", "N_specimens")
castSE <- cbind(castvariances$Species, sqrt(castvariances[,-1])/sqrt(castN$N_specimens))
castSE_logs <- cbind(castvariance_logs$Species, sqrt(castvariance_logs[,-1])/sqrt(castN[,-1]))
names(castSE)[1] <- "Species"
names(castSE_logs)[1] <- "Species"

#Make overview of data and heatmap
Ses = castSE[,-1]
Ses[Ses == 0] = NA
morphdata = cbind(castN, castmeans[,-1], Ses)
morphdata = morphdata[,c(1:3,33,4,34,5,35,6,36,7,37,8,38,9,39,10,40,11,41,12,42,13,43,14,44,15,45,16,46,17,47,18,48,19,49,20,50,21,51,22,52,23,53,24,54,25,55,26,56,27,57,28,58,29,59,30,60,31,61,32,62)]
names(morphdata) = str_replace_all(names(morphdata),".1","_SE")
#write.csv(morphdata, "characterdata.csv")

heatdata = as.matrix(castmean_logs[,-1])
rownames(heatdata) = castmean_logs$Species
heatdata[is.nan(heatdata)]<- -1
#hcolors = grDevices::terrain.colors(20)
library(colorspace)
hcolors <- colorRampPalette(c("#e0ecf4","#9ebcda","#8856a7"), bias=1)
hcolors <- colorRampPalette(c("#ffeda0","#feb24c","#f03b20"), bias=1)
hcolors <- colorRampPalette(c("yellow","purple"), bias=1)
#hcolors <- sequential_hcl(17, h = 245, c = c(40, 75, 0), l = c(30, 95), power = 1)
hcolors <- c(rep("#000000FF",10), hcolors(30))
heatmap(heatdata, scale = "column", cexCol = 0.2, col=hcolors, keep.dendro = T)

#Load phylogenetic tree
consensus = read.nexus("~/tentilla_morph/Supplementary_materials/TreeBase/RB_constrained_timetree/TimeTree_siphs_mcmc_MAP.tre") %>% drop.tip(56:61)
consensus$tip.label = str_replace_all(consensus$tip.label,"_"," ")
#Switch Nanomia bijuga for Nanaomia sp
consensus$tip.label[which(consensus$tip.label == "Nanomia bijuga")] <- "Nanomia sp"

#Prune quant matrix to tree species
matrix = castnumbers[which(!is.na(castnumbers$Species)),]
matrix_logs = castlogs[which(!is.na(castlogs$Species)),]
sharedspp = matrix$Species[which(matrix$Species %in% consensus$tip.label)]
sharedmatrix = matrix[which(matrix$Species %in% sharedspp),]
sharedlogs = matrix_logs[which(matrix$Species %in% sharedspp),]
sharedmeans = castmeans[which(castmeans$Species %in% sharedspp),]
sharedmean_logs = castmean_logs[which(castmean_logs$Species %in% sharedspp),]
sharedvars = castvariances[which(castvariances$Species %in% sharedspp),]
sharedvar_logs = castvariance_logs[which(castvariance_logs$Species %in% sharedspp),]

#Prune categorical matrix to tree species
cat_rowNAs = apply(categorical[,-1], 1, function(x) sum(is.na(x)))
catmatrix = categorical[which(cat_rowNAs<1),] %>% .[,-1]
sharedspp_cat = rownames(catmatrix)[which(rownames(catmatrix) %in% consensus$tip.label)]
sharedcategorical = catmatrix[which(rownames(catmatrix) %in% sharedspp_cat),]

#Prune tree to quant matrix species
nodatatipnames = consensus$tip.label[which(!(consensus$tip.label %in% sharedmatrix$Species))]
nodatatips = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames)]
prunedtree = drop.tip(consensus, nodatatips)
prunedtree$tip.label
plot(prunedtree)
#ultram = chronos(prunedtree)
ultram = prunedtree
plot(ultram)

#Prune tree to categorical matrix species
nodatatipnames_cat = consensus$tip.label[which(!(consensus$tip.label %in% sharedspp_cat))]
nodatatips_cat = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames_cat)]
cat_tree = drop.tip(consensus, nodatatips_cat)
cat_tree$tip.label
plot(cat_tree)
#ultram_cat = chronos(cat_tree)
ultram_cat = cat_tree
plot(ultram_cat)
sharedcategorical = sharedcategorical[match(ultram_cat$tip.label, rownames(sharedcategorical)),]
cprunedmatrix = sharedcategorical[which(rownames(sharedcategorical)%in%rownames(sharedmatrix)),] %>% .[which(sapply(.,function(x) length(unique(x)))>1)]
sharedbinary = sharedcategorical
sharedbinary$Haploneme.type = as.character(sharedbinary$Haploneme.type)
sharedbinary[sharedbinary=="Isorhizas"] = 0
sharedbinary[sharedbinary=="Anisorhizas"] = 1
sharedbinary$Haploneme.type = as.numeric(sharedbinary$Haploneme.type)
sharedbinary$Heteroneme.type = as.character(sharedbinary$Heteroneme.type)
sharedbinary[sharedbinary=="Stenotele"] = 0
sharedbinary[sharedbinary=="Microbasic mastigophore"] = 1
sharedbinary[sharedbinary=="Eurytele"] = NA
sharedbinary$Heteroneme.type = as.numeric(sharedbinary$Heteroneme.type)

#Purely numerical character sets without slide or species
Q_sharedmatrix = sharedmatrix[,c(-1,-2,-3)]
rownames(sharedmeans) = sharedmeans$Species
Q_sharedmeans = sharedmeans[,-1]
rownames(sharedvars) = sharedvars$Species
Q_sharedvars = sharedvars[,-1]
Q_sharedmean_logs = sharedmean_logs[,-1]
rownames(Q_sharedmean_logs) = sharedmean_logs$Species

## Retrieve diet data ##
GC = read.csv("~/tentilla_morph/Supplementary_materials/Dryad/literature_diet_data.tsv", header = T, sep='\t')
GC$Prey.type = factor(GC$Prey.type, levels=unique(GC$Prey.type))
GC$Siphonophore.species = as.character(GC$Siphonophore.species)
#Fix typos#
GC$Siphonophore.species[which(GC$Siphonophore.species == "Nanomia bijuga")] <- "Nanomia sp"
GC$Siphonophore.species[which(GC$Siphonophore.species == "Rhizophysa eyesenhardti")] <- "Rhizophysa eysenhardtii"
GC$Siphonophore.species[which(GC$Siphonophore.species == "Agalma okeni")] <- "Agalma okenii"
GC = GC[,1:2]
names(GC)<-c("species", "character")
#write.csv(GC, "gutcontentliteraturereview.csv")
GC = split(GC,GC$character)
nrowGC = purrr::map(GC,nrow) %>% as.numeric()
GC = GC[which(nrowGC>0)]
GC = purrr::map(GC,unique)
diet= matrix(ncol=length(GC),nrow=length(unique(castlogs$Species))) %>% as.data.frame()
names(diet) = names(GC)
rownames(diet) = unique(castlogs$Species)
for(E in GC){
  print(E$species)
  for(S in E$species){
    diet[which(rownames(diet) == S),as.character(unique(E$character))] = 1
  }
}
diet[is.na(diet)] <- 0
diet = diet[which(rowSums(diet)>0),which(colSums(diet)<nrow(diet))]
diet=diet[,-c(3,7,14:17)]

#Add personal observations of the authors
cladeB <- matrix(rep(c(0,1,0,0,0,0,0,0,0,0,0),3),nrow=11, ncol = 3) %>% t() %>% as.data.frame()
rownames(cladeB) = c("Erenna richardi", "Erenna sirena", "Stephanomia amphytridis")
krilleaters = matrix(rep(c(0,0,0,1,0,0,1,0,0,0,0),4),nrow=11, ncol = 4) %>% t() %>% as.data.frame()
rownames(krilleaters) = c("Praya dubia", "Resomia ornicephala", "Lychnagalma utricularia", "Bargmannia amoena")
gelateaters = matrix(rep(c(0,0,0,0,0,0,0,0,0,1,0),1),nrow=11, ncol = 1) %>% t() %>% as.data.frame()
names(gelateaters) = names(diet)
rownames(gelateaters) = c("Apolemia rubriversa")
names(gelateaters)=names(diet)
names(krilleaters)=names(diet)
names(cladeB)=names(diet)
diet = rbind(diet, cladeB, krilleaters, gelateaters)
#Decapod diet column actually encompasses decapods, krill, and mysids (large crustaceans/shrimp like animals)

# Retrive ROV annotation data
VARS <- read.csv("~/tentilla_morph/Supplementary_materials/Dryad/raw_ROV_data.csv")
VARS_curated = VARS[which(VARS$Siphonophore.concept %in% ultram$tip.label | VARS$Siphonophore.concept=="Nanomia bijuga"),]
VARS_cast = acast(VARS_curated, Siphonophore.concept~Prey.taxonomy, fun.aggregate = length)

#Prune morphological matrix to species in diet
dprunedmatrix = sharedmeans[which(sharedmeans$Species%in%rownames(diet)),]
dprunedmatrix_logs = sharedmean_logs[which(sharedmean_logs$Species%in%rownames(diet)),]
#Prune tree to diet species
dprunedTree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(diet))))

hypdiet = c("Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Large crustacean", "Mixed", "Mixed", "Mixed", "Large crustacean", "Large crustacean", "Mixed", "Small crustacean", "Large crustacean", "Fish", "Fish", "Fish", "Large crustacean", "Gelatinous", "Fish", "Fish", "Fish")
names(hypdiet) = dprunedTree$tip.label

#Soft bodied vs hard bodied prey
soft_hard = cbind((diet$`Copepod diet`+diet$`Crustacean diet` + diet$`Decapod diet`+ diet$`Amphipod diet` + diet$`Ostracod diet`), (diet$`Fish diet`+diet$`Chaetognath diet`+diet$`Gelatinous diet`+diet$`Mollusc diet`+diet$`Polychaete diet`))
colnames(soft_hard) = c("Hard", "Soft")
rownames(soft_hard) = rownames(diet)
soft_hard <- as.data.frame(soft_hard)
soft_hard[soft_hard>0]<-1
softORhard = as.data.frame(soft_hard$Hard+soft_hard$Soft)
rownames(softORhard)=rownames(soft_hard)
names(softORhard) = "Type"
softORhard[softORhard==2]<-"Both"
softORhard[which(soft_hard$Hard==1 & soft_hard$Soft==0),]<-"Hard"
softORhard[which(soft_hard$Hard==0 & soft_hard$Soft==1),]<-"Soft"
softORhard

#PCA prep
Q_sharedmean_logs = sharedmean_logs[,-1]
#Q_sharedmean_logs = castmean_logs[,-1]
#rownames(Q_sharedmean_logs) = castmean_logs$Species
rownames(Q_sharedmean_logs) = sharedmean_logs$Species
raw_matrix = Q_sharedmean_logs[,c(1:2,4:20)] %>% .[which(!is.na(rowSums(.))),]
raw_matrix_notf = Q_sharedmean_logs[,c(1:2,4:20)] %>% .[,c(1:7,12:17)] %>% .[which(!is.na(rowSums(.))),]
raw_matrix_NaZeroes = Q_sharedmean_logs[,c(1:2,4:20)]
raw_matrix_NaZeroes[is.na(raw_matrix_NaZeroes)]<-0
compound_matrix = Q_sharedmean_logs[which(!is.na(rowSums(Q_sharedmean_logs))),c(3,21:31)]

temp <- castlogs[,c(-1,-3)]
temp[is.na(temp)] <- 0
fullraw_nozeroes <- aggregate(. ~ Species, data=temp, FUN = mean)
rownames(fullraw_nozeroes) <- fullraw_nozeroes$Species

#hypdiet but for all species not only those in tree
## ALSO reconsidering Forskalia interpretation as small crustacean specialist, and assigning all Forskalia species
hypdiet_full = c("Mixed","Mixed","Gelatinous","Gelatinous","Mixed","Large crustacean","Small crustacean","Small crustacean","Small crustacean","Small crustacean","Fish","Fish","Mixed","Small crustacean","Large crustacean","Small crustacean","Large crustacean","Fish","Large crustacean","Large crustacean","Large crustacean","Fish","Fish","Small crustacean","Small crustacean","Fish","Small crustacean","Small crustacean","Mixed", "Mixed", "Mixed")
names(hypdiet_full) = c("Agalma elegans", "Agalma okenii", "Apolemia rubriversa", "Apolemia uvaria", "Athorybia rosacea", "Bargmannia amoena", "Bassia bassensis", "Chelophyes appendiculata", "Cordagalma ordinatum", "Diphyes dispar", "Erenna richardi", "Erenna sirena", "Forskalia tholoides", "Hippopodius hippopus", "Lychnagalma utricularia", "Marrus orthocanna", "Nanomia sp", "Physalia physalis", "Praya dubia", "Resomia ornicephala", "Resomia persica", "Rhizophysa eysenhardtii", "Rhizophysa filiformis", "Rosacea cymbiformis", "Sphaeronectes koellikeri", "Stephanomia amphytridis", "Stephanophyes superba", "Sulculeolaria quadrivalvis", "Forskalia formosa", "Forskalia asymmetrica", "Forskalia edwardsii")
#dpruned_full = fullraw_nozeroes[which(fullraw_nozeroes$Species%in%names(hypdiet_full)),]
#dpruned_full <- data.frame(dpruned_full, "Diet" = hypdiet_full)
dpruned_full = data.frame(fullraw_nozeroes, "Diet" = hypdiet_full[match(rownames(fullraw_nozeroes), names(hypdiet_full))], stringsAsFactors = F)

### COMPARATIVE ANALYSES ###

phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmean_logs[,-1])))
#phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmeans[,-1])))
names(phylosignals) = c("K", "P", "Ntaxa")
for(i in 2:ncol(sharedmean_logs)){
  CH_I=as.numeric(sharedmean_logs[,i])
  #CH_I=as.numeric(sharedmeans[,i])
  names(CH_I) = sharedmean_logs$Species
  #names(CH_I) = sharedmeans$Species
  SE_I = as.numeric(castSE[,i]) %>% log()
  SE_I[which(SE_I==-Inf)] = 0
  names(SE_I) = castSE$Species
  CH_I= CH_I[!is.na(CH_I)]
  SE_I= SE_I[!is.na(SE_I)]
  SE_I = SE_I[which(names(SE_I)%in%names(CH_I))]
  CH_I = CH_I[which(names(CH_I)%in%names(SE_I))]
  treeI = drop.tip(ultram,which(!(ultram$tip.label %in% names(CH_I))))
  class(treeI) = "phylo"
  rownames(phylosignals)[i-1] <- names(sharedmean_logs)[i]
  PSIG <- phylosig(treeI, CH_I, se = SE_I, test=T)
  phylosignals[i-1,1] <- PSIG$K
  phylosignals[i-1,2] <- PSIG$P
  phylosignals[i-1,3] <- length(CH_I)
  phylosignals$K = round(phylosignals$K,3)
}
#write.csv(phylosignals, "PhylosignalsWSE_log.csv")

#Model support
AICdf = as.data.frame(matrix(ncol=6,nrow=ncol(Q_sharedmean_logs)))
colnames(AICdf) = c("Variable", "white_noise", "starBM", "BM", "EB", "OU")
for(c in 1:ncol(Q_sharedmean_logs)){
  startree <- rescale(ultram, "lambda", 0)
  C = Q_sharedmean_logs[,c]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  startree = drop.tip(startree, which(!(startree$tip.label %in% names(C))))
  Cse = castSE[,c+1] %>% log() %>% abs()
  Cse[which(Cse == Inf)] <- 0
  names(Cse) = castSE$Species
  Cse = Cse[which(names(Cse) %in% names(C))]
  model_matrix = matrix("NA", nrow = 5, ncol = 3)
  colnames(model_matrix) = c("aicc","aicc_best","dAICc")
  row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
  for(j in 1:dim(model_matrix)[1]){
    if(j==2){
      temp_model = fitContinuous(startree, C, model="BM", SE = Cse)$opt
    }
    else{
      temp_model = fitContinuous(Ctree, C, model=row.names(model_matrix)[j], SE = Cse)$opt
    }
    model_matrix = apply(model_matrix,2, as.numeric)
    row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
    model_matrix[j, "aicc"] <- temp_model$aicc
  }
  model_matrix[,"aicc_best"] <- min(model_matrix[,"aicc"])
  model_matrix[,"dAICc"] <- model_matrix[, "aicc"] - model_matrix[j, "aicc_best"]
  print(names(Q_sharedmean_logs)[c])
  string_c <- c(names(Q_sharedmean_logs)[c], model_matrix[,3])
  names(string_c) = colnames(AICdf)
  AICdf[c,] <- string_c
}
AICdf[,2:6] = apply(AICdf[,2:6], 2, as.numeric) %>% apply(2, round3)
#write.csv(AICdf, "log_model_support.csv")

#Model adequacy
worthy_models = AICdf[which(AICdf$white_noise != 0 & AICdf$starBM != 0),]
MAD = as.data.frame(matrix(ncol = 6, nrow = nrow(worthy_models)))
names(MAD) = c("msig", "cvar", "svar", "sasr", "shgt", "dcfd")
rownames(MAD) = worthy_models$Variable
for(m in 1:nrow(worthy_models)){
  C = Q_sharedmean_logs[,which(names(Q_sharedmean_logs) == worthy_models$Variable[m])]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  C = C[match(Ctree$tip.label, names(C))]
  class(Ctree)="phylo"
  FC <- fitContinuous(Ctree,C,model=worthy_models$Best_model[m])
  UTC <- make_unit_tree(FC)
  picstat_data <- calculate_pic_stat(UTC)
  sim <- simulate_char_unit(UTC)
  picstat_sim <- calculate_pic_stat(sim)
  compare_pic_stat(picstat_data, picstat_sim) %>% .$p.values -> MAD[m,]
}
phy_models = data.frame(worthy_models,MAD)
phy_models[,7:12] = apply(phy_models[,7:12], 2,round3)
#write.csv(phy_models, "model_adequacy.csv")

#Nematocyst shape evolution - phylomorphospace figure
het_el = sharedmean_logs$heteroneme_elongation
names(het_el) = sharedmean_logs$Species
hap_el = sharedmean_logs$haploneme_elongation
names(hap_el) = sharedmean_logs$Species
het_el = het_el[which(!is.na(het_el))]
hap_el = hap_el[which(!is.na(hap_el))]
hap_el = hap_el[which(names(hap_el) %in% names(het_el))]
het_el = het_el[which(names(het_el) %in% names(hap_el))]
het_el = het_el[match(names(hap_el), names(het_el))]
elon_data = cbind(het_el, hap_el) %>% as.data.frame()
names(elon_data) = c("Heteroneme elongation (log um)", "Haploneme elongation (log um)")
elon_tree = drop.tip(ultram, which(!(ultram$tip.label %in% names(hap_el))))
phylomorphospace(elon_tree, elon_data, label="horizontal")

## SIMMAPS used to make categorical evolution figure ##
par(ask=F)
tentilla = sharedcategorical$Tentilla
names(tentilla) = rownames(sharedcategorical)
prox_het = sharedcategorical$Proximal.heteronemes
names(prox_het) = rownames(sharedcategorical)
desmo = sharedcategorical$Desmonemes
names(desmo) = rownames(sharedcategorical)
rhopalo = sharedcategorical$Rhopalonemes
names(rhopalo) = rownames(sharedcategorical)
dyn_cnido = sharedcategorical$Dynamic.cnidoband
names(dyn_cnido) = rownames(sharedcategorical)
elastic = sharedcategorical$Elastic.strand
names(elastic) = rownames(sharedcategorical)
distal_desmo = sharedcategorical$Distal.cnidoband.desmonemes
names(distal_desmo) = rownames(sharedcategorical)
coiled = sharedcategorical$Coiled.tentilla
names(coiled) = rownames(sharedcategorical)
heterotype = sharedcategorical$Heteroneme.type
heterotype=as.character(heterotype)
names(heterotype) = rownames(sharedcategorical)
haplotype = sharedcategorical$Haploneme.type
haplotype=as.character(haplotype)
names(haplotype) = rownames(sharedcategorical)
rhopalotype = sharedcategorical$Rhopaloneme.type
rhopalotype=as.character(rhopalotype)
names(rhopalotype) = rownames(sharedcategorical)

Simmap_list = list()
###SIMMAP Tentilla:
make.simmap(ultram_cat, tentilla, nsim = 100) -> tentilla_sim
Simmap_list[[1]] <- tentilla_sim
plotTree(ultram_cat, lwd = 4)
tentilla_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(tentilla_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(tentilla_sim)

###SIMMAP Proximal Heteronemes:
make.simmap(ultram_cat, prox_het, nsim = 100) -> prox_het_sim
Simmap_list[[2]] <- prox_het_sim
plotTree(ultram_cat, lwd = 4)
prox_het_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Present", "Absent")
nodelabels(pie=(describe.simmap(prox_het_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(prox_het_sim)

###SIMMAP Desmonemes:
make.simmap(ultram_cat, desmo, nsim = 100) -> desmo_sim
Simmap_list[[3]] <- desmo_sim
plotTree(ultram_cat, lwd = 4)
desmo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(desmo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(desmo_sim)

###SIMMAP Rhopalonemes:
make.simmap(ultram_cat, rhopalo, nsim = 100) -> rhopalo_sim
Simmap_list[[4]] <- rhopalo_sim
plotTree(ultram_cat, lwd = 4)
rhopalo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(rhopalo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(rhopalo_sim)

###SIMMAP Dynamic Cnidoband:
make.simmap(ultram_cat, dyn_cnido, nsim = 100) -> dyn_cnido_sim
Simmap_list[[5]] <- dyn_cnido_sim
plotTree(ultram_cat, lwd = 4)
dyn_cnido_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(dyn_cnido_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(dyn_cnido_sim)

###SIMMAP Elastic Strand:
make.simmap(ultram_cat, elastic, nsim = 100) -> elastic_sim
Simmap_list[[6]] <- elastic_sim
plotTree(ultram_cat, lwd = 4)
elastic_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(elastic_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(elastic_sim)

###SIMMAP Distal CB Desmonemes:
make.simmap(ultram_cat, distal_desmo, nsim = 100) -> distal_desmo_sim
Simmap_list[[7]] <- distal_desmo_sim
plotTree(ultram_cat, lwd = 4)
distal_desmo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(distal_desmo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(distal_desmo_sim)

###SIMMAP Tentilla Coiledness:
make.simmap(ultram_cat, coiled, nsim = 100) -> coiled_sim
Simmap_list[[8]] <- coiled_sim
plotTree(ultram_cat, lwd = 4)
coiled_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(coiled_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(coiled_sim)

###SIMMAP Heteroneme type:
heterotype = heterotype[heterotype!=""]
HTtree = drop.tip(ultram_cat, which(!(ultram_cat$tip.label %in% names(heterotype))))
make.simmap(HTtree, heterotype, nsim = 100) -> heterotype_sim
Simmap_list[[9]] <- heterotype_sim
plotTree(HTtree, lwd = 4)
heterotype_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red", "green")
names(colors) = c("Birhopaloid", "Microbasic mastigophore", "Stenotele")
nodelabels(pie=(describe.simmap(heterotype_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)

###SIMMAP Haploneme type:
haplotype = haplotype[haplotype!=""]
HTtree = drop.tip(ultram_cat, which(!(ultram_cat$tip.label %in% names(haplotype))))
make.simmap(HTtree, haplotype, nsim = 100) -> haplotype_sim
Simmap_list[[10]] <- haplotype_sim
plotTree(HTtree, lwd = 4)
haplotype_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Anisorhizas", "Isorhizas")
nodelabels(pie=(describe.simmap(haplotype_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(haplotype_sim)

###SIMMAP Rhopaloneme type:
rhopalotype = rhopalotype[rhopalotype!=""]
HTtree = drop.tip(ultram_cat, which(!(ultram_cat$tip.label %in% names(rhopalotype))))
make.simmap(HTtree, rhopalotype, nsim = 100) -> rhopalotype_sim
Simmap_list[[10]] <- rhopalotype_sim
plotTree(HTtree, lwd = 4)
rhopalotype_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Anacrophores", "Acrophores")
nodelabels(pie=(describe.simmap(rhopalotype_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(rhopalotype_sim)

## Character correlations ##
C = sharedmatrix[,c(-1,-3)] %>% .[which(!is.na(rowSums(.[,-1]))),]
Cspecies = C$Species %>% as.character()
C = as.matrix(C[,-1])
rownames(C) = Cspecies
Ctree = drop.tip(ultram,which(!(ultram$tip.label %in% rownames(C))))
class(Ctree) = "phylo"
PICi = Rcontrast(tree=Ctree, X=C, path="phylip-3.695/exe", cleanup=TRUE)

#Correlation visualizations R2
phy_corr = PICi$VarA.Correlations
intra_corr = PICi$VarE.Correlations
phy_corr[upper.tri(phy_corr)] <- NA
combicorr = phy_corr
combicorr[upper.tri(combicorr)] = intra_corr[upper.tri(intra_corr)]
rownames(combicorr) = colnames(C)
colnames(combicorr) = colnames(C)
corrplot(combicorr, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs intraspecific correlations
PICOLS = combicorr
PICOLS[upper.tri(PICOLS)]=cor(C)[upper.tri(cor(C))]
corrplot(PICOLS, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs regular correlations for figure

#Scatterplot phylo vs regular correlations
cbind(as.vector(phy_corr), as.vector(cor(C))) %>% as.data.frame()->phyreg
names(phyreg)<-c("Phylo", "Reg")
ggplot(phyreg, aes(x=Reg, y=Phylo, color=(Phylo+Reg)/2)) + geom_point() + geom_hline(yintercept = 0)  + geom_vline(xintercept = 0) + theme_bw()
abline(h=0)
abline(v=0)

## PCA ##

#Using simple characters
#PCA(raw_matrix_NaZeroes) -> Pca_raw
PCA(raw_matrix_notf) -> Pca_raw
Pca_raw %>% fviz_contrib(choice="var", axes=1, sort.val="desc")
Pca_raw %>% fviz_contrib(choice="var", axes=2, sort.val="desc")
Pca_raw %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #For figure A
Pca_raw %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, axes=c(3,4)) 
raw_tree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(raw_matrix_notf))))
phylomorphospace(tree=raw_tree, Pca_raw$ind$coord[,1:2], label = "horizontal", xlab = "PC1", ylab = "PC2") #For figure B
multiPhylosignal(Pca_raw$ind$coord, raw_tree)
physignal(Pca_raw$ind$coord, raw_tree)

#Using Morphometric characters
PCA(compound_matrix) -> Pca_compound
Pca_compound %>% fviz_contrib(choice="var", axes=1, sort.val="desc")
Pca_compound %>% fviz_contrib(choice="var", axes=2, sort.val="desc")
Pca_compound %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
compound_tree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(compound_matrix))))
phylomorphospace(tree=compound_tree, Pca_compound$ind$coord[,1:2], label = "horizontal", xlab = "PC1", ylab = "PC2")
physignal(Pca_compound$ind$coord, compound_tree)
multiPhylosignal(Pca_compound$ind$coord, compound_tree)

#PhyloPCA on simple character
PPCA_raw = phyl.pca(raw_tree, raw_matrix_notf)
PPCA_compound = phyl.pca(compound_tree, compound_matrix)
PPCA_raw %>% biplot(main=paste(signif(summary(PPCA_raw)$importance[2,1]*100,3),"%"), ylab=paste(signif(summary(PPCA_raw)$importance[2,2]*100,3),"%"), cex = .6, expand =1)
PPCA_compound %>% biplot(main=paste(signif(summary(PPCA_compound)$importance[2,1]*100,3),"%"), ylab=paste(signif(summary(PPCA_compound)$importance[2,2]*100,3),"%"), cex = .6, expand =0.4)
phylomorphospace(raw_tree, PPCA_raw$S, label = "horizontal")

## SIMPLE ANALYSES OF KINEMATIC DATA ##

kine = read.csv("~/tentilla_morph/Scrap/RCode_dependencies/Kinematics.tsv", sep='\t', header=T) #%>% .[which(apply(., 1, function(x) sum(is.na(x)))<ncol(.)-2),-2]
kineWNA = kine[which(kine$Species %in% ultram$tip.label),]
rownames(kineWNA)=kineWNA$Specimen
rownames(kine)=kine$Specimen
#kineWNA = kineWNA[,-1]
#kine=kine[,-1]
kinetree = drop.tip(ultram, which(!(ultram$tip.label %in% kine$Species)))
kinetree$edge.length = 200*kinetree$edge.length
#kine_clean = kineWNA[,which(colSums(kineWNA)>1)]
#names(kine) = c("ADS", "MDS", "HeDS","HeMDS", "HSDS", "HFL", "HaDS")
kine_byspp = aggregate(. ~ kine$Species, data = kine[,c(-1,-2)], mean.na, na.action = na.pass)
names(kine_byspp)[1] <- "Species"
kinemorph <- castmeans[which(castmeans$Species %in% kine_byspp$Species),] %>% cbind(kine_byspp[which(kine_byspp$Species %in% castmeans$Species),])
plot(kinemorph$Cnidoband.free.length..um., kinemorph$Average.CB.discharge.speed..mm.s.)
calys = kinemorph$Species[c(5,9,12,15,16,17,18)]
euphys = kinemorph$Species[which(!(kinemorph$Species %in% calys))]
relspeed = data.frame(kinemorph$Species, c(kinemorph$Average.CB.discharge.speed..mm.s./kinemorph$Cnidoband.free.length..um.))
names(relspeed) = c("Species", "Speed/Length")
relspeed <- relspeed[which(!(is.na(relspeed$`Speed/Length`))),]
t.test(relspeed[which(relspeed$Species %in% calys),2], relspeed[which(relspeed$Species %in% euphys),2])
t.test(kine$Average.CB.discharge.speed..mm.s.[which(kine$Species %in% calys)], kine$Average.CB.discharge.speed..mm.s.[which(kine$Species %in% euphys)])
t.test(castnumbers$Cnidoband.free.length..um.[which(castnumbers$Species %in% calys)], castnumbers$Cnidoband.free.length..um.[which(castnumbers$Species %in% euphys)])
lm(kinemorph$Average.CB.discharge.speed..mm.s.~ kinemorph$Cnidoband.free.length..um.) %>% summary()
lm(kinemorph$Heteroneme.discharge.speed.MAX..mm.s.~ kinemorph$Heteroneme.volume..um3.) %>% summary()
t.test(kine$Heteroneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% calys)], kine$Heteroneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% euphys)])
t.test(kine$Haploneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% calys)], kine$Haploneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% euphys)])
t.test(kine$Haploneme.discharge.speed.AVG..mm.s., kine$Heteroneme.discharge.speed.AVG..mm.s.)
cor(kinemorph[,c(-1,-32)], use="pairwise.complete.obs") %>% .[c(1:30),c(31:41)] %>% corrplot(diag=F, tl.cex = 0.4, tl.col="black")

### Morphospace analyses ###
hypfood <- data.frame(hypdiet, hypdiet, stringsAsFactors = F)
hypfood_full <- data.frame(hypdiet_full, hypdiet_full, stringsAsFactors = F)
#morphood <- raw_matrix_NaZeroes[which(rownames(raw_matrix_NaZeroes) %in% rownames(hypfood)),]
#morphood <- raw_matrix_notf[which(rownames(raw_matrix_notf) %in% rownames(hypfood)),]
#hypfood <- hypfood[match(rownames(morphood), rownames(hypfood)),]
#morphospace_data <- cbind(morphood, hypfood[,1])
morphospace_data = cbind(raw_matrix_NaZeroes, hypfood_full[match(rownames(raw_matrix_NaZeroes), rownames(hypfood_full)),1])
names(morphospace_data)[ncol(morphospace_data)] <- "Diet"
morphospace_data$Diet <- as.character(morphospace_data$Diet) %>% str_replace_all(" ", "_")

#Test for morphospatial overlap between feeding guilds
pca_morph <- prcomp(morphospace_data[,-ncol(morphospace_data)])
autoplot(pca_morph, data = morphospace_data, frame = TRUE, label=T, colour="Diet", loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black")+theme_bw()
autoplot(pca_morph, data = morphospace_data, frame = TRUE, label=T, colour="Diet")+theme_bw()
phylomorphospace(drop.tip(ultram, which(!(ultram$tip.label %in% rownames(morphospace_data)))), pca_morph$x[,1:2], label="horizontal")

#MANOVA for diet
LMdata <- cbind(pca_morph$x, morphospace_data$Diet) %>% as.data.frame(stringsAsFactors = F)
LMdata[,-ncol(LMdata)] <- apply(LMdata[,-ncol(LMdata)], 2, as.numeric)
names(LMdata)[ncol(LMdata)] <- "Diet"
MAN1 <- manova(cbind(PC1, PC2) ~ Diet, data = LMdata)
MAN1 %>% summary.manova() %>% .$SS %>% .$Diet -> MSS1
var_exp = (MSS1[1,2]+MSS1[2,1])/sum(MSS1)

#phyl.MANOVA
LMtree <- drop.tip(ultram, which(!(ultram$tip.label %in% rownames(LMdata))))
dat <- LMdata[,1:2]
D <- LMdata$Diet %>% as.factor()
names(D) <- rownames(LMdata)
PMAN <- aov.phylo(dat ~ D, phy = LMtree, nsim=50, test="Wilks")
PMAN  %>% summary.manova() %>% .$SS %>% .$group -> pMSS
var_exp_phyl = (pMSS[1,2]+pMSS[2,1])/sum(pMSS)

#geomorph disparity
Di <- D[!is.na(D)]
dati <- dat[which(rownames(dat) %in% names(Di)),]
tree_i <- drop.tip(ultram, which(!(ultram$tip.label %in% rownames(dati))))
morphol.disparity(dati ~ 1, groups = Di)
gdf <- geomorph.data.frame(shape = as.matrix(dati), phy = tree_i, grp = Di)
ANOVA <- procD.pgls(shape ~ grp, phy = phy, data = gdf)
morphol.disparity(f1 = ANOVA, groups = ~grp, data = gdf, iter = 999)

###FULL PCA DIET ###
dpruned_full <- dpruned_full[which(!(dpruned_full$Species %in% c("Nectadamas richardi", "Thermopalia taraxaca", "Halistemma foliacea", "Halistemma transliratum", "Halistemma cupulifera", "Resomia dunni", "Resomia persica", "Cardianecta parchelion", "Forskalia tholoides"))),]
dpruned_full_raw <- dpruned_full[,c(1:21,33)]
pca_fullmorph <- prcomp(dpruned_full_raw[,c(-1,-ncol(dpruned_full_raw))])
autoplot(pca_fullmorph, data = dpruned_full_raw, frame = TRUE, label=T, colour="Diet", loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black")+theme_bw()

autoplot(pca_fullmorph, data = dpruned_full_raw, frame = TRUE, label=T, colour="Diet")+theme_bw()

phylomorphospace(drop.tip(ultram, which(!(ultram$tip.label %in% rownames(dpruned_full_raw)))), pca_fullmorph$x[which(rownames(pca_fullmorph$x) %in% ultram$tip.label),1:2], label="horizontal")

#full MANOVA
full_LMdata <- cbind(pca_fullmorph$x, dpruned_full$Diet) %>% as.data.frame(stringsAsFactors = F)
full_LMdata[,-ncol(full_LMdata)] <- apply(full_LMdata[,-ncol(full_LMdata)], 2, as.numeric)
names(full_LMdata)[ncol(full_LMdata)] <- "Diet"
full_MAN1 <- manova(cbind(PC1, PC2) ~ Diet, data = full_LMdata)
full_MAN1 %>% summary.manova() %>% .$SS %>% .$Diet -> full_MSS1
full_var_exp = (full_MSS1[1,2]+full_MSS1[2,1])/sum(full_MSS1)

#phyl.MANOVA full
full_LMtree <- drop.tip(ultram, which(!(ultram$tip.label %in% rownames(full_LMdata))))
full_dat <- full_LMdata[,1:2]
full_D <- full_LMdata$Diet %>% as.factor()
names(full_D) <- rownames(full_LMdata)
full_PMAN <- aov.phylo(full_dat ~ full_D, phy = full_LMtree, nsim=50, test="Wilks")
full_PMAN  %>% summary.manova() %>% .$SS %>% .$group -> full_pMSS
full_var_exp_phyl = (full_pMSS[1,2]+full_pMSS[2,1])/sum(full_pMSS)

#geomorph disparity
full_Di <- full_D[!is.na(full_D)]
full_Di <- full_Di[which(names(full_Di) %in% rownames(full_dat) & names(full_Di) %in% ultram$tip.label)]
full_dati <- full_dat[which(rownames(full_dat) %in% names(full_Di) & rownames(full_dat) %in% ultram$tip.label),]
full_tree_i <- drop.tip(ultram, which(!(ultram$tip.label %in% rownames(full_dati))))
morphol.disparity(full_dati ~ 1, groups = full_Di)
full_gdf <- geomorph.data.frame(shape = as.matrix(full_dati), phy = full_tree_i, grp = full_Di)
full_ANOVA <- procD.pgls(shape ~ grp, phy = phy, data = full_gdf)
morphol.disparity(f1 = full_ANOVA, groups = ~grp, data = full_gdf, iter = 999)

##SURFACE##

surfaceALL <- function(data){
  Tree <- nameNodes(drop.tip(ultram, which(!(ultram$tip.label %in% rownames(data)))))
  olist <- convertTreeData(Tree, data)
  otree<-olist[[1]]
  odata<-olist[[2]]
  fwd<-surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, verbose = FALSE, plotaic = FALSE)
  k<-length(fwd)
  fsum<-surfaceSummary(fwd)
  bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
  bsum<-surfaceSummary(bwd)
  kk<-length(bwd)
  print("N regimes BWD")
  bsum$n_regimes %>% return()
  par(mfrow=c(1,2))
  surfaceTreePlot(Tree, bwd[[kk]], labelshifts = T) %>% return()
  surfaceTraitPlot(data, bwd[[kk]], whattraits = c(1,2)) %>% return()
  newsim<-surfaceSimulate(Tree, type="hansen-fit", hansenfit=fwd[[k]]$fit, shifts=fwd[[k]]$savedshifts, sample_optima=TRUE)
  newout<-runSurface(Tree, newsim$dat, only_best = TRUE)
  newsum<-surfaceSummary(newout$bwd)
  newkk<-length(newout$bwd)
  print("N regimes SIM")
  newsum$n_regimes %>% return()
}

raw_matrix_NaZeroes[,6:7]
pcaphenos$ind$coord[,1:2] %>% as.data.frame()
PCA(raw_matrix_NaZeroes)$ind$coord[,1:2] %>% as.data.frame()

surfaceALL(raw_matrix_NaZeroes[-which(rownames(raw_matrix_NaZeroes) %in% c("Apolemia lanosa", "Apolemia rubriversa")),6:7]) #haploneme shape
surfaceALL(pcaphenos$ind$coord[,1:3] %>% as.data.frame()) #trait values PCA 1 2 without TF, excluding some species
surfaceALL(PCA(raw_matrix_NaZeroes)$ind$coord[,1:2] %>% as.data.frame()) #ALL TRAITS PC1 and 2, all SPP in tree, NA to zeroes ALL TRAITS PC1 and 2, all SPP in MORPH dataset, NA to zeroes
surfaceALL(raw_matrix_NaZeroes[,8:11]) #desmonemes rhopalonemes?

#Phenotypic integration test

M1_het <- sharedmean_logs[c(2:6)]
rownames(M1_het) <- sharedmean_logs$Species
M1_het <- M1_het[which(!is.na(rowSums(M1_het))),]

M2_hap <- sharedmean_logs[c(8:9)]
rownames(M2_hap) <- sharedmean_logs$Species
M2_hap <- M2_hap[which(!is.na(rowSums(M2_hap))),]

M3_tf <- sharedmean_logs[c(10:13)]
rownames(M3_tf) <- sharedmean_logs$Species
M3_tf <- M3_tf[which(!is.na(rowSums(M3_tf))),]

M4_tent <- sharedmean_logs[c(14:21)]
rownames(M4_tent) <- sharedmean_logs$Species
M4_tent <- M4_tent[which(!is.na(rowSums(M4_tent))),]

phint <- function(MA,MB,Tr){
  a <- MA[which(rownames(MA) %in% rownames(MB)),]
  b <- MB[which(rownames(MB) %in% rownames(a)),]
  t <- drop.tip(Tr, which(!(Tr$tip.label %in% rownames(a) & Tr$tip.label %in% rownames(b))))
  phylo.integration(A = as.matrix(a), A2 = as.matrix(b), phy = t) %>% return()
}

table.phint <- function(modules, tree){
  I <- length(modules)
  R <- data.frame(matrix(ncol=I, nrow=I))
  rownames(R) = names(modules)
  names(R) = names(modules)
  P <- data.frame(matrix(ncol=I, nrow=I))
  rownames(P) = names(modules)
  names(P) = names(modules)
  for(i in 1:I){
    for(j in 1:I){
      if(i != j){
        ph_ij <- phint(modules[[i]], modules[[j]], Tr = tree)
        R[i,j] <- ph_ij$r.pls
        P[i,j] <- ph_ij$P.value
      }
      print(R);print(P)
    }
    print(R); print(P)
  }
  Rm <- melt(cbind(names(modules),R))
  names(Rm) <- c("M1", "M2", "rPLS")
  Pm <- melt(cbind(names(modules),P))
  names(Pm) <- c("M1", "M2", "P.value")
  return(full_join(Rm,Pm))
}

MODS <- list(M1_het, M2_hap, M3_tf, M4_tent)
names(MODS) <- c("Heteronemes", "Haplonemes", "Terminal filament", "Tentillum")
MODSphint <- table.phint(MODS, tree = ultram)

ordering <- c("Tentillum","Heteronemes","Haplonemes", "Terminal filament")
ggplot(MODSphint, aes(factor(M1, levels = ordering), factor(M2, levels = ordering))) + 
  geom_tile(aes(fill = rPLS), colour = "black") + 
  scale_fill_distiller(palette = "Blues", direction=1) + 
  geom_text(aes(label = P.value %>% round(3))) +
  labs(y= "", x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))