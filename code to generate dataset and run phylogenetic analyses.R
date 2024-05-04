library(ape)
library(phytools)
library(geiger)
library(brms)
library(rcartocolor)

################################################################################################################################
############################################# get the tree and data ############################################################
################################################################################################################################

# TREE
tree_full <- read.tree("rolland_2014.txt")
tree.full.species<-c("Artibeus_jamaicensis", "Cebus_capucinus", "Macaca_fascicularis", "Myotis_bechsteinii", "Orcinus_orca", "Crocuta_crocuta", "Cryptomys_damarensis", "Eptesicus_fuscus", "Pan_troglodytes", "Canis_lupus", "Cercopithecus_mitis", "Macaca_mulatta", "Papio_cynocephalus", "Marmota_flaviventris", "Vulpes_vulpes", "Ctenodactylus_gundi", "Cynomys_ludovicianus")
ii<-sapply(tree.full.species,grep,tree_full$tip.label)
ii
tree.species.analysis<-tree_full$tip.label[ii]
tree.species.analysis
tree<-drop.tip(tree_full, setdiff(tree_full$tip.label,tree.species.analysis))
is.ultrametric(tree) # the tree is currently not ultrametric
tree_UM<-force.ultrametric(tree, method="extend") # we force the tree to be ultrametric
is.ultrametric(tree_UM) # the tree is ultrametric
is.binary.tree(tree_UM) # the tree is bifurcating
tree<- tree_UM # rename tree for simplicity
plot(tree)
kinship_composition <- c("mixrelated", "related", "related", "related", "related", "mixrelated", "mixrelated", "related", "mixrelated", "mixrelated", "mixrelated", "mixrelated", "mixrelated", "related", "related", "related", "related" ) # create a vector with the state of the kinship composition for each species (in the order of the species according to the phylogeny)
names(kinship_composition) <- tree$tip.label  # name each state of the character according  to the species in the phylogeny
cols<-setNames(c("#44AA99","#332288"),c("mixrelated", "related"))
kinship_composition # Confirm that everything was assigned correctly

################################################################################################################################
##################################################### stochastic maps ##########################################################
################################################################################################################################

##################################################### model choice #############################################################
################################################################################################################################

SYM<-fitDiscrete(tree,kinship_composition, model="SYM")
SYM # SYM and ER are the same here, so we exclude ER.
ARD<-fitDiscrete(tree,kinship_composition, model="ARD")
ARD

aic.vals<-setNames(c(SYM$opt$aicc,ARD$opt$aicc),
                   c("SYM","ARD"))
aic.vals
aic.w(aic.vals)
# We choose the SYM model

##################################################### model #############################################################
##########################################################################################################################

mtrees_SYM<- make.simmap(tree, kinship_composition, model= "SYM", nsim=10000, Q="mcmc", prior= list(use.empirical=TRUE))
mtrees_SYM
pd_SYM<- summary(mtrees_SYM,plot=FALSE)
pd_SYM

##################################################### figure 1 ###########################################################
##########################################################################################################################

windowsFonts(A = windowsFont("Arial Unicode MS"))
par(family = "A", ps=12)
obj<-densityMap(mtrees_SYM,lwd=4,outline=TRUE)
svg("figure1.svg")
n<-length(obj$cols) #lenght of our colour ramp 
obj$cols[1:n]<-colorRampPalette(c("#44AA99", "#332288"), space="Lab")(n)
plot(obj,par(family = "A", ps=12))
plot(obj, no.margin=TRUE, edge.width=2)
nodelabels(pie=pd_SYM$ace,piecol=c("#44AA99", "#332288"), cex=0.4)
legend("bottomright", 
       legend=c("Related", "Mix-related"), 
       pch=c(21,21),
       pt.bg=c("#332288","#44AA99"),
       inset=0.01,
       box.lty=0,
       par(family = "A", ps=12))
dev.off()

##################################################### figure S1 (A) ###########################################################
##########################################################################################################################

mtrees_SYM_1<- mtrees_SYM[sample(1)]
windowsFonts(A = windowsFont("Arial Unicode MS"))
svg(file="figureS1_A.svg")
par(family = "A", ps=12)
dotTree(tree, kinship_composition, ftype="i", colors=cols, legend=FALSE)
par(family = "A", ps=12)
legend("bottomleft", 
       legend=c("Related -> Mix-related", "Mix-related -> Related"), 
       pch=c(22,22),
       pt.bg=c("#44AA99", "#332288"),
       inset=.001, pt.cex=2,
       box.lty=0)
legend("bottomright", 
       legend=c("Related", "Mix-related"), 
       pch=c(21,21),
       pt.bg=c("#332288","#44AA99"),
       inset=.001, pt.cex=2,
       box.lty=0)
nulo<-sapply(mtrees_SYM_1,markChanges, cols)
dev.off()

##################################################### figure S1 (B) ###########################################################
##########################################################################################################################

mtrees_SYM_100<- mtrees_SYM[sample(100)]
windowsFonts(A = windowsFont("Arial Unicode MS"))
svg(file="figureS1_B.svg")
par(family = "A", ps=12)
dotTree(tree, kinship_composition, ftype="i", colors=cols, legend=FALSE)
par(family = "A", ps=12)
legend("bottomleft", 
       legend=c("Related -> Mix-related", "Mix-related -> Related"), 
       pch=c(22,22),
       pt.bg=c("#44AA99", "#332288"),
       inset=.001, pt.cex=2,
       box.lty=0)
legend("bottomright", 
       legend=c("Related", "Mix-related"), 
       pch=c(21,21),
       pt.bg=c("#332288","#44AA99"),
       inset=.001, pt.cex=2,
       box.lty=0)
nulo<-sapply(mtrees_SYM_100,markChanges, cols)
dev.off()

##################################################### figure S1 (C) ###########################################################
##########################################################################################################################

windowsFonts(A = windowsFont("Arial Unicode MS"))
par(family = "A", ps=12)
svg(file="figure2.svg")
dotTree(tree, kinship_composition, ftype="i", colors=cols, legend=FALSE)
par(family = "A", ps=12)
legend("bottomleft", 
       legend=c("Related -> Mix-related", "Mix-related -> Related"), 
       pch=c(22,22),
       pt.bg=c("#44AA99", "#332288"),
       inset=.001, pt.cex=2,
       box.lty=0)
legend("bottomright", 
       legend=c("Related", "Mix-related"), 
       pch=c(21,21),
       pt.bg=c("#332288","#44AA99"),
       inset=.001, pt.cex=2,
       box.lty=0)
nulo<-sapply(mtrees_SYM,markChanges, cols)
dev.off()

#################################################### Mammal tree for figure S2 ###########################################
##########################################################################################################################

windowsFonts(A = windowsFont("Arial Unicode MS"))
svg(file="figureS2.svg")
tree_full <- read.tree("rolland_2014.txt")
tree.full.species<-c("Artibeus_jamaicensis", "Cebus_capucinus","Orcinus_orca", "Crocuta_crocuta", "Cryptomys_damarensis", "Zaglossus_bruijni","Glironia_venusta","Caenolestes_caniventer","Dromiciops_gliroides","Notoryctes_caurinus","Sarcophilus_harrisii","Isoodon_auratus","Dendrolagus_inustus","Orycteropus_afer","Geogale_aurita","Elephantulus_intufi","Dendrohyrax_dorsalis","Elephas_maximus","Trichechus_inunguis","Calyptophractus_retusus","Choloepus_didactylus","Anathana_ellioti","Cynocephalus_volans","Pentalagus_furnessi","Desmana_moschata","Manis_crassicaudata","Rhinoceros_unicornis") # We use a species per Order of mammals
ii<-sapply(tree.full.species,grep,tree_full$tip.label)
ii
tree.species.figs2<-tree_full$tip.label[ii]
tree.species.figs2
tree_figures2<-drop.tip(tree_full, setdiff(tree_full$tip.label,tree.species.figs2))
new_tiplabels <- c("Monotremata", "Didelphimorphia", "Paucituberculata", "Peramelemorphia", "Notoryctemorphia", "Dasyuromorpha", "Microbiotheria", "Diprotodontia", "Proboscidea", "Sirenia", "Hyracoidea", "Tubulidentata", "Macroscelidea", "Afrosoricida", "Cingulata", "Pilosa", "Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Artiodactyla", "Chiroptera", "Scandentia", "Dermoptera", "Primata", "Lagomorpha", "Rodentia") # Rename the tips of the tree for the Order that each species in the tree represents
tree_figures2$tip.label <- new_tiplabels
plot(tree_figures2)
dev.off()

################################################################################################################################
##################################################### BRMS models ##########################################################
################################################################################################################################

##################################################### get the data #############################################################
################################################################################################################################

species_data_female<-read.csv("kinship_composition.csv", header=T) # data are provided in a supplementary table, for illustration purposes we assume the data are available as a csv file
str(species_data_female)
species_data_female$kinship_composition<-as.character(species_data_female$kinship_composition)
species_data_female$litter_size<-as.character(species_data_female$litter_size)
species_data_female$phylogeny<-as.character(species_data_female$phylogeny)
species_data_female$unit_size<-as.numeric(as.character(species_data_female$unit_size))

################################################# Prepare tree #################################################################
################################################################################################################################

# note: the tree is called "full" but it doesn't include Artibeus_jamaicensis because we do not have female Artibeus_jamaicensis data.
tree_full <- read.tree("rolland_2014.txt")
tree.full.species<-c("Cebus_capucinus", "Macaca_fascicularis", "Myotis_bechsteinii", "Orcinus_orca", "Crocuta_crocuta", "Cryptomys_damarensis", "Eptesicus_fuscus", "Pan_troglodytes", "Canis_lupus", "Cercopithecus_mitis", "Macaca_mulatta", "Cuon_alpinus", "Papio_cynocephalus", "Marmota_flaviventris", "Vulpes_vulpes", "Ctenodactylus_gundi", "Cynomys_ludovicianus")
ii<-sapply(tree.full.species,grep,tree_full$tip.label)
ii
tree.species.analysis<-tree_full$tip.label[ii]
tree.species.analysis
final_tree_full<-drop.tip(tree_full, setdiff(tree_full$tip.label,tree.species.analysis))
final_tree_full
is.ultrametric(final_tree_full) # the tree currently is not ultrametric
tree_UM_full<-force.ultrametric(final_tree_full, method="extend")
is.ultrametric(tree_UM_full) # the tree is ultrametric

################################################# Prepare phylogeny for model ##################################################
################################################################################################################################

species_full<-tree_UM_full
A <- ape::vcv.phylo(species_full)

######################################################### MODELS ###############################################################
################################################################################################################################

######################################################## Unit size #############################################################
################################################################################################################################

prior1<- get_prior(kinship_composition ~ unit_size + (1|gr(phylogeny, cov = A)), data=species_data_female, family=bernoulli(), data2= list(A = A))
prior1
unit_size_model <- brm(
  kinship_composition ~ unit_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c( 
    prior(normal(0,10), class="b", coef="unit_size"),
    prior(student_t(3,0,2.5), "Intercept"), # from get_prior
    prior(student_t(3,0,2.5), "sd") # from get_prior
  ),
  control = list(adapt_delta = 0.999, max_treedepth=15)
)

plot(unit_size_model) 
pp1= pp_check(unit_size_model, ndraw=1000)
pp1 
summary(unit_size_model) 

###
# Check if priors were appropriate

unit_size_model_prior1 <- brm(
  kinship_composition ~ unit_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c(
    prior(normal(0,5), class="b", coef="unit_size"),
    prior(student_t(3,0,1.25), "Intercept"),
    prior(student_t(3,0,1.25), "sd") 
  ),
  control = list(adapt_delta = 0.999, max_treedepth=15)
)
summary(unit_size_model_prior1)

unit_size_model_prior2 <- brm(
  kinship_composition ~ unit_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c(
    prior(normal(0,20), class="b", coef="unit_size"),
    prior(student_t(3,0,5), "Intercept"), 
    prior(student_t(3,0,5), "sd")
  ),
  control = list(adapt_delta = 0.999, max_treedepth=15)
)
summary(unit_size_model_prior2)

######################################################## Litter size #################################################################
################################################################################################################################

prior2<- get_prior(kinship_composition ~ litter_size + (1|gr(phylogeny, cov = A)), data=species_data_female, family=bernoulli(), data2= list(A = A))
prior2
litter_size_model <- brm(
  kinship_composition ~ litter_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c(
    prior(normal(0,10), class="b", coef="litter_size1"),
    prior(student_t(3,0,2.5), "Intercept"), # from get_prior
    prior(student_t(3,0,2.5), "sd")), # from get_prior
  control = list(adapt_delta = 0.999, max_treedepth=15)
)

plot(litter_size_model) 
pp2= pp_check(litter_size_model, ndraw=1000)
pp2 
summary(litter_size_model)

###
# Check if priors were appropriate

litter_size_model_prior1 <- brm(
  kinship_composition ~ litter_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c( 
    prior(normal(0,5), class="b", coef="litter_size1"),
    prior(student_t(3,0,1.25), "Intercept"), 
    prior(student_t(3,0,1.25), "sd")), 
  control = list(adapt_delta = 0.999, max_treedepth=15)
)
summary(litter_size_model_prior1)

litter_size_model_prior2 <- brm(
  kinship_composition ~ litter_size + (1|gr(phylogeny, cov = A)), 
  data = species_data_female, 
  family = bernoulli(),
  data2 = list(A = A), 
  prior= c(
    prior(normal(0,20), class="b", coef="litter_size1"),
    prior(student_t(3,0,5), "Intercept"), 
    prior(student_t(3,0,5), "sd")), 
  control = list(adapt_delta = 0.999, max_treedepth=15)
)
summary(litter_size_model_prior2)