# -------------------------------------------------
# step 1, estimating transition rates, per species
# -------------------------------------------------

#load packages and function
source("R/packages.R")
source("R/transition_rates_function.R")

# 1 - Load phylogenetic trees
tree_list <- tree <- read.nexus(here ("Data","PhylogenyTraits","Sigmodontinae_413species100Trees.trees"))
## corrigir os nomes
tree_list <- lapply (tree_list, function (i)
	{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# 2 - Load trait data
atributos_imp <- readRDS (here ("Data","PhylogenyTraits","trait_data.rda"))

# 3 - defining diets
# Insectivorous
insEat <- atributos_imp [atributos_imp$Diet_Inv >= 50 &
	atributos_imp$Diet_PlantO <= 49 & 
	atributos_imp$FrGr <= 49,]
# Plant eaters
plantEat <- atributos_imp [atributos_imp$Diet_PlantO >= 50 & 
	atributos_imp$Diet_Inv <= 49 &
	atributos_imp$FrGr <= 49,]
# Fruit eaters
fruitEat <- atributos_imp [atributos_imp$FrGr >= 50 & 
	atributos_imp$Diet_Inv <= 49 & 
	atributos_imp$Diet_PlantO <= 49,]

# discrete scle
diet_disc<-matrix (c (rep ("fruitEat", length (rownames(fruitEat))),
                      rep ("plantEat", length (rownames(plantEat))) , 
	                    rep ("insEat", length (rownames(insEat)))), 
                   dimnames =list(c(rownames(fruitEat), 
                                    rownames(plantEat), 
                                    rownames(insEat)),"diet")) 

## species not included in the previous categories
general <- atributos_imp [which ((rownames(atributos_imp) %in% rownames(diet_disc)) == FALSE),]
# bind them
diet_disc <- rbind (diet_disc, matrix (rep ("general", length(rownames (general))), dimnames =list(rownames(general), "diet")))
# traits as vector
trait <- as.vector (diet_disc)
names(trait) <- rownames(diet_disc) # set names
trait <- trait [order (names(trait))] # order

# number of simulations (ancestral character estimation) to run
nsim <- 2

# run transition rates function 
cl <- makeCluster(4) ## number of cores
# export packages to the cores
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(daee))

# export your data and function
clusterExport(cl, c("trait",
                    "tree_list",
                    "tip.based.trait.evo",
                    "nsim"))

## you must apply nsim to 1 one phylogeny each time. Thus, the list argument contains the phylogenies as class multiphlylo
transitions_estimates_ages<- parLapply(cl, as.list(seq(1,length (tree_list))), function (i)
  
  tip.based.trait.evo (tree=tree_list[[i]], 
                       trait=trait, 
                       nsim=nsim,
                    
                    # the three methods developed in this study 
                    method = c("transition_rates", 
                               "last_transition_time", 
                               "stasis_time")) 
    )

stopCluster(cl)

# save estimated metrics
saveRDS (transitions_estimates_ages, 
         here ("Output","res_step1_transitions_estimates_ages.rda"))

# descriptive statistics (reported in the results)
mean_across_species <- lapply (transitions_estimates_ages, function (i)
	do.call(rbind,lapply (i, function (k) colMeans(k))))

apply(do.call(rbind,mean_across_species),2,mean)
apply(do.call(rbind,mean_across_species),2,sd)

