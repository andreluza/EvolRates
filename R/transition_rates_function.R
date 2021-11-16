# FUNCTION TO ESTIMATE THE THREE METRICS OF TRAIT TRANSITION RATES PRESENTED IN LUZA ET AL. 2021. ECOLOGY AND EVOLUTION
# the function was developed by Andre L. Luza and Vanderlei J. Debastiani


### Funtion to calculate tip-based metrics of trait evolution ###

## Arguments

# tree - a phylogenetic tree as an object of class "phylo"
# traits - a named vector containing the tip states for a discretely valued character. The names must match the tip labels of the tree.
# nsim - number of simulations to stochastic character maps.
# method - tip-based metric, partial match to "transition_rates", "last_transition_time" and "stasis_time".

## Value

# A list (length equal to nsim) with tip-based metrics by species.

tip.based.trait.evo <- function(tree, trait, nsim = 1, method = c("transition_rates", "last_transition_time", "stasis_time")) {
  ## Internal function
  adjacency.tree <- function(tree){
    temp <- ape::prop.part(tree)
    result <- matrix(0, nrow = length(tree$tip), ncol = length(temp), dimnames = list(tree$tip.label, tree$node.label))
    for(i in 1:ncol(result)){
      result[temp[[i]],i] <- 1
    }
    return(result)	
  }
  tree <- daee::node.tree(tree)$tree
  n.sp <- ape::Ntip(tree)
  n.no <- ape::Nnode(tree)
  spxnode <- adjacency.tree(tree)
  spxnode.edge <- matrix(NA, n.sp, n.no, dimnames = list(rownames(spxnode), colnames(spxnode)))
  spxnode.length <- matrix(NA, n.sp, n.no, dimnames = list(rownames(spxnode), colnames(spxnode)))
  for(i in 1:n.no){
    temp <- daee::tree.label.info(tree, colnames(spxnode)[i])
    spxnode.edge[spxnode[,i]==1, i] <- temp$edge
    spxnode.length[spxnode[,i]==1, i] <- temp$edge.length
  }
  sp.edge <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "Edge"))
  sp.length <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "Length"))
  for(i in 1:n.sp){
    temp <- daee::tree.label.info(tree, rownames(spxnode)[i])
    sp.edge[i,1] <- temp$edge  
    sp.length[i,1] <- temp$edge.length  
  }
  ## stochastic mapping of discrete traits via Bayesian inference 
  result.scm <-  phytools::make.simmap(tree, trait, model = "SYM", type = "discrete", nsim = nsim, message = FALSE, pi = "estimated", Q = "mcmc")
  ## Internal function 
  f.int <- function(scm, spxnode, sp.edge, trait, n.sp, n.no, method){
    states <- setNames(sapply(scm$maps, function(x) names(x)[1]), scm$edge[, 1])
    states <- states[as.character(Ntip(scm) + 1:scm$Nnode)]
    sp.node.trait <- sweep(spxnode, 2, STATS = cbind(states), function(x, z) ifelse(x==1, z, NA))
    METHOD <- c("transition_rates", "last_transition_time", "stasis_time")
    method <- pmatch(method, METHOD)
    if (any(is.na(method))) {
      stop("\n Invalid method \n")
    }
    result <- data.frame(row.names = rownames(spxnode))
    if(any(method==1)){
      transitions <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "transition"))
      for(i in 1:n.sp){
        base.temp <- trait[rownames(spxnode)[i]]
        n.tra.temp <- 0
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            if(sp.node.trait[i,j]!=base.temp){
              n.tra.temp <- n.tra.temp+1
              base.temp <- sp.node.trait[i,j]
            }
          }
        }
        transitions[i,1] <- n.tra.temp
      }
      total.nodes <- cbind(apply(spxnode, 1, sum))
      prop.transitions <- transitions/total.nodes
      colnames(prop.transitions) <- "prop.transitions"
      result$transitions <- transitions
      result$total.nodes <- total.nodes
      result$prop.transitions <- prop.transitions
    }
    if(any(method==2)){
      time.transition <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "last.transition.time"))
      for(i in 1:n.sp){
        base.temp <- trait[rownames(spxnode)[i]]
        time.tra.temp <- sp.length[i,1]
        go <- TRUE
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            if(sp.node.trait[i,j]!=base.temp){
              go <- FALSE
            } else{
              if(go){
                time.tra.temp <- time.tra.temp+spxnode.length[i,j]
              }
            }
          }
        }
        time.transition[i,1] <- time.tra.temp
      }
      result$last.transition.time <- time.transition
    }
    if(any(method==3)){
      scm.maps.max <- sapply(scm$maps, function(x) x[which(x==max(x))], simplify = FALSE)
      stasis.time <- matrix(NA, n.sp, 1, dimnames=list(rownames(spxnode), "stasis.time"))
      for(i in 1:n.sp){
        time.sta.temp <- scm.maps.max[[sp.edge[i,1]]]
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            time.sta.temp <- max(time.sta.temp, scm.maps.max[[spxnode.edge[i,j]]])
          }
        }
        stasis.time[i,1] <- time.sta.temp
      }
      result$stasis.time <- stasis.time
    }
    return(result)
  }
  if(nsim==1){
    RES <- sapply(list(result.scm), f.int, spxnode = spxnode, sp.edge = sp.edge, trait = trait, n.sp = n.sp, n.no = n.no, method = method, simplify = FALSE)
  } else{
    RES <- sapply(result.scm, f.int, spxnode = spxnode, sp.edge = sp.edge, trait = trait, n.sp = n.sp, n.no = n.no, method = method, simplify = FALSE)
  }
  return(RES)
}

### Working Example ###

## Phylogenetic tree
#tree <- read.tree(text = "(s1:0.9,((s2:0.5,(s3:0.1,s4:0.1)n4:0.4)n3:0.1,(s5:0.5,s6:0.5)n5:0.1)n2:0.3)n1;")
## Trait states
#trait <- c("General", "Fruit", "Fruit", "Fruit", "Plant", "Plant")
#names(trait) <- tree$tip.label
#
## Plot
#color <- c("gray", "white", "black")
#plot.phylo(tree, label.offset = 0.05, show.node.label = TRUE)
#tiplabels(pch = 22, bg = color[as.numeric(as.factor(trait))], cex = 2, adj = 0.5)
#legend("topleft", fill = color, legend = levels(as.factor(trait)), bty = "n")
## Tip-based metrics
#results <- tip.based.trait.evo(tree, trait, nsim = 2, method = c("transition_rates", "last_transition_time", "stasis_time"))
#results
#