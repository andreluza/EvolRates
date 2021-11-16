
# function to check where uncertainty is larger (either phylogeny or ancestral reconstruction)
# based on the resampling of data within and between phylogenetic trees


fc_tr <- function (min,max,nsample,data){
  seq1<- seq (min,max)
  seq2<-seq (min,10000,100)
  s1 <- sample (seq1,nsample,replace=T)
  s2 <- sample (seq2,nsample,replace=T)
  out <- data.frame (
    sd_within =sd(sapply (data[s1],"[[",14)),
    sd_between = sd(sapply (data[s2],"[[",14)));
  out
  
}

# neighborhood function
mst.nb <-
  function(filt1){
    mymst <- ade4::mstree(filt1)
    return(neig2nb(mymst))
  }
