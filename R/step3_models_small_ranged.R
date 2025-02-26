
# ---------------------------------
#        GLMMM
# ---------------------------------

#load packages and functions
source("R/packages.R")
source("R/other_functions.R")

# load modeling data
load(here ("Output","small_ranges_output","res_step2_data_for_GLMM_small_ranged.RData"))

## settings
ctrl <- lmeControl(opt='optim',
                   maxIter =  1e+08, 
                   msMaxIter =  1e+08)

# too hard have estimates for each phylogeny
# lets get a subset of 1000 samples
sample_estimates <- sample(seq(1,length(dataTR)),2000)

# and save it to have the record
save(sample_estimates, 
     file=here ("Output","small_ranges_output","res_step3_sample_estimates_small_ranged.RData"))
# load if needed in the future
load(file=here ("Output","small_ranges_output","res_step3_sample_estimates_small_ranged.RData"))

# obtain subsets of the complete datasets
dataTR_sub <- dataTR[sample_estimates]

# subtly jitter coordinates to avoid errors in GLMM
dataTR_sub <- lapply (dataTR_sub, function (i) {
  
      i$coordinates.V2_jitter <- jitter(i$coordinates.V2,1)
      ;
      i
})

# ncores to use throughout analyses
ncores<- 6

# complete model with spatial autocorrelation term

cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataTR_sub","ctrl"))

teste_exp_nugget  <- parLapply (cl,dataTR_sub, function (i) {

	tryCatch(
		lme ((TR) ~ position*tipo_hab+
			            neigh_area + 
		              over_FOR+
		              over_AB+ 
				          mata_atl + 
		              andes_eco, 
		     
			 data = i,
			 control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			 correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=T)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)

stopCluster (cl)

table(unlist(lapply(teste_exp_nugget,is.character)))
save(teste_exp_nugget,file=here ("Output","small_ranges_output","res_step3_model_exp_nugget_TR_small_ranged.RData"))

# --------------------------------------------
# MODELO WITH AUTOCORRELATION BUT NO NUGGET EFFECT


cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataTR_sub","ctrl"))

teste_exp  <- parLapply (cl, dataTR_sub, function (i) {

	tryCatch(
		lme ((TR) ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=F)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)

stopCluster (cl)
# check convergence issues
# table (unlist (lapply (teste_exp, is.character)))
save(teste_exp,file=here ("Output","small_ranges_output","res_step3_model_exp_TR_small_ranged.RData"))

### modelo sem autoc espacial

cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataTR_sub","ctrl"))

no_aut <- parLapply (cl, dataTR_sub, function (i) {

	tryCatch(
		lme ((TR) ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML"),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)

stopCluster (cl)
# check convergence issues
#table (unlist (lapply (no_aut, is.character)))
save(no_aut, file=here ("Output","small_ranges_output","res_step3_modelNoaut_TR_small_ranged.RData"))

###################################
#### STASIS TIME

# obtain subsets of the complete datasets
dataST_sub <- dataST[sample_estimates]

# subtly jitter coordinates to avoid errors in GLMM
dataST_sub <- lapply (dataST_sub, function (i) {
  
  i$coordinates.V2_jitter <- jitter(i$coordinates.V2,1)
  ;
  i
})

# modelo com autocorrelacao e nugget
cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataST_sub","ctrl"))

teste_exp_nugget <- parLapply (cl, dataST_sub, function (i) {

	tryCatch(
		lme (ST ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=T)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)

stopCluster (cl)

## no convergence problem!
table (unlist (lapply (teste_exp_nugget, is.null)))
save (teste_exp_nugget, file=here ("Output","small_ranges_output","res_step3_teste_exp_nugget_ST_small_ranged.RData"))

## MODELO COM AUTOCORRELACAO SEM EFEITO NUGGET
cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataST_sub","ctrl","sample_estimates"))

teste_exp  <- parLapply (cl, dataST_sub, function (i) {

	tryCatch(
		lme (ST ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=F)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
	}
)
stopCluster (cl)
save(teste_exp,file=here ("Output","small_ranges_output","res_step3_model_exp_ST_small_ranged.RData"))

### modelo sem autoc espacial
cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataST_sub","ctrl"))

no_aut<- parLapply (cl, dataST_sub, function (i) {

	tryCatch(
		lme (ST ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML"),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)
stopCluster (cl)
	
save(no_aut, file=here ("Output","small_ranges_output","res_step3_modelNoaut_ST_small_ranged.RData"))

###################################
#### LAST TRANSITION TIME

# obtain subsets of the complete datasets
dataLT_sub <- dataLT[sample_estimates]

# subtly jitter coordinates to avoid errors in GLMM
dataLT_sub <- lapply (dataLT_sub, function (i) {
  
  i$coordinates.V2_jitter <- jitter(i$coordinates.V2,1)
  ;
  i
})

# run

cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataLT_sub","ctrl"))

teste_exp_nugget <- parLapply (cl, dataLT_sub, function (i) {

	tryCatch(
		lme (LT ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=T)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
   }
)
stopCluster(cl)

## no convergence problem!
table (unlist (lapply (teste_exp_nugget, is.null)))
save (teste_exp_nugget, file=here ("Output","small_ranges_output","res_step3_teste_exp_nugget_LT_small_ranged.RData"))


## MODELO COM AUTOCORRELACAO SEM EFEITO NUGGET

cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataLT_sub","ctrl"))

teste_exp  <- parLapply (cl, dataLT_sub, function (i) {

	tryCatch(
		lme (LT ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML",
			correlation = corExp (form = ~coordinates.V1+coordinates.V2_jitter, nugget=F)),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
	}
)
stopCluster (cl)

table (unlist (lapply (teste_exp, is.null)))
save(teste_exp,file=here ("Output","small_ranges_output","res_step3_model_exp_LT_small_ranged.RData"))

### modelo sem autoc espacial
cl <- makeCluster(ncores) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(nlme))

# export your data and function
clusterExport(cl, c("dataLT_sub","ctrl"))

no_aut<- parLapply (cl, dataLT_sub, function (i) {

	tryCatch(
		lme (LT ~ position*tipo_hab+
			    neigh_area+ over_FOR+over_AB+ 
				mata_atl + andes_eco, 
			 data = i,control=ctrl,random = ~1|ecoreg, 
			 method= "REML"),#, 
            error = function(e) return ("NULL")#{print(e); print("retrying...")}
        )
    }
)
stopCluster(cl)

table (unlist (lapply (no_aut, is.null)))	
save(no_aut, file=here ("Output","small_ranges_output","res_step3_modelNoaut_LT_small_ranged.RData"))

# end