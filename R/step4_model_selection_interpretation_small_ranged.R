# ---------------------------------------------- #
# model selection and interpretation (Density plots)
# ---------------------------------------------- #

#load packages and function
source("R/packages.R")
source("R/other_functions.R")

# load aTR models
load (here ("Output","small_ranges_output","res_step3_model_exp_nugget_TR_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_model_exp_TR_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_modelNoaut_TR_small_ranged.RData"))

# remove nulls in nugget
teste_exp_nugget_sub <- teste_exp_nugget[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
teste_exp_sub<-teste_exp[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
no_aut_sub<-no_aut[which(unlist(lapply(teste_exp_nugget,is.character))==F)]

# rm nulls in exponential
teste_exp_subB<- teste_exp_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
teste_exp_nugget_subB <- teste_exp_nugget_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
no_aut_subB<-no_aut_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]

# comparing the exponential with nugget with the same without nugget effect
comp <- lapply(seq(1,length (no_aut_subB)), function (i) 
  anova (no_aut_subB [[i]], teste_exp_subB[[i]],
         test=F))

#save(comp,file="comparacao_noAut_exp_TR.RData")

# compare AIC
table (unlist(lapply (seq(1,length (comp)), function (i) 
  comp[[i]]$AIC[1] < comp[[i]]$AIC[2])))/length (comp)
table(lapply (seq (1,length (comp)), function (i) 
  comp [[i]]$AIC[1] - comp [[i]]$AIC[2]) > 2)

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,sd)

# compare nugget and wihtout nugget
compB <- lapply(seq(1,length (teste_exp_nugget_subB)), function (i) 
  anova (
  teste_exp_subB [[i]], teste_exp_nugget_subB[[i]],test=F))

# save(compB,file="comparacao_exp_nugget_rates.RData")
# load("comparacao_exp_nugget_rates.RData")
## quantas vezes sem auto teve AIC menor do que com autoc
# lower AIC
test1 <- table (unlist(lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC[1] < compB[[i]]$AIC[2])))
test1/length(compB)
# diff AIC
test2<-table (lapply (seq (1,length (compB)), function (i) 
  compB [[i]]$AIC[1] - compB [[i]]$AIC[2]) > 2)  
test2/length(compB)
## fica com modelo exponencial com efeito nugget

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC)),2,sd)

################## obter valores medios dos parametros

## VALORES DOS PARAMETROS
## fixed effects

coef_table <- lapply(teste_exp_subB, function (i) summary (i)$tTable)
## mean
list.mean <- Reduce("+", coef_table) / length(coef_table)
round(list.mean,3)
## squared mean
list.squared.mean <- Reduce("+", lapply(coef_table, "^", 2)) / length(coef_table)
## variance
list.variance <- list.squared.mean - list.mean^2
## standard deviation
list.sd <- sqrt(list.variance)
round(list.sd,3)

### random effect
# sd
rdn_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  as.matrix (VarCorr(teste_exp_subB[[i]])))
mean (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
sd (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
# residuals
mean (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))
sd (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))

## 
rdn_mean <- mean (as.numeric(unlist(rdn_table)))
rdn_sd <- sd(as.numeric(unlist(rdn_table )))
### spatial structure
space_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  coef(teste_exp_subB[[i]]$modelStruct$corStruct, unconstrained = F))
# range of spatial model
(space_mean_range <-  (unlist(lapply (seq (1,length (space_table)), function (i) 
  
  space_table [[i]][1]))
  ))
space_mean_range <- (space_mean_range[is.infinite(space_mean_range) == F])
mean(space_mean_range)
sd(space_mean_range[c(1:100)])

# in km scale
mean(space_mean_range)/1000
sd(space_mean_range[c(1:100)])/1000
# nugget effect
#(space_mean_nugget <- mean (unlist(lapply (seq (1,length (space_table)), function (i) space_table [[i]][2]))))
#(space_sd_nugget <- sd(unlist(lapply (seq (1,length (space_table)), function (i) space_table [[i]][2]))))

# ------------------------------------------------------
# data for density plots

todos_os_parametros <- lapply (as.list(seq(1,nrow(list.mean))), function (k) ### for each line in model output
  lapply (as.list (seq(1,length (coef_table))), function (i) ## and for each simulation (ACE x phylogeny)
    round (coef_table[[i]][k,1],3))) ### get the intercept (average tip based metric)

## average values of parameters (same as list.mean)
# mean
#lapply(
#  lapply(todos_os_parametros, unlist),
#  mean)
## sd
#lapply(
#  lapply(todos_os_parametros, unlist),
#  sd)
## do all overlap zero?
#data.frame (q5=unlist(lapply(
#  lapply(todos_os_parametros, unlist),
#  quantile, .05)),
#  q95=unlist(lapply(
#    lapply(todos_os_parametros, unlist),
#    quantile, .95)))
#
## cores para o histograma; o ultimo numero eh transparencia
#red: rgb(1,0,0,1)
#green: rgb(0,1,0,1)
#blue: rgb(0,0,1,1)
#gray: rgb(0.5,0.5,0.5,1)
#yellow: rgb(1,0.9,0,1)
#purple: rgb(1,0,1,1)
#orange: rgb(1,0.5,0,1)

## params
variaveis <- c("Intercept",
               "Ecotone",
               "Forest", 
               "Neighbor area",
               "Forest neighbors",	
               "Open-habitat neighbors",
               "Atlantic Rainforest",
               "Andes",
               "Ecotone-Forest")

lapply (seq(2,9), function (i) {
  
  df.plot <- rbind (data.frame (Value = unlist(todos_os_parametros[[1]]),
                                Parameter  = variaveis[1],
                                grp.mean = mean(unlist(todos_os_parametros[[1]]))),
                    data.frame (Value=unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])),
                                Parameter = "Ecotone",
                                grp.mean = mean(unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])))))
  
  #### plot density
  
  a <- ggplot(df.plot, aes(x=Value, color=Parameter, fill=Parameter)) +
    geom_density(size=1,alpha=0.1) +
    scale_fill_manual(values=c("Intercept"="gray", "Ecotone"="black")) + 
    scale_colour_manual(values=c("Intercept"="gray", "Ecotone"="black"))+
    theme_classic()  
  
  b <- a+ theme (legend.title = element_blank(),
                 #legend.position = "none",
                 legend.position = c(.95, .95),
                 legend.justification = c("right", "top"),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 plot.margin=unit(c(.2,1,.1,1),"cm")) +
    xlab ("Assemblage Transition Rates (aTR)") + 
    ylab("Density") + 
    scale_x_continuous(breaks = seq(0,0.68,0.15),
                       limits = c(0,0.68)) 
  
  ## boxplot
  
  p <- ggboxplot(df.plot, x = "Parameter", y = "Value",
                 color = "Parameter", palette =c("gray", "black"),
                 add = "jitter", add.params = list(alpha=0.05),size=0.5)
  
  p1 <- p + theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank(),
                  plot.margin=unit(c(.2,1,.1,1),"cm")) +
    scale_y_continuous(limits=c(0, 0.68)) + coord_flip() 
  
  
  # plot and save each density plot
  pdf (here ("Output","Figures_small_ranges", 
             paste ("compTR",variaveis[i],".pdf",sep="")), 
       width=4,height=3,family="serif")
  
  grid.arrange(p1,b,ncol=4,nrow = 4,
               layout_matrix = rbind (c(1,1,1,1),
                                      c(2,2,2,2),
                                      c(2,2,2,2),
                                      c(2,2,2,2)));
  
 
  
  dev.off()
  
})


#############################################
######## STASIS TIME

# load aTR models
load (here ("Output","small_ranges_output","res_step3_teste_exp_nugget_ST_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_model_exp_ST_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_modelNoaut_ST_small_ranged.RData"))

# remove nulls in nugget
teste_exp_nugget_sub <- teste_exp_nugget[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
teste_exp_sub<-teste_exp[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
no_aut_sub<-no_aut[which(unlist(lapply(teste_exp_nugget,is.character))==F)]

# rm nulls in exponential
teste_exp_subB<- teste_exp_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
teste_exp_nugget_subB <- teste_exp_nugget_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
no_aut_subB<-no_aut_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]

# comparing the exponential with nugget with the same without nugget effect
comp <- lapply(seq(1,length (no_aut_subB)), function (i) 
  anova (no_aut_subB [[i]], teste_exp_subB[[i]],
         test=F))

#save(comp,file="comparacao_noAut_exp_TR.RData")

# compare AIC
test1<- table (unlist(lapply (seq(1,length (comp)), function (i) 
  comp[[i]]$AIC[1] < comp[[i]]$AIC[2])))
test1/length(comp)
test2<- table(lapply (seq (1,length (comp)), function (i) 
  comp [[i]]$AIC[1] - comp [[i]]$AIC[2]) > 2)
test2/length(comp)

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,sd)

# compare nugget and wihtout nugget
compB <- lapply(seq(1,length (teste_exp_nugget_subB)), function (i) anova (
  teste_exp_subB [[i]], teste_exp_nugget_subB[[i]],test=F))

# save(compB,file="comparacao_exp_nugget_rates.RData")
# load("comparacao_exp_nugget_rates.RData")
## quantas vezes sem auto teve AIC menor do que com autoc
test1<-table (unlist(lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC[1] < compB[[i]]$AIC[2])))
test1/length(compB)

test2<-table (lapply (seq (1,length (compB)), function (i) 
  compB [[i]]$AIC[1] - compB [[i]]$AIC[2]) > 2)
test2/length(compB)

## fica com modelo exponencial com efeito nugget

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) 
    compB[[i]]$AIC)),2,sd)

################## obter valores medios dos parametros

## VALORES DOS PARAMETROS
## fixed effects

coef_table <- lapply(teste_exp_subB, function (i) summary (i)$tTable)
## mean
list.mean <- Reduce("+", coef_table) / length(coef_table)
round(list.mean,3)
## squared mean
list.squared.mean <- Reduce("+", lapply(coef_table, "^", 2)) / length(coef_table)
## variance
list.variance <- list.squared.mean - list.mean^2
## standard deviation
list.sd <- sqrt(list.variance)
round(list.sd,3)

### random effect
# sd
rdn_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  as.matrix (VarCorr(teste_exp_subB[[i]])))
mean (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
sd (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
# residuals
mean (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))
sd (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))

## 
rdn_mean <- mean (as.numeric(unlist(rdn_table)))
rdn_sd <- sd(as.numeric(unlist(rdn_table )))
### spatial structure
### spatial structure
space_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  coef(teste_exp_subB[[i]]$modelStruct$corStruct, unconstrained = F))
# range of spatial model
(space_mean_range <-  (unlist(lapply (seq (1,length (space_table)), function (i) 
  
  space_table [[i]][1]))
))
space_mean_range <- (space_mean_range[is.infinite(space_mean_range) == F])
mean(space_mean_range)
sd(space_mean_range[c(1:100)])

# in km scale
mean(space_mean_range)/1000
sd(space_mean_range[c(1:100)])/1000
# nugget effect

# ------------------------------------------------------
# data for density plots

todos_os_parametros <- lapply (as.list(seq(1,nrow(list.mean))), function (k) ### for each line in model output
  lapply (as.list (seq(1,length (coef_table))), function (i) ## and for each simulation (ACE x phylogeny)
    round (coef_table[[i]][k,1],3))) ### get the intercept (average tip based metric)

## params
variaveis <- c("Intercept",
               "Ecotone",
               "Forest", 
               "Neighbor area",
               "Forest neighbors",	
               "Open-habitat neighbors",
               "Atlantic Rainforest",
               "Andes",
               "Ecotone-Forest")

lapply (seq(2,9), function (i) {
  
  df.plot <- rbind (data.frame (Value = unlist(todos_os_parametros[[1]]),
                                Parameter  = variaveis[1],
                                grp.mean = median(unlist(todos_os_parametros[[1]]))),
                    data.frame (Value=unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])),
                                Parameter = "Ecotone",
                                grp.mean = median(unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])))))
  
  #### plot density
  
  a <- ggplot(df.plot, aes(x=Value, color=Parameter, fill=Parameter)) +
    geom_density(size=1,alpha=0.1) +
    scale_fill_manual(values=c("Intercept"="gray", "Ecotone"="black")) + 
    scale_colour_manual(values=c("Intercept"="gray", "Ecotone"="black"))+
    theme_classic()  
  
  b <- a+ theme (legend.title = element_blank(),
                 #legend.position = "none",
                 legend.position = c(.95, .95),
                 legend.justification = c("right", "top"),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 plot.margin=unit(c(.2,1,.1,1),"cm")) +
    xlab ("Assemblage Stasis Time (aST)") + 
    ylab("Density") + 
    scale_x_continuous(breaks = seq(0,5,0.5),
                       limits = c(0,5)) 
  
  ## boxplot
  
  p <- ggboxplot(df.plot, x = "Parameter", y = "Value",
                 color = "Parameter", palette =c("gray", "black"),
                 add = "jitter", add.params = list(alpha=0.05),size=0.5)
  
  p1 <- p + theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank(),
                  plot.margin=unit(c(.2,1,.1,1),"cm")) +
    scale_y_continuous(limits=c(0, 5)) + coord_flip() 
  
  
  pdf (here ("Output","Figures_small_ranges",paste ("compST",variaveis[i],".pdf",sep="")), 
       width=5,height=4,family="serif")
  grid.arrange(p1,b,ncol=4,nrow = 4,
               layout_matrix = rbind (c(1,1,1,1),
                                      c(2,2,2,2),
                                      c(2,2,2,2),
                                      c(2,2,2,2)))
  
  dev.off()
  
})

# ------------------------------------------------------------ #
# LT

# load aTR models
load (here ("Output","small_ranges_output","res_step3_teste_exp_nugget_LT_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_model_exp_LT_small_ranged.RData"))
load (here ("Output","small_ranges_output","res_step3_modelNoaut_LT_small_ranged.RData"))

# remove nulls in nugget
teste_exp_nugget_sub <- teste_exp_nugget[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
teste_exp_sub<-teste_exp[which(unlist(lapply(teste_exp_nugget,is.character))==F)]
no_aut_sub<-no_aut[which(unlist(lapply(teste_exp_nugget,is.character))==F)]

# rm nulls in exponential
teste_exp_subB<- teste_exp_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
teste_exp_nugget_subB <- teste_exp_nugget_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]
no_aut_subB<-no_aut_sub[which(unlist(lapply(teste_exp_sub,is.character))==F)]

# comparing the exponential with nugget with the same without nugget effect
comp <- lapply(seq(1,length (no_aut_subB)), function (i) 
  anova (no_aut_subB [[i]], teste_exp_subB[[i]],
         test=F))

#save(comp,file="comparacao_noAut_exp_TR.RData")

# compare AIC
test1<-table (unlist(lapply (seq(1,length (comp)), function (i) 
  comp[[i]]$AIC[1] < comp[[i]]$AIC[2])))
test1/length(comp)

test2<-table (lapply (seq (1,length (comp)), function (i) 
  comp [[i]]$AIC[1] - comp [[i]]$AIC[2]) > 2)
test2/length(comp)

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (teste_exp_subB)), function (i) 
  comp[[i]]$AIC)),2,sd)

# compare nugget and wihtout nugget
compB <- lapply(seq(1,length (teste_exp_nugget_subB)), function (i) anova (
  teste_exp_subB [[i]], teste_exp_nugget_subB[[i]],test=F))

# save(compB,file="comparacao_exp_nugget_rates.RData")
# load("comparacao_exp_nugget_rates.RData")
## quantas vezes sem auto teve AIC menor do que com autoc
test1<-table (unlist(lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC[1] < compB[[i]]$AIC[2])))
test1/length(compB)

test2<-table (lapply (seq (1,length (compB)), function (i) 
  compB [[i]]$AIC[1] - compB [[i]]$AIC[2]) > 2)
test2/length(compB)
## fica com modelo exponencial com efeito nugget

## valores medios dos parametros
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) 
  compB[[i]]$AIC)),2,mean)
## desvio padrao
apply (do.call (rbind.data.frame, lapply (seq(1,length (compB)), function (i) compB[[i]]$AIC)),2,sd)

################## obter valores medios dos parametros

## VALORES DOS PARAMETROS
## fixed effects

coef_table <- lapply(teste_exp_subB, function (i) summary (i)$tTable)
## mean
list.mean <- Reduce("+", coef_table) / length(coef_table)
round(list.mean,3)
## squared mean
list.squared.mean <- Reduce("+", lapply(coef_table, "^", 2)) / length(coef_table)
## variance
list.variance <- list.squared.mean - list.mean^2
## standard deviation
list.sd <- sqrt(list.variance)
round(list.sd,3)

### random effect
# sd
rdn_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  as.matrix (VarCorr(teste_exp_subB[[i]])))
mean (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
sd (unlist(int_sd<- lapply (rdn_table, function (i) as.numeric (i[1,2]))))
# residuals
mean (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))
sd (unlist(res_sd<- lapply (rdn_table, function (i) as.numeric (i[2,2]))))

## 
rdn_mean <- mean (as.numeric(unlist(rdn_table)))
rdn_sd <- sd(as.numeric(unlist(rdn_table )))
### spatial structure
### spatial structure
space_table <- lapply(seq(1,length (teste_exp_subB)), function (i) 
  coef(teste_exp_subB[[i]]$modelStruct$corStruct, unconstrained = F))
# range of spatial model
(space_mean_range <-  (unlist(lapply (seq (1,length (space_table)), function (i) 
  
  space_table [[i]][1]))
))
space_mean_range <- (space_mean_range[is.infinite(space_mean_range) == F])
mean(space_mean_range)
sd(space_mean_range[c(1:100)])

# in km scale
mean(space_mean_range)/1000
sd(space_mean_range[1000:1020],na.rm=T)/1000
# nugget effect
# ------------------------------------------------------
# data for density plots

todos_os_parametros <- lapply (as.list(seq(1,nrow(list.mean))), function (k) ### for each line in model output
  lapply (as.list (seq(1,length (coef_table))), function (i) ## and for each simulation (ACE x phylogeny)
    round (coef_table[[i]][k,1],3))) ### get the intercept (average tip based metric)

## params
variaveis <- c("Intercept",
               "Ecotone",
               "Forest", 
               "Neighbor area",
               "Forest neighbors",	
               "Open-habitat neighbors",
               "Atlantic Rainforest",
               "Andes",
               "Ecotone-Forest")

lapply (seq(2,9), function (i) {
  
  df.plot <- rbind (data.frame (Value = unlist(todos_os_parametros[[1]]),
                                Parameter  = variaveis[1],
                                grp.mean = median(unlist(todos_os_parametros[[1]]))),
                    data.frame (Value=unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])),
                                Parameter = "Ecotone",
                                grp.mean = median(unlist(todos_os_parametros[[1]])+ (unlist(todos_os_parametros[[i]])))))
  
  #### plot density
  
  a <- ggplot(df.plot, aes(x=Value, color=Parameter, fill=Parameter)) +
    geom_density(size=1,alpha=0.1) +
    scale_fill_manual(values=c("Intercept"="gray", "Ecotone"="black")) + 
    scale_colour_manual(values=c("Intercept"="gray", "Ecotone"="black"))+
    theme_classic()  
  
  b <- a+ theme (legend.title = element_blank(),
                 #legend.position = "none",
                 legend.position = c(.95, .95),
                 legend.justification = c("right", "top"),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 plot.margin=unit(c(.2,1,.1,1),"cm")) +
    xlab ("Assemblage Last Transition Time (aLT)") + 
    ylab("Density") + 
    scale_x_continuous(breaks = seq(0,13,1),
                       limits = c(0,13)) 
  
  ## boxplot
  
  p <- ggboxplot(df.plot, x = "Parameter", y = "Value",
                 color = "Parameter", palette =c("gray", "black"),
                 add = "jitter", add.params = list(alpha=0.05),
                 size=0.5)
  
  p1 <- p + theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank(),
                  plot.margin=unit(c(.2,1,.1,1),"cm")) +
    scale_y_continuous(limits=c(0, 13)) + coord_flip() 
  
  
  pdf (here ("Output","Figures_small_ranges",paste ("compLT",variaveis[i],".pdf",sep="")), 
       width=5,height=4,family="serif")
  grid.arrange(p1,b,ncol=4,nrow = 4,
               layout_matrix = rbind (c(1,1,1,1),
                                      c(2,2,2,2),
                                      c(2,2,2,2),
                                      c(2,2,2,2)))
  
  dev.off()
  
})

# end
