
# -----------------------------------------
# Step 4, mapping and supporting information
# -----------------------------------------

#load packages and functions
source("R/packages.R")
source("R/other_functions.R")

# load ecoregions
ecor<- readOGR(dsn=here("Data","Environment","Ecorregions", "terr-ecoregions-TNC"),
               layer="tnc_terr_ecoregions")
ecor <- ecor [which (ecor@data$WWF_REALM == "NT"), ]
ecor <- ecor [order (ecor@data$ECO_CODE, decreasing=F), ]
# transform into lambert projection
SA_ecor_EA <- spTransform(ecor,
                          CRS = "+proj=laea +x_0=0 +y_0=0 +lon_0=-56 +lat_0=-15 +datum=WGS84" )

## load modeling data
load (here ("Output","res_step2_data_for_GLMM.RData"))

# load PD data
load(here ("Output","res_step2_phylogenetic_div.RData"))

# relationship between area and TR, per position in the ecoregion
plot(
	gArea (ecor,byid=T) [match (dataTR[[1]]$ecoreg,
	                            SA_ecor_EA@data$ECO_CODE)],
	dataTR[[1]]$TR,
	col = dataTR[[1]]$position)

### maps
# subset of ecoregions we could have data
SA_ecor3 <- SA_ecor_EA [which(SA_ecor_EA$ECO_CODE %in% dataTR[[1]]$ecoreg),]

# --------------------------------------------------------
# data for mapping
data_plot_TR <- cbind (dataTR [[1]][,-which(colnames(dataTR [[1]]) == "TR")],
                       PD = apply (sapply (pd_list,"[[","PD"),1,mean),
                       TRm = apply(sapply (dataTR,"[[","TR"),1,mean),
	                     TRsd = apply(sapply (dataTR,"[[","TR"),1,sd),
                       STm = apply(sapply (dataST,"[[","ST"),1,mean),
                       STsd = apply(sapply (dataST,"[[","ST"),1,sd),
                       LTm = apply(sapply (dataLT,"[[","LT"),1,mean),
                       LTsd = apply(sapply (dataLT,"[[","LT"),1,sd)
	        )

# correlation between tip-based metrics
cor (data_plot_TR[,c("TRm","STm","LTm")])
# significance
cor.test (data_plot_TR[,"TRm"], data_plot_TR[,"STm"])
cor.test (data_plot_TR[,"TRm"], data_plot_TR[,"LTm"])
cor.test (data_plot_TR[,"STm"], data_plot_TR[,"LTm"])


# ------------------------------------------------
# relationship between species richness and tip based metrics
# TR
m1TR <- gls (TRm ~ riqueza, data = data_plot_TR,
	correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))
# ST
m1ST <- gls (STm ~ riqueza, data = data_plot_TR,
             correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))
# LT
m1LT <- gls (LTm ~ riqueza, data = data_plot_TR,
             correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))

# ----------------------------------------------------
# relationship between PD and tip based metrics	
# TR
m1TR_PD <- gls (TRm ~ PD, data = data_plot_TR,
	correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))
# ST
m1ST_PD <- gls (STm ~ PD, data = data_plot_TR,
                correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))
# ST
m1LT_PD <- gls (LTm ~ PD, data = data_plot_TR,
                correlation = corExp (form = ~coordinates.V1+jitter(coordinates.V2,1), nugget=T))

# bind the residuals in the DF  to map

# bind richness TR residuals in DF (Only for Species richness)
data_plot_TR <- cbind (data_plot_TR, 
                       TRres = m1TR$residuals,
                       STres = m1ST$residuals,
                       LTres = m1LT$residuals)

#----------------------------------------------------------------

# south america map for all maps/ panels
south_america_map <- ggplot() + 
  geom_polygon (data=SA_ecor3, aes(x=long, y=lat, group=group),size = 0.3, 
                fill="gray85", colour="gray55",alpha=0.5) +
  coord_fixed (ratio = 1) +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())


# map TR (FIg 5A, B)
b.mean <- south_america_map + geom_point (data=data_plot_TR, 
                              aes (x=coordinates.V1, y=coordinates.V2,
                                   colour=TRm),
                              size = 2,
                              alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Average aTR",
                         limits = c(min(data_plot_TR$TRm),
                                    max(data_plot_TR$TRm))
  )
# label
b.mean <- b.mean + annotate(geom="text", x=2333214, y=2999892 , label="A",
                          color="black",size=8) 
# leGend pos
b.mean <- b.mean+theme(legend.position = c(0.9, 0.25), 
                       legend.key.width=unit(0.55,"cm"),
                       legend.key.height=unit(0.45,"cm"),
                       legend.key.size = unit(0.4,"cm"),
                       legend.text = element_text(color = "black", size = 7),
                       strip.text = element_text(size=5),
                       legend.title = element_text (size=7))
#  SD of TR

b.sd <- south_america_map + geom_point (data=data_plot_TR, 
                                        aes (x=coordinates.V1, y=coordinates.V2,colour=TRsd),
                                        size = 2,
                                        alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="SD aTR",
                         limits = c(min(data_plot_TR$TRsd),
                                    max(data_plot_TR$TRsd))
                         )

# label
b.sd <- b.sd + annotate(geom="text", x=2333214, y=2999892 , label="B",
                            color="black",size=8) 
# leGend pos
b.sd <- b.sd  + theme(legend.position = c(0.9, 0.25), 
                      legend.key.width=unit(0.55,"cm"),
                      legend.key.height=unit(0.45,"cm"),
                      legend.key.size = unit(0.4,"cm"),
                      legend.text = element_text(color = "black", size = 7),
                      strip.text = element_text(size=5),
                      legend.title = element_text (size=7))

#---------------------------- 
# map statis time (fig 5C,D)

b.mean.ST <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                         y=coordinates.V2,
                                                         colour=STm),
                                 size = 2,
                                 alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Average aST",
                         limits = c(min(data_plot_TR$STm),
                                    max(data_plot_TR$STm))
  )
# label
b.mean.ST <- b.mean.ST + annotate(geom="text", x=2333214, y=2999892 , label="C",
                            color="black",size=8) 
# leGend pos
b.mean.ST <- b.mean.ST +theme(legend.position = c(0.9, 0.25), 
                              legend.key.width=unit(0.55,"cm"),
                              legend.key.height=unit(0.45,"cm"),
                              legend.key.size = unit(0.4,"cm"),
                              legend.text = element_text(color = "black", size = 7),
                              strip.text = element_text(size=5),
                              legend.title = element_text (size=7))

# map of SD

b.sd.ST <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                                   y=coordinates.V2,
                                                                   colour=STsd),
                                           size = 2,
                                           alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="SD aST",
                         limits = c(min(data_plot_TR$STsd),
                                    max(data_plot_TR$STsd))
  )
# label
b.sd.ST <- b.sd.ST + annotate(geom="text", x=2333214, y=2999892 , label="D",
                                  color="black",size=8) 
# leGend pos
b.sd.ST <- b.sd.ST +theme(legend.position = c(0.9, 0.25), 
                          legend.key.width=unit(0.55,"cm"),
                          legend.key.height=unit(0.45,"cm"),
                          legend.key.size = unit(0.4,"cm"),
                          legend.text = element_text(color = "black", size = 7),
                          strip.text = element_text(size=5),
                          legend.title = element_text (size=7))

# ------------------------------------------
## map of LT (figs 5E,F)

b.mean.LT <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                                     y=coordinates.V2,
                                                                     colour=LTm),
                                             size = 2,
                                             alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Average aLT",
                         limits = c(min(data_plot_TR$LTm),
                                    max(data_plot_TR$LTm))
  )
# label
b.mean.LT <- b.mean.LT + annotate(geom="text", x=2333214, y=2999892 , label="E",
                                  color="black",size=8) 
# leGend pos
b.mean.LT <- b.mean.LT +theme(legend.position = c(0.9, 0.25), 
                              legend.key.width=unit(0.55,"cm"),
                              legend.key.height=unit(0.45,"cm"),
                              legend.key.size = unit(0.4,"cm"),
                              legend.text = element_text(color = "black", size = 7),
                              strip.text = element_text(size=5),
                              legend.title = element_text (size=7))

# map of SD

b.sd.LT <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                                   y=coordinates.V2,
                                                                   colour=LTsd),
                                           size = 2,
                                           alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="SD aST",
                         limits = c(min(data_plot_TR$LTsd),
                                    max(data_plot_TR$LTsd))
  )
# label
b.sd.LT <- b.sd.LT + annotate(geom="text", x=2333214, y=2999892 , label="F",
                              color="black",size=8) 
# leGend pos
b.sd.LT <- b.sd.LT +theme(legend.position = c(0.9, 0.25), 
                          legend.key.width=unit(0.55,"cm"),
                          legend.key.height=unit(0.45,"cm"),
                          legend.key.size = unit(0.4,"cm"),
                          legend.text = element_text(color = "black", size = 7),
                          strip.text = element_text(size=5),
                          legend.title = element_text (size=7))

# arrange maps for the main text

pdf (here ("Output","Figures","Fig5.pdf"), 
     width=6,height=11,family="serif")
# arrange
grid.arrange(b.mean, b.sd,                              
             b.mean.ST,b.sd.ST,
             b.mean.LT, b.sd.LT,
             ncol = 12, nrow =12, 
             layout_matrix = rbind(c(rep(1,6), rep(2,6)),
                                   c(rep(1,6), rep(2,6)),
                                   c(rep(1,6), rep(2,6)),
                                   c(rep(1,6), rep(2,6)),
                                   c(rep(3,6), rep(4,6)),
                                   c(rep(3,6), rep(4,6)),
                                   c(rep(3,6), rep(4,6)),
                                   c(rep(3,6), rep(4,6)),
                                   c(rep(5,6), rep(6,6)),
                                   c(rep(5,6), rep(6,6)),
                                   c(rep(5,6), rep(6,6)),
                                   c(rep(5,6), rep(6,6))
             ))
dev.off()


# -------------------------------------------------------
# SUPPORTING INFORMATION MAPS
# richness and tip based metrics

# map of richness
b.riq <- south_america_map+geom_point (data=dataTR[[1]], 
                                       aes (x=coordinates.V1, y=coordinates.V2,colour=riqueza),
                                       size = 2,
                                       alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Richness",
                         limits = c(min(data_plot_TR$riqueza),
                                    max(data_plot_TR$riqueza))
  )

b.riq <- b.riq + annotate(geom="text", x=2333214, y=2999892 , label="A",
                          color="black",size=8) 

b.riq <- b.riq +theme(legend.position = c(0.9, 0.25), 
                      legend.key.width=unit(0.55,"cm"),
                      legend.key.height=unit(0.45,"cm"),
                      legend.key.size = unit(0.4,"cm"),
                      legend.text = element_text(color = "black", size = 7),
                      strip.text = element_text(size=5),
                      legend.title = element_text (size=7))

## RIchness - aTR
b.res.TR <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                     y=coordinates.V2,
                                                     colour=TRres),
                                            size = 2,
                                            alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Res(aTR~SR)",
                         limits = c(min(data_plot_TR$TRres),
                                    max(data_plot_TR$TRres))
  )
# label
b.res.TR <- b.res.TR + annotate(geom="text", x=2333214, y=2999892 , label="B",
                              color="black",size=8) 
# leGend pos
b.res.TR <- b.res.TR +theme(legend.position = c(0.9, 0.25), 
                            legend.key.width=unit(0.55,"cm"),
                            legend.key.height=unit(0.45,"cm"),
                            legend.key.size = unit(0.4,"cm"),
                            legend.text = element_text(color = "black", size = 7),
                            strip.text = element_text(size=5),
                            legend.title = element_text (size=7))
## RIchness - aST
b.res.ST <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                                    y=coordinates.V2,
                                                                    colour=STres),
                                            size = 2,
                                            alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Res(aST~SR)",
                         limits = c(min(data_plot_TR$STres),
                                    max(data_plot_TR$STres))
  )
# label
b.res.ST <- b.res.ST + annotate(geom="text", x=2333214, y=2999892 , label="C",
                                color="black",size=8) 
# leGend pos
b.res.ST <- b.res.ST +theme(legend.position = c(0.9, 0.25), 
                            legend.key.width=unit(0.55,"cm"),
                            legend.key.height=unit(0.45,"cm"),
                            legend.key.size = unit(0.4,"cm"),
                            legend.text = element_text(color = "black", size = 7),
                            strip.text = element_text(size=5),
                            legend.title = element_text (size=7))

# res (aLT ~ richness)
b.res.LT <- south_america_map + geom_point (data=data_plot_TR, aes (x=coordinates.V1, 
                                                                    y=coordinates.V2,
                                                                    colour=LTres),
                                            size = 2,
                                            alpha=0.5) + 
  scale_colour_viridis_c(option="magma",
                         direction=-1,
                         #midpoint = 12
                         name="Res(aLT~SR)",
                         limits = c(min(data_plot_TR$LTres),
                                    max(data_plot_TR$LTres))
  )
# label
b.res.LT <- b.res.LT + annotate(geom="text", x=2333214, y=2999892 , label="D",
                                color="black",size=8) 
# leGend pos
b.res.LT <- b.res.LT +theme(legend.position = c(0.9, 0.25), 
                            legend.key.width=unit(0.55,"cm"),
                            legend.key.height=unit(0.45,"cm"),
                            legend.key.size = unit(0.4,"cm"),
                            legend.text = element_text(color = "black", size = 7),
                            strip.text = element_text(size=5),
                            legend.title = element_text (size=7))


# arrange maps of supp infor
pdf (here ("Output","Figures","FigS12.pdf"), 
           width=6,height=8,family="serif")
grid.arrange(b.riq, b.res.TR,                              
             b.res.ST,b.res.LT, 
             ncol = 4, nrow = 8, 
             layout_matrix = rbind(c(1,1,2,2),
                                   c(1,1,2,2),
                                   c(1,1,2,2),
                                   c(1,1,2,2),
                                   c(3,3,4,4),
                                   c(3,3,4,4),
                                   c(3,3,4,4),
                                   c(3,3,4,4)
             ))
dev.off()


pdf (file = here ("Output","Figures","rich_PD_tipbasedindex.pdf"), width = 8, height=6)
par(mfrow=c(2,3),mar=c(4,4,4,4))
## TR
with (data_plot_TR,plot (riqueza,TRm,pch=19, 
                         xlab= "", ylab="Assemblage transition rates (aTR)"))
abline (lm(TRm~riqueza,data_plot_TR,),lwd=2,col="gray30")
summary(m1TR)
text (30,0.4,"p<0.05")
## ST
with (data_plot_TR,plot (riqueza,STm,pch=19, 
                         xlab= "Richness", ylab="Assemblage stasis time (aST)"))
with (data_plot_TR,abline (m1ST,lwd=2,col="gray30"))
summary (m1ST)
text (30,4.15,"p<0.001")
## LT
with (data_plot_TR,plot (riqueza,LTm,pch=19, 
                         ylab="Assemblage last transition time (aLT)"))
with (data_plot_TR,abline (m1LT,lwd=2,col="gray30"))
summary(m1LT)
text (30,8.5,"p<0.01")

## relationship between PD and tip-based metrics

## TR
with (data_plot_TR,plot (PD,TRm,pch=19, 
                         xlab= "", ylab="Assemblage transition rates (aTR)"))
with (data_plot_TR,abline (m1TR_PD,lwd=2,col="gray30"))
summary (m1TR_PD)
text (30,0.4,"p>0.05")
## ST
with (data_plot_TR,plot (PD,STm,pch=19, 
                         xlab= "Phylogenetic Diversity", ylab="Assemblage stasis time (aST)"))
with (data_plot_TR,abline (m1ST_PD,lwd=2,col="gray30"))
summary (m1ST_PD)
text (30,4.15,"p>0.05")
## LT
with (data_plot_TR,plot (PD,LTm,pch=19, 
                         xlab= "", ylab="Assemblage last transition time (aLT)"))
with (data_plot_TR,abline (m1LT_PD,lwd=2,col="gray30"))
summary(m1LT_PD)
text (30,8.5,"p<0.001")

dev.off()

# ------------------------------------------------------ #
# supporting information
# relationship betwen PD and position
anova(lm(PD ~ position, data = data_plot_TR))
PD_TR <- ggplot (data=data_plot_TR, aes(x=position,y=PD)) +
  geom_boxplot(fill="gray75") + theme_classic() +
  xlab("Point position")+
  ylab("Phylogenetic diversity (PD)")+
  annotate("text", x=1.5,y=31,
           label = "F=0.16\np=0.69")

ggsave(filename = here ("Output","Figures","SuppInfo_PD_position.pdf"),
       family="serif",
       width = 3,
       height= 4,dpi = 300)



# -------------------------------- #
# run function to check where uncertainty varies more (phylogeny or ancestral reconstruction)
# the basis of uncertainty
# -------------------------------- #
n_phylo<-100 # upham et al. 2019
# dataframe with combinations of phylogenies and iterations of ancestral reconstruction
# 100 phylogenies x 1000 reconstructions per phylogeny
df_vals <- data.frame (s1 = seq (1,length(dataTR),n_phylo),
                       s2 = seq (0,length(dataTR),n_phylo)[-1])# minus the zero

# number of randomizations
niter <- 100

# transition rates
resTR <- lapply (seq(1,niter ), function (k)
		do.call (rbind,lapply (seq(1,nrow(df_vals)), function (i)

			fc_tr (df_vals [i,1],df_vals [i,2],10,
				data=dataTR)
		))
	)

# sd within and between phylogenies (which one generally is higher)
# the proportion in which sd_within was lower than sd between
table (
	do.call(rbind,resTR)[,1]<do.call(rbind,resTR)[,2]
)/length(dataTR)

# stasis time 
resST <- lapply (seq(1,niter ), function (k)
		do.call (rbind,lapply (seq(1,nrow(df_vals)), function (i)

			fc_tr (df_vals [i,1],df_vals [i,2],10,
				data=dataST)
		))
	)
# proportion
table (
	do.call(rbind,resST)[,1]<do.call(rbind,resST)[,2]
)/length(dataST)

# last trasition time
resLT <- lapply (seq(1,niter ), function (k)
		do.call (rbind,lapply (seq(1,nrow(df_vals)), function (i)

			fc_tr (df_vals [i,1],df_vals [i,2],10,
				data=dataLT)
		))
	)
# proportion
table (
	do.call(rbind,resLT)[,1]<do.call(rbind,resLT)[,2]
)/length(dataLT)


## correlacao entre indices

cor (data.frame(data_plot_TR$riqueza,
                     data_plot_TR$PD,
                     data_plot_TR$TRm,
                     data_plot_TR$STm,
                     data_plot_TR$LTm))

# averages
mean(data_plot_TR$TRm)
sd(data_plot_TR$TR)

mean(data_plot_ST$ST)
sd(data_plot_ST$ST)

mean(data_plot_LT$LT)
sd(data_plot_LT$LT)


# end