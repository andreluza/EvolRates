# ----------------------------------------------- #
# step 2- up to assemblage scale
# organizing evolution data and ecoregion variables to GLMM analysis 
# ------------------------------------------------ #

#load packages and function
source("R/packages.R")
source("R/other_functions.R")

# load ecoregions
ecor<- readOGR(dsn=here("Data","Environment","Ecorregions", "terr-ecoregions-TNC"),
                        layer="tnc_terr_ecoregions")
ecor <- ecor [which (ecor@data$WWF_REALM == "NT"), ] # subsetting to get Neotropical ecoregions
ecor <- ecor [order (ecor@data$ECO_CODE, decreasing=F), ] # ordering
#proj4string(ecor) <- crs(ecor) # setting projection

# load environment data (ecoregion properties)
env<-read.csv (here ("Data","Environment","env.csv"),h=T,sep=";",row.names=1)
env <- env [order (rownames (env), decreasing=F),]
env <- env [which (env$WWF_REALM == "NT"),]

# bind env into ecoregion dataframe
ecor@data <- cbind (ecor@data, Forest=env$Forest)
ecor@data <- cbind (ecor@data, count=seq(1:length(rownames(ecor@data)))) # count of ecoregions
# plot(gBuffer (gCentroid (ecor, byid=T), byid=T, width=1))
# sum(gIsValid(gBuffer (gCentroid (ecor, byid=T), byid=T, width=1), byid=TRUE)==FALSE)

# ---------------------------------------------------------
# load Sigmodontinae distribution data (shapefiles from Patton et al. 2015)
rodent_data <- readOGR(dsn=here ("Data","Occurrence"),layer="Sigmodontinae SA")
crs (rodent_data) <- crs (ecor) # reference coordinate
rodent_data@data$binomial <- gsub (" ","_", rodent_data@data$binomial) # edit species names to match phylogeny and traits

# binding different polygons from the same species (disjunct distribution)
a<-unionSpatialPolygons(rodent_data, 
                        IDs=rodent_data@data$binomial, # bind polygons by species ID
                        threshold=NULL, 
                        avoidGEOS=FALSE, 
                        avoidUnaryUnion=FALSE) 
nc90_df <- as(rodent_data, "data.frame")[!duplicated(rodent_data@data$binomial), ]# identify and remove duplicates
row.names(nc90_df) <- paste(nc90_df$binomial)# names
rodent_data <- SpatialPolygonsDataFrame(a, nc90_df) # back into spatialPolygons
# check
# plot(rodent_data,
#     col = rgb(0.5,0.1,0.1,alpha=0.05),
#     border=rgb(0.1,0.1,0.1,alpha=0.2))# [grep ("Akodon", rodent_data@data$binomial),])

## ------------------------------------------------------------------------ #
# load evolutionary ancestral reconstruction (trait) per species
## ------------------------------------------------------------------------ #

evol_rates <- readRDS (here ("Output","res_step1_transitions_estimates_ages.rda"))

## average of tip-based metrics
# number of nodes per species
nnodes <- (lapply (evol_rates, function (i) sapply(i,"[[", "diet_total_nodes")))
nnodes <- do.call (cbind,nnodes)
mean(apply (nnodes,1,mean))
sd(apply (nnodes,1,mean))
# number of diet transitions per species
prop_trans <-  (lapply (evol_rates, function (i) sapply(i,"[[", "diet_transitions")))
prop_trans <- do.call (cbind,prop_trans)
mean(apply (prop_trans,1,mean))
sd(apply (prop_trans,1,mean))
# transition rates
Trates <- (lapply (evol_rates, function (i) sapply(i,"[[", "diet_transitions_prop")))
Trates <- do.call (cbind,Trates)
mean(apply (Trates,1,mean))
sd(apply (Trates,1,mean))
# stasis time
Stime <- (lapply (evol_rates, function (i) sapply(i,"[[", "stasis_time")))
Stime <- do.call (cbind,Stime)
mean(apply (Stime,1,mean))
sd(apply (Stime,1,mean))
# last transition time
Ltime <- (lapply (evol_rates, function (i) sapply(i,"[[", "time_last_transitions")))
Ltime <- do.call (cbind,Ltime)
mean(apply (Ltime,1,mean))
sd(apply (Ltime,1,mean))

## obtaining only # Transition rates
transition_rates <- lapply(evol_rates, function (i) 
  lapply (i, function (k)
    
    k[,"diet_transitions_prop"])
  
  )
# melt to have a dataframe
transition_rates_df <- do.call (cbind.data.frame, transition_rates)

## obtaining only # Stasis time
stasis_time <- lapply(evol_rates, function (i) 
  lapply (i, function (k)
    
    k[,"stasis_time"])
  )
# melt to have a dataframe
stasis_time_df <- do.call (cbind.data.frame, stasis_time)

## obtaining only # Last transition  time
time_last_transitions <- lapply(evol_rates, function (i) 
  
  lapply (i, function (k)
    
    k[,"time_last_transitions"])
  
  )
# melt to have a dataframe
time_last_transitions_df <- do.call (cbind.data.frame, time_last_transitions)

# ----------------------------------------------------
# matching species names in phylogeny, trait and occurrence data
# transition rates data
rownames(transition_rates_df) [which (rownames(transition_rates_df) %in% rodent_data@data$binomial)]
row.names(transition_rates_df) [which (row.names(transition_rates_df) == "Akodon_serrensis")] <- "Castoria_angustidens"
row.names(transition_rates_df) [which (row.names(transition_rates_df) == "Paralomys_gerbillus")] <- "Phyllotis_gerbillus"
row.names(transition_rates_df) [which (row.names(transition_rates_df) == "Abrothrix_illuteus")] <- "Abrothrix_illutea"
row.names(transition_rates_df) [which (row.names(transition_rates_df) == "Akodon_latebricola")] <- "Neomicroxus_latebricola"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Abrothrix_andinus")] <- "Abrothrix_andina" 
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Abrothrix_lanosus")] <- "Abrothrix_lanosa" 
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Phyllotis_wolffsohni")] <- "Tapecomys_wolffsohni" 
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Akodon_bogotensis")] <- "Neomicroxus_bogotensis"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Sigmodontomys_aphrastus")] <- "Tanyuromys_aphrastus"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Chelemys_macronyx")] <- "Paynomys_macronyx"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Oligoryzomys_eliurus")] <- "Oligoryzomys_utiaritensis"
row.names(transition_rates_df) [which (row.names(transition_rates_df) ==  "Kunsia_fronto")] <- "Gyldenstolpia_fronto"

# stasis time data
row.names(stasis_time_df) [which (row.names(stasis_time_df) == "Akodon_serrensis")] <- "Castoria_angustidens"
row.names(stasis_time_df) [which (row.names(stasis_time_df) == "Paralomys_gerbillus")] <- "Phyllotis_gerbillus"
row.names(stasis_time_df) [which (row.names(stasis_time_df) == "Abrothrix_illuteus")] <- "Abrothrix_illutea"
row.names(stasis_time_df) [which (row.names(stasis_time_df) == "Akodon_latebricola")] <- "Neomicroxus_latebricola"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Abrothrix_andinus")] <- "Abrothrix_andina" 
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Abrothrix_lanosus")] <- "Abrothrix_lanosa" 
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Phyllotis_wolffsohni")] <- "Tapecomys_wolffsohni" 
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Akodon_bogotensis")] <- "Neomicroxus_bogotensis"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Sigmodontomys_aphrastus")] <- "Tanyuromys_aphrastus"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Chelemys_macronyx")] <- "Paynomys_macronyx"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Oligoryzomys_eliurus")] <- "Oligoryzomys_utiaritensis"
row.names(stasis_time_df) [which (row.names(stasis_time_df) ==  "Kunsia_fronto")] <- "Gyldenstolpia_fronto"

# last transition time data
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) == "Akodon_serrensis")] <- "Castoria_angustidens"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) == "Paralomys_gerbillus")] <- "Phyllotis_gerbillus"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) == "Abrothrix_illuteus")] <- "Abrothrix_illutea"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) == "Akodon_latebricola")] <- "Neomicroxus_latebricola"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Abrothrix_andinus")] <- "Abrothrix_andina" 
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Abrothrix_lanosus")] <- "Abrothrix_lanosa" 
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Phyllotis_wolffsohni")] <- "Tapecomys_wolffsohni" 
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Akodon_bogotensis")] <- "Neomicroxus_bogotensis"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Sigmodontomys_aphrastus")] <- "Tanyuromys_aphrastus"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Chelemys_macronyx")] <- "Paynomys_macronyx"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Oligoryzomys_eliurus")] <- "Oligoryzomys_utiaritensis"
row.names(time_last_transitions_df) [which (row.names(time_last_transitions_df) ==  "Kunsia_fronto")] <- "Gyldenstolpia_fronto"
# check
# rownames(time_last_transitions_df) [which (rownames(time_last_transitions_df) %in% rodent_data@data$binomial )]
# lapply (lapply (as.list (rownames(time_last_transitions_df) [which (rownames(time_last_transitions_df) %in% rodent_data@data$binomial != T)]),
#	             function (i) unlist(strsplit (i,"_"))[2]), function (i)
                       #	rodent_data@data$binomial [grep(i, rodent_data@data$binomial)]) # species not in one dataset
# rownames(time_last_transitions_df) [which (rownames(time_last_transitions_df) %in% rodent_data@data$binomial != T)]
rownames (transition_rates_df) <- rownames(time_last_transitions_df)

# finally matching occurrence and diet transtition datasets
stasis_time_df <- stasis_time_df [rownames(stasis_time_df) %in% rodent_data@data$binomial,]
transition_rates_df <- transition_rates_df [rownames(transition_rates_df) %in% rodent_data@data$binomial,]
time_last_transitions_df <- time_last_transitions_df [rownames(time_last_transitions_df) %in% rodent_data@data$binomial,]
# check
# rownames (transition_rates_df) == rownames(time_last_transitions_df)
# rownames (stasis_time_df) == rownames(time_last_transitions_df)

# ------------------------------------------------------------- #
# range maps for spcies with ancestral character estimates
# ------------------------------------------------------------- #

rodent_data  <- rodent_data  [which (rodent_data@data$binomial %in% rownames (time_last_transitions_df)),]

plot(rodent_data,
     col = rgb(0.5,0.1,0.1,alpha=0.05),
     border=rgb(0.1,0.1,0.1,alpha=0.2))# [grep ("Akodon", rodent_data@data$binomial),])

# ------------------------------------------------ #
# 		start  spatial analysis
#      cropping ecoregions according to rodent data
# ------------------------------------------------- #

## SA_ecor <- crop (ecor, rodent_data)

SA_ecor <- over (rodent_data,ecor,returnList=T)#
SA_ecor <- ecor [which (ecor@data$ECO_CODE %in% unique(unlist(lapply (SA_ecor, function (i) unique(i$ECO_CODE))))),]

#plot(SA_ecor)
#plot(rodent_data [grep ("Akodon", rodent_data@data$binomial),],col="green",add=T)
#plot(SA_ecor)
#plot(rodent_data [grep ("Oligoryzomys", rodent_data@data$binomial),],col="green",add=T)

# model raster from worldclim
#worldRaster <- getData('worldclim',var='alt',res=2.5)
worldRaster <- raster (here("wc2-5","alt.bil"))
# crop the world raster based on ecoregions
SA_raster <- crop (worldRaster,SA_ecor)
# disagregate to get the res of 0.25
SA_raster <-aggregate(SA_raster, fact=6.3) # 
# transform into equal area proj
SA_raster_EA<-projectRaster (SA_raster, 
               crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-56 +lat_0=-15 +datum=WGS84" )
# ecoregions too
SA_ecor_EA <- spTransform(SA_ecor,
                       CRS = "+proj=laea +x_0=0 +y_0=0 +lon_0=-56 +lat_0=-15 +datum=WGS84" )

## obter as coordenadas de cada celula
points_raster <- as.data.frame(coordinates (SA_raster_EA))
points_raster <- SpatialPointsDataFrame(coords = as.data.frame (points_raster), 
                                        data = points_raster,
                               proj4string = crs (SA_raster_EA))

# check if maps overlap
#plot(points_raster,cex=0.1,pch=19)
#plot(SA_ecor_EA,border="red",add=T)

## sobrepor pontos com as ecoregioes
over_pts_ecor <- lapply (as.list(SA_ecor_EA@data$ECO_CODE), function (i) 
  over (points_raster, SA_ecor_EA[which (SA_ecor_EA@data$ECO_CODE == i),]))
## remover o que n?o sobrepoe
over_pts_ecor <- lapply (over_pts_ecor, function (i) 
  points_raster [which (i$ECO_CODE != "NA"),])

#
## list of ecoregions
SA_ecor_list <- lapply (as.list (seq(1,length (SA_ecor_EA))), function (i) SA_ecor_EA[i,])
# unlist(lapply (SA_ecor_list, function (i) i$ECO_CODE)) == SA_ecor_EA@data$ECO_CODE

# measuring point distance to the ecoregion edge

################################
## parallel-processing settings
require (parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(raster))

# export your data and function
clusterExport(cl, c("SA_ecor_list",
                    "over_pts_ecor"))

distances_pts_to_edge <- parLapply (cl, seq (1,length (SA_ecor_list)), function (i)

  pointDistance(over_pts_ecor[[i]], geom(SA_ecor_list [[i]])[,c("x","y")])

)

stopCluster (cl)

# save
# saveRDS (distances_pts_to_edge, here ("Output","res_step2_distances_pts_to_edge.rda"))
## load if needed (all points and their distance to ecoregion edges)
distances_pts_to_edge <- readRDS (here ("Output",
                                        "res_step2_distances_pts_to_edge.rda"))

# removing ecoregions with too few points
# npoints per ecoregion
subsetting_data <- unlist(lapply(over_pts_ecor, function (i)
  nrow (i@coords)
))
over_pts_ecor2 <- over_pts_ecor[which(subsetting_data >= 20)] # subsetting points
distances_pts_to_edge2 <- distances_pts_to_edge[which(subsetting_data >= 20)] # subsetting distances
SA_ecor2 <- SA_ecor_EA[which(subsetting_data >= 20),] # subsetting ecoregions

## points with minimum distance from the edge
points_min_dist <- lapply (seq (1,length(over_pts_ecor2)), function (i) 
  # subsetting and putting in order
  as(over_pts_ecor2 [[i]][order(apply (distances_pts_to_edge2[[i]],1,max),decreasing=T),][1:10,],
   'SpatialPoints')) # to then transform into spatial points

# buffer around
buffer_min <- lapply (points_min_dist, function (i) 
	gBuffer(i, byid=T,width=13200))
names (buffer_min) <- SA_ecor2@data$ECO_CODE

## points of max distance from the edge
points_max_dist <- lapply (seq (1,length(over_pts_ecor2)), function (i) 
  # subsetting and putting in order
  as(over_pts_ecor2 [[i]][order(apply (distances_pts_to_edge2[[i]],1,min),decreasing=T),][1:10,],
     'SpatialPoints')) # to then transform into spatial points
# check
# apply (distances_pts_to_edge2[[i]],1,max)[order(apply (distances_pts_to_edge2[[i]],1,max),decreasing=T)][1:10]

# buffer around
buffer_max <- lapply (points_max_dist, function (i) 
  gBuffer(i, byid=T,width=13200))
names (buffer_max) <- SA_ecor2@data$ECO_CODE

# # check plot
# plot(SA_ecor2[1,])
# plot(points_min_dist[[i]],add=T,col="red",pch=3)
# plot(buffer_min[[i]],add=T,col="black",pch=1)
# plot(points_max_dist[[1]],add=T,pch=3)
# plot(buffer_max[[i]],add=T,col="black",pch=1)

## points of max distance from the edge
points_all_dist <- lapply (seq (1,length(over_pts_ecor2)), function (i) 
  # subsetting and putting in order
  as(over_pts_ecor2 [[i]][order(apply (distances_pts_to_edge2[[i]],1,min),decreasing=T),],
     'SpatialPoints')) # to then transform into spatial points

# buffer around all points
buffer_t <- lapply (points_all_dist, function (i) 
  gBuffer(i, byid=T,width=13200))
names (buffer_t) <- SA_ecor2@data$ECO_CODE

###### superimposing rodent data on points

# transform rodent data into lambert equal area projection
rodent_data_EA <- spTransform(rodent_data,
                              CRS = "+proj=laea +x_0=0 +y_0=0 +lon_0=-56 +lat_0=-15 +datum=WGS84" )


## exploracao da relacao entre a area da ecoregiao e o range dos roedores
# sobrepor os shapes

overlap_pol <- over (SA_ecor_EA, rodent_data_EA,returnList = TRUE)

## pegar a area da ecoregiao e a area do range das sp
area_ecor <- as.list(gArea(SA_ecor_EA,byid=T))#[which(unlist(lapply (overlap_pol,nrow)) != 0)])

## colocar em um DF o numero de sp com distribuicao menor do que a ecoregiao

area_rodent <- lapply (overlap_pol, function (i) 
    
  rodent_data_EA[which(rodent_data_EA@data$binomial %in% rownames(i)),]
  
  )
# area in lambert proj

area_rodent<- lapply (area_rodent, function (i) 
  
  gArea(i,byid=T)
  
)

# number of spp with range smaller than ecoregion area
sp_smaller_range <- lapply (seq(1,length(area_rodent)), function (i)
  cbind(smaller=table (area_rodent[[i]] <= area_ecor[[i]])[2],
        total=sum(table (area_rodent[[i]] <= area_ecor[[i]]))))

## proporcao de sp com o range menor do que a ecoregiao
prop_smaller <- unlist(lapply (seq(1,length (sp_smaller_range)), function (i)
  sp_smaller_range[[i]][,'smaller']/sp_smaller_range[[i]][,'total']))

## trocando NAs (i.e., nenhuma sp sobrepoe, ou nenhuma sp tem range menor do que a area da ecoregiao)
prop_smaller [is.na(prop_smaller)] <- 0
SA_ecor_EA@data <- cbind (SA_ecor_EA@data, 
                          prop_smaller=prop_smaller)

## 
require(ggplot2)

cores <- data.frame (cores=prop_smaller,NM_MUNICIP=SA_ecor_EA@data$ECO_ID_U)
f.mun<-fortify(SA_ecor_EA, region="ECO_ID_U")
f.mun<- cbind (f.mun, prop= cores [match (f.mun$id, cores$NM_MUNICIP),]$cores)

a.mean <- ggplot() + geom_polygon (data=SA_ecor_EA, aes(x=long, y=lat, group=group),
                                   size = 0.3, fill="gray85", colour="gray45",alpha=0.5) +
  coord_fixed (ratio = 1)+#xlim = c(-80, -34),  ylim = c(-58, 10), ) +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),axis.ticks.y=element_blank())
a.mean

b.mean <- a.mean+geom_polygon (data=f.mun, aes (x=long, y=lat,group=group,fill=prop),colour = "black", size=0.1) + 
  scale_fill_gradient2 (low='white',mid='red',high='darkred', midpoint = 0.5,
                        name="Proportion",
                        na.value="white",limits=c(0,1))+
  ggtitle ("Proportion of species with range size\n        smaller than ecoregion area")
b.mean

ggsave(filename = here ("Output",
                        "Figures_small_ranges",
                        "proportion_small_range.pdf"),family="serif")

###  the average geographic range size of assemblages in each bioregion


# -------------------------------------------------------------------------------
## descriptive statistics about range size, useful for sensitivity analysis
est_descr_ranges <-summary(gArea (rodent_data_EA,byid=T))

## especies de roedores com ranges menores do que o quartil de 25%
## (i.e., especies com range pequeno)
sp_small_range <- gArea (rodent_data_EA,byid=T)

# number of spp with range smaller than ecoregion area
sp_smaller_quartile <- lapply (seq(1,length(area_rodent)), function (i)
  cbind(smaller=table (area_rodent[[i]] <= est_descr_ranges[2])[2],
        total=sum(table (area_rodent[[i]] <= est_descr_ranges[2]))))

## proporcao de sp com o range menor do que a ecoregiao
prop_smaller_quartile<- unlist(lapply (seq(1,length (sp_smaller_quartile)), function (i)
  sp_smaller_quartile[[i]][,'smaller']/sp_smaller_quartile[[i]][,'total']))

## trocando NAs (i.e., nenhuma sp sobrepoe, ou nenhuma sp tem range menor do que a area da ecoregiao)
prop_smaller_quartile [is.na(prop_smaller_quartile)] <- 0
SA_ecor_EA@data <- cbind (SA_ecor_EA@data, 
                          prop_smaller_quartile=prop_smaller_quartile)
#

cores.b <- data.frame (cores=prop_smaller_quartile,
                       NM_MUNICIP=SA_ecor_EA@data$ECO_ID_U)
f.mun.b<-fortify(SA_ecor_EA, 
                 region="ECO_ID_U")
f.mun.b<- cbind (f.mun.b, 
                 prop_smaller_quartile= cores [match (f.mun$id, 
                                                      cores$NM_MUNICIP),]$cores)

require(ggplot2)
a.mean <- ggplot() + geom_polygon (data=SA_ecor_EA, aes(x=long, y=lat, group=group),
                                   size = 0.3, fill="gray85", colour="gray45",alpha=0.5) +
  coord_fixed (ratio = 1) +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),axis.ticks.y=element_blank())
a.mean

b.mean <- a.mean+geom_polygon (data=f.mun.b, 
                               aes (x=long, y=lat,group=group,
                                    fill=prop_smaller_quartile),
                               colour = "black", size=0.1) + 
  scale_fill_gradient2 (low='white',mid='red',high='darkred', midpoint = 0.12,
                        name="Proportion",
                        na.value="white",limits=c(0,max(prop_smaller_quartile)))+
  ggtitle ("Proportion of species with range size\n        smaller than the 1st quartile")
b.mean

ggsave(filename = here ("Output",
                        "Figures_small_ranges",
                        "proportion_small_range_quartile.pdf"),family="serif")


##########################################

# superimpose rodent and ecoregion data
sp_by_point_max <- lapply (buffer_max, function (i) over (i,rodent_data_EA, returnList = TRUE))
names(sp_by_point_max) <- SA_ecor2$ECO_CODE
sp_by_point_min <- lapply (buffer_min, function (i) over (i,rodent_data_EA, returnList = TRUE))
names(sp_by_point_min) <- SA_ecor2$ECO_CODE

### a identidade das ecoregioes e dos pontos 
# cores
nomes_ecor_max <- unlist(lapply (seq (1,length (names(sp_by_point_max ))), function (k) 
			         lapply (seq (1,length (sp_by_point_max[[k]])), function (i)
			            paste(names(sp_by_point_max)[k], i,sep="."))))
# ecotones
nomes_ecor_min <- unlist(lapply (seq (1,length (names(sp_by_point_min))), function (k) 
			lapply (seq (1,length (sp_by_point_min[[k]])), function (i)
			paste(names(sp_by_point_min)[k], i,sep="."))))

### selecionar sp exclusivas de cada porcao da eecoregiao
# core spp
core_species <- lapply (seq (1,length(sp_by_point_max)), function (i) # cada ecoregiao
			lapply (seq (1,10), function (k) # cada ponto
				tryCatch (			
				  sp_by_point_max[[i]][[k]] [which(sp_by_point_max[[i]][[k]]$binomial %in% sp_by_point_min[[i]][[k]]$binomial==F),]

					, error = function(e) return ("NULL")

      	  		)
			)
	)

# ecotone spp
ecotone_species <- lapply (seq (1,length(sp_by_point_min)), function (i) # cada ecoregiao
			  lapply (seq (1,10), function (k) # cada ponto
				tryCatch (			
				  sp_by_point_min[[i]][[k]] [which(sp_by_point_min[[i]][[k]]$binomial %in% sp_by_point_max[[i]][[k]]$binomial==F),]

					, error = function(e) return ("NULL")

      	  		)
			)
	)

## ajustar nomes
names (core_species)<- names(sp_by_point_max) 
names (ecotone_species)<- names(sp_by_point_min) 
# species id
cs<-unique(unlist (lapply (core_species, function (i) sapply(i, "[", "binomial"))))
es<-unique(unlist (lapply (ecotone_species, function (i) sapply(i, "[", "binomial"))))
# very few species are exclusive from either habitat
table(es %in% cs)
table(cs %in% es)

# removing points lacking species
# core
filtering_core_species <- lapply (core_species, function (i) {
	subset_dados <- which (unlist(lapply (i,nrow))>=1)
	filtragem_dados <- i [subset_dados]
	; filtragem_dados
	}
)

# ecotone
filtering_ecotone_species <- lapply (ecotone_species, function (i) {
	subset_dados <- which (unlist(lapply (i,nrow))>=1)
	filtragem_dados <- i [subset_dados]
	; filtragem_dados
	}
)

## id of ecoregions with data (ecor and points with more than 1spp)
## core
id_ecor_pts_core <- lapply (core_species, function (k) ## for each ecoregion
      which (
		lapply (k, function (i) ## each point
			which (nrow(i)>=1)) ## check if it has data
	==1))

## remove core points lacking data
id_ecor_pts_core  <- id_ecor_pts_core [lapply (id_ecor_pts_core,length) >0]

### ecotone
id_ecor_pts_ecotone <- lapply (ecotone_species, function (k) ## for each ecor
	which (
		lapply (k, function (i) ## for each point
			which (nrow(i)>=1)) ## check if it has data
	==1))

##  remove ecotone points lacking data
id_ecor_pts_ecotone <- id_ecor_pts_ecotone [lapply (id_ecor_pts_ecotone,length) >0]

# filtering out points lacking data
filtered_core_species <- filtering_core_species [lapply (filtering_core_species,length) >0] # core
filtered_ecotone_species <- filtering_ecotone_species [lapply (filtering_ecotone_species,length) >0] # ecotone

## list of species of each habtiat
list_core_spp <- unique(unlist(sapply (filtered_core_species , function (i) sapply (i, "[", "binomial"))))
list_ecotone_spp <- unique(unlist(sapply (filtered_ecotone_species, function (i) sapply (i, "[", "binomial"))))

# only 26 spp are exclusive of either habitat
table(
  list_core_spp %in% list_ecotone_spp
)

# statistics per habitat
# number of nodes
stat<- lapply (evol_rates, function (i) 
  lapply (i, function (k){
  
          dados <- k
          stat <- data.frame (meanCore=mean(dados [which(rownames(dados) %in% list_core_spp),"diet_total_nodes"]),
                              sdCore = sd(dados [which(rownames(dados) %in% list_core_spp),"diet_total_nodes"]),
                            meanEcotone=mean(dados [which(rownames(dados) %in% list_ecotone_spp),"diet_total_nodes"]),
                            sdEcotone =sd(dados [which(rownames(dados) %in% list_ecotone_spp),"diet_total_nodes"]))
          ;
          stat
}))

# average number of nodes per position in the ecoregion
sim_means <- do.call (rbind, lapply(stat, function (i) apply (do.call (rbind,i),2,mean)))
apply(sim_means,2,mean)

# diet_transitions_prop

stat<- lapply (evol_rates, function (i) 
  lapply (i, function (k){
    
    dados <- k
    stat <- data.frame (meanCore=mean(dados [which(rownames(dados) %in% list_core_spp),"diet_transitions_prop"]),
                        sdCore = sd(dados [which(rownames(dados) %in% list_core_spp),"diet_transitions_prop"]),
                        meanEcotone=mean(dados [which(rownames(dados) %in% list_ecotone_spp),"diet_transitions_prop"]),
                        sdEcotone =sd(dados [which(rownames(dados) %in% list_ecotone_spp),"diet_transitions_prop"]))
    ;
    stat
  }))

# average transition rates per position in the ecoregion
sim_means <- do.call (rbind, lapply(stat, function (i) apply (do.call (rbind,i),2,mean)))
apply(sim_means,2,mean)

# stasis_time

stat<- lapply (evol_rates, function (i) 
  lapply (i, function (k){
    
    dados <- k
    stat <- data.frame (meanCore=mean(dados [which(rownames(dados) %in% list_core_spp),"stasis_time"]),
                        sdCore = sd(dados [which(rownames(dados) %in% list_core_spp),"stasis_time"]),
                        meanEcotone=mean(dados [which(rownames(dados) %in% list_ecotone_spp),"stasis_time"]),
                        sdEcotone =sd(dados [which(rownames(dados) %in% list_ecotone_spp),"stasis_time"]))
    ;
    stat
  }))

# average stasis time
sim_means <- do.call (rbind, lapply(stat, function (i) apply (do.call (rbind,i),2,mean)))
apply(sim_means,2,mean)

# time_last_trans 

stat<- lapply (evol_rates, function (i) 
  lapply (i, function (k){
    
    dados <- k
    stat <- data.frame (meanCore=mean(dados [which(rownames(dados) %in% list_core_spp),"time_last_transitions"]),
                        sdCore = sd(dados [which(rownames(dados) %in% list_core_spp),"time_last_transitions"]),
                        meanEcotone=mean(dados [which(rownames(dados) %in% list_ecotone_spp),"time_last_transitions"]),
                        sdEcotone =sd(dados [which(rownames(dados) %in% list_ecotone_spp),"time_last_transitions"]))
    ;
    stat
  }))

# average last transition time
sim_means <- do.call (rbind, lapply(stat, function (i) apply (do.call (rbind,i),2,mean,na.rm=T)))
apply(sim_means,2,mean)

###################### 
# scale up species-level rates to point per ecoregion scale
# always considering different phylogenies and simulations

## transition rates
mean_tran_rates_min  <- lapply (seq(1,ncol(transition_rates_df)), function (m)
				  lapply (filtered_ecotone_species, function (k)
				    lapply (seq(1,length(k)), function (i)
					mean (transition_rates_df [which (rownames(transition_rates_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_tran_rates_minDF  <- lapply(mean_tran_rates_min,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_tran_rates_minECOT.rda"))

#transition rates max
mean_tran_rates_max  <- lapply (seq(1,ncol(transition_rates_df)), function (m)
				  lapply (filtered_core_species, function (k)
				    lapply (seq(1,length(k)), function (i)
					mean (transition_rates_df [which (rownames(transition_rates_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_tran_rates_maxDF  <- lapply(mean_tran_rates_max,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_tran_rates_maxCORE.rda"))

##########################
# stasis time min
mean_stasis_time_min  <- lapply (seq(1,ncol(stasis_time_df )), function (m)
				lapply (filtered_ecotone_species, function (k)
				lapply (seq(1,length(k)), function (i)
					mean (stasis_time_df [which (rownames(stasis_time_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_stasis_time_minDF <- lapply(mean_stasis_time_min,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_stasis_time_minECOT.rda"))
#
## stasis time max
mean_stasis_time_max  <- lapply (seq(1,ncol(stasis_time_df )), function (m)
				lapply (filtered_core_species, function (k)
				lapply (seq(1,length(k)), function (i)
					mean (stasis_time_df [which (rownames(stasis_time_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_stasis_time_maxDF <- lapply(mean_stasis_time_max,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_stasis_time_maxCORE.rda"))

#time last transition min
mean_last_trans_min  <- lapply (seq(1,ncol(time_last_transitions_df)), function (m)
				lapply (filtered_ecotone_species, function (k)
				lapply (seq(1,length(k)), function (i)
					mean (time_last_transitions_df [which (rownames(time_last_transitions_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_last_trans_minDF  <- lapply(mean_last_trans_min ,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_last_trans_minECOT.rda"))
#

mean_last_trans_max  <- lapply (seq(1,ncol(time_last_transitions_df)), function (m)
				  lapply (filtered_core_species, function (k)
				     lapply (seq(1,length(k)), function (i)
					 mean (time_last_transitions_df [which (rownames(time_last_transitions_df) %in% 
						k[[i]][,'binomial']),m],na.rm=T))))

mean_last_trans_maxDF  <- lapply(mean_last_trans_max,unlist)

saveRDS (mean_tran_rates_minDF, here ("Output","res_step2_mean_last_trans_maxCORE.rda"))

# --------------------------------------------------------------
# calculate phylogenetic diversity

tree_list <- tree <- read.nexus(here ("Data","PhylogenyTraits","Sigmodontinae_413species100Trees.trees"))

## corrigir os nomes
tree_list <- lapply (tree_list, function (i)
     {i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

## PD of core species
pd_core <-lapply (tree_list, function (tree)

  unlist(lapply (filtered_core_species, function (i)
    unlist(lapply (i, function (k)
  
    sum(tree$edge.length [which(tree$tip.label %in% k$binomial)])

  ))))
  )

## PD of ecotone spp
pd_ecotone <-lapply (tree_list, function (tree)
  
  unlist(lapply (filtered_ecotone_species, function (i)
    unlist(lapply (i, function (k)
      
      sum(tree$edge.length [which(tree$tip.label %in% k$binomial)])
      
    ))))
)

# put these data into a list
pd_list <- lapply (seq(1,length (pd_core)), function (i)
  
  rbind(data.frame (PD=pd_core[[i]],
                  position="core"),
      data.frame(PD=pd_ecotone[[i]],
                 position="ecotone"))
)

#lapply (pd_list,dim)

save (pd_list, file = here ("Output","res_step2_phylogenetic_div.RData"))

############################################################ 

# covariates of ecoregion points and neighborhood 

## come?ar a montar o dataframe com nomes e dados de evolu??o
#list_names <- lapply (as.list (SA_ecor2@data$ECO_CODE), function (i) rep (i,10))
#data_analysis <- data.frame (ecoreg=unlist(lapply (list_names, function (j) lapply (as.list(seq(1,10)), function (i) paste (j,i,sep=".")[1]))))

data_analysis <- rbind (data.frame (point=melt(id_ecor_pts_core)[,1],
			ecoreg = melt(id_ecor_pts_core)[,2],
			position="core"),
	data.frame (point=melt(id_ecor_pts_ecotone)[,1],
			ecoreg = melt(id_ecor_pts_ecotone)[,2],
			position="ecotone")
)


### filtering ecoregions with data for each position
# core
buffer_max_filtered <- buffer_max [which (names(buffer_max) %in% 
		data_analysis[which(data_analysis$position == "core"),]$ecoreg == T)]
buffer_min_filtered <- buffer_min [which (names(buffer_min) %in% 
		data_analysis[which(data_analysis$position == "ecotone"),]$ecoreg == T)]

### now removing points lacking data
buffer_max_filtered <- lapply (seq(1,length(id_ecor_pts_core)), function (i)
	
	buffer_max_filtered [[i]] [id_ecor_pts_core[[i]]]

)
names(buffer_max_filtered ) <- names(id_ecor_pts_core)

# Ecotones
buffer_min_filtered <- lapply (seq(1,length(id_ecor_pts_ecotone)), function (i)
	
	buffer_min_filtered [[i]] [id_ecor_pts_ecotone[[i]]]

)
names(buffer_min_filtered ) <- names(id_ecor_pts_ecotone)
#
##########################################
## area media das ecoregioes vizinhas dos pontos
# primeiro obter a area de todas as ecoregioes

areas <- lapply (SA_ecor2@polygons, function (i)
			
			i@area
		)

names (areas) <- SA_ecor2@data$ECO_CODE

# neighborhood
## if the neighbor is different or not of each ecoregion
## over to obtain ecoregion in each buffer around each point (considering core position)
nb_ecor_max <- lapply (buffer_max, function (i) 
	
	over(i,SA_ecor2, returnList=T)
)
# names (nb_ecor_max) == names(buffer_max )

# melt into a dataframe
nb_ecor_max_type <- lapply (nb_ecor_max, function (i) 

  	do.call (rbind.data.frame,i))

## Ordering ecoregion ids
## IDs in the correct order
nb_ecor_max_type <- do.call (rbind.data.frame, nb_ecor_max_type) [order (do.call (rbind.data.frame, nb_ecor_max_type)$count, decreasing=F),]
# ecoregion = neighbor ecoregion
nb_ecor_max_type <- cbind(nb_ecor_max_type, 
				number=substr(rownames(nb_ecor_max_type),8,12), 
				ecoregion=substr(rownames(nb_ecor_max_type),1,6))

# find the habitat of each ecoregion within the buffer
hab_max <- lapply (seq (1,nrow (nb_ecor_max_type)), function (i)
	paste(ecor@data [which (ecor@data$ECO_CODE %in% nb_ecor_max_type$ecoregion[i]) ,"Forest"]))
nb_ecor_max_type <- cbind(nb_ecor_max_type, 
                          hab_neigh=unlist (hab_max))

## saber qual ponto sobrepoe ecoregiao ABERTA ou FLORESTAL
ecor_over_max <- table (nb_ecor_max_type$number, nb_ecor_max_type$hab_neigh)
ecor_over_max <- matrix (ecor_over_max, ncol=2, byrow=F)
rownames (ecor_over_max) <- rownames(table (nb_ecor_max_type$number, nb_ecor_max_type$hab_neigh))
## botando em ordem dos pontos nas ecoregions
ecor_over_max <- ecor_over_max [match(unique (nb_ecor_max_type$number), rownames(ecor_over_max)),]
colnames(ecor_over_max) <- c("over_AB", "over_FOR")
#
ecor_over_max <- ecor_over_max[match (substr(rownames (nb_ecor_max_type), 8,15), rownames(ecor_over_max)),]
rownames(ecor_over_max) == substr(rownames (nb_ecor_max_type), 8,15)


## 
nb_ecor_max_type <- cbind(nb_ecor_max_type,ecor_over_max)

# dados das vizinhas
unique_count <- unique (nb_ecor_max_type$count)
compar <- lapply (seq(1,length (unique_count)), function (i)
	nb_ecor_max_type [which (nb_ecor_max_type$count == unique_count[i]),c("ECO_CODE", "ecoregion")])
vizinhas <- lapply (compar, function (i) 
	data.frame(area=ifelse (as.character (i [,1]) != as.character (i[,2]),  paste(i[,-1]),"NA"),
	row.names = rownames(i)))

### ID DAS VIZINHAS
id_vizinhas <-lapply (vizinhas, function (i) i[,"area"])

area_viz <- lapply (id_vizinhas, function (i) 
		lapply (seq(1,length (i)), function (k)
	
		 ifelse (which (names (areas) %in% as.character (i[k])),	
		        areas [which (names (areas) %in% as.character (i[k]))],
	      	  paste(i[k]))))
## 
area_viz <- lapply (area_viz, function (i) 
		
			lapply (seq(1,length(i)),  function (k) {
	
				ifelse (length(i[[k]])<1, NA, i[[k]])

		}
	))

nb_ecor_max_type <- cbind (nb_ecor_max_type, 
                           area_vizinho=unlist(area_viz))

## soma da ?rea dos vizinhos de cada ponto

sum_area_viz_max <- with(nb_ecor_max_type,
	aggregate (area_vizinho,by=list(number),FUN=sum,na.rm=T))
sum_over_ABE_viz_max <- with(nb_ecor_max_type,
	aggregate (over_AB,by=list(number),FUN=mean,na.rm=T))
sum_over_FOR_viz_max <- with(nb_ecor_max_type,
	aggregate (over_FOR,by=list(number),FUN=mean,na.rm=T))
sum_area_viz_max <- data.frame(Group.1=sum_area_viz_max [,1],
	area_viz=sum_area_viz_max [,2],
	over_AB=sum_over_ABE_viz_max[,2] ,
	over_FOR=sum_over_FOR_viz_max[,2] )
		
# sum_area_viz_max <- sum_area_viz_max [unique(nb_ecor_max_type$number),]

#------------------------------------------
#  neighbor of ecotone points
nb_ecor_min <- lapply (buffer_min, function (i) 
	over(i,SA_ecor2, returnList=T)
)

# transform list into DF
nb_ecor_min_type <- lapply (nb_ecor_min, function (i) 
	do.call (rbind.data.frame,i))

## ordenamento porque no meio do count de cada ecoregiao tem a id da regiao vizinha,
## assim, a ID da regiao vizinha vai ficar no lugar correto
nb_ecor_min_type <- do.call (rbind.data.frame, nb_ecor_min_type) [order (do.call (rbind.data.frame, nb_ecor_min_type)$count, decreasing=F),]
# ecoregion = ecoregiao vizinha
nb_ecor_min_type <- cbind(nb_ecor_min_type, 
				number=substr(rownames(nb_ecor_min_type),8,12), 
				ecoregion=substr(rownames(nb_ecor_min_type),1,6))

# a que habitat pertecem as ecoregioes dentro dos buffers

hab_min <- lapply (seq (1,nrow (nb_ecor_min_type)), function (i)
	paste(ecor@data [which (ecor@data$ECO_CODE %in% nb_ecor_min_type$ecoregion[i]) ,"Forest"]))

nb_ecor_min_type <- cbind(nb_ecor_min_type, hab_neigh=unlist (hab_min))

## saber o numero de ecoregioes vizinhas que sobrepoe, por tipo de hab ABERTA ou FLORESTAL
ecor_over_min <- table (nb_ecor_min_type$number, nb_ecor_min_type$hab_neigh)
ecor_over_min <- matrix (ecor_over_min, ncol=2, byrow=F)
rownames (ecor_over_min) <- rownames(table (nb_ecor_min_type$number, nb_ecor_min_type$hab_neigh))
## botando em ordem dos pontos nas ecoregions
ecor_over_min <- ecor_over_min [match(unique (nb_ecor_min_type$number), rownames(ecor_over_min)),]
colnames(ecor_over_min) <- c("over_AB", "over_FOR")

##
#
ecor_over_min <- ecor_over_min[match (substr(rownames (nb_ecor_min_type), 8,15), rownames(ecor_over_min)),]
table(rownames(ecor_over_min) == substr(rownames (nb_ecor_min_type), 8,15))

## 
nb_ecor_min_type <- cbind(nb_ecor_min_type,ecor_over_min)

## area media das vizinhas
# dados das vizinhas
unique_count <- unique (nb_ecor_min_type$count)
compar <- lapply (seq(1,length (unique_count)), function (i)
	nb_ecor_min_type [which (nb_ecor_min_type$count == unique_count[i]),c("ECO_CODE", "ecoregion")])

vizinhas <- lapply (compar, function (i) 
	data.frame(area=ifelse (as.character (i [,1]) != as.character (i[,2]),  paste(i[,-1]),"NA"),
	row.names = rownames(i)))

### ID DAS VIZINHAS
id_vizinhas <-lapply (vizinhas, function (i) i[,"area"])

area_viz <- lapply (id_vizinhas, function (i) lapply (seq(1,length (i)), function (k)
	 ifelse (which (names (areas) %in% as.character (i[k])),	
	        areas [which (names (areas) %in% as.character (i[k]))],
	        paste(i[k]))))
## 
area_viz <- lapply (area_viz, function (i) 
		
			lapply (seq(1,length(i)),  function (k) {
	
				ifelse (length(i[[k]])<1, NA, i[[k]])

		}
	)
)

nb_ecor_min_type <- cbind (nb_ecor_min_type, area_vizinho=unlist(area_viz))

## soma da ?rea dos vizinhos de cada ponto
sum_area_viz_min <- with(nb_ecor_min_type,
	aggregate (area_vizinho,by=list(number),FUN=sum,na.rm=T))
sum_over_ABE_viz_min <- with(nb_ecor_min_type,
	aggregate (over_AB,by=list(number),FUN=mean,na.rm=T))
sum_over_FOR_viz_min <- with(nb_ecor_min_type,
	aggregate (over_FOR,by=list(number),FUN=mean,na.rm=T))
sum_area_viz_min <- data.frame(Group.1=sum_area_viz_min [,1],
	area_viz=sum_area_viz_min [,2],
	over_AB=sum_over_ABE_viz_min[,2] ,
	over_FOR=sum_over_FOR_viz_min[,2] )
			
# sum_area_viz_min <- sum_area_viz_min [unique(nb_ecor_min_type$number),]

###  add habitat que cada ponto sobrepoe

#########################################
# nomes das ecoregioes, pontos, e vizinhos
sum_area_viz_max <- cbind (sum_area_viz_max, 
	ecoreg = SA_ecor2 [match (substr (sum_area_viz_max$Group.1,1,3), 
		rownames (SA_ecor2@data)),]$ECO_CODE,
	point = rep (seq(1,10),length(buffer_max)))

sum_area_viz_min <- cbind (sum_area_viz_min, 
	ecoreg =  SA_ecor2 [match (substr (sum_area_viz_min$Group.1,1,3), 
			rownames (SA_ecor2@data)),]$ECO_CODE,
	point = rep (seq(1,10),length(buffer_min)))

### colar os dados com caracteristicas das ecoregioes e seus vizinhos
# agora pegar as linhas (ponto, ecoreg) destas tabelas dos vizinhos para as quais temos dados

subset_neighborhood_data_max <- lapply (seq(1,length(id_ecor_pts_core)), function (i) {

	subset_ecoreg <- sum_area_viz_max [which(sum_area_viz_max$ecoreg %in% names(id_ecor_pts_core)[i]),]
	pontos_para_sel <- id_ecor_pts_core [[i]]	
	subset_points <- subset_ecoreg[which(subset_ecoreg$point %in% pontos_para_sel),]
	; subset_points 
	}
)


subset_neighborhood_data_min <- lapply (seq(1,length(id_ecor_pts_ecotone)), function (i) {

	subset_ecoreg <- sum_area_viz_min [which(sum_area_viz_min$ecoreg %in% names(id_ecor_pts_ecotone)[i]),]
	pontos_para_sel <- id_ecor_pts_ecotone [[i]]	
	subset_points <- subset_ecoreg[which(subset_ecoreg$point %in% pontos_para_sel),]
	; subset_points 
	}
)

## colar coordenadas
coord_max <- do.call(rbind.data.frame, lapply (buffer_max_filtered, coordinates))
coord_min <- do.call(rbind.data.frame, lapply (buffer_min_filtered, coordinates))

### colar os dados no DF de analise

data_analysis <- cbind(data_analysis,
		coordinates = rbind (coord_max,	
						            coord_min)
	)
	
## colar a area media das ecoregioes vizinhas
data_analysis<- cbind (data_analysis, 
				neigh_area=c(do.call(rbind,
							subset_neighborhood_data_max
							)$area_viz,
						do.call(rbind,
							subset_neighborhood_data_min
							)$area_viz)
			)

## colar o numero medio de vizinhas com veg aberta
data_analysis<- cbind (data_analysis, 
				over_AB=c(do.call(rbind,
							subset_neighborhood_data_max
							)$over_AB,
						do.call(rbind,
							subset_neighborhood_data_min
							)$over_AB)
			)

## colar o numero medio de vizinhas com veg florestal
data_analysis<- cbind (data_analysis, 
				over_FOR=c(do.call(rbind,
							subset_neighborhood_data_max
							)$over_FOR,
						do.call(rbind,
							subset_neighborhood_data_min
							)$over_FOR)
			)


## colar o tipo de habitat da ecoregiao foco
data_analysis<- cbind (data_analysis, 
	tipo_hab=SA_ecor2@data [match (substr (data_analysis$ecoreg, 1,6), 
			SA_ecor2@data$ECO_CODE),"Forest"])

## colar a interacao entre position e habitat type
data_analysis <- cbind(data_analysis, 
			int_pos_hab = 
				paste(data_analysis$position,
					data_analysis$tipo_hab,sep=""))

#### ver se a ecoregiao esta na mata atlantica
# from https://github.com/LEEClab/ATLANTIC-limits/blob/master/limites_wgs94.rar

mata_atlantica<-readOGR(dsn=here ("Data","Environment","MAtl"),layer="limite_ma_wwf_200ecprregions_wgs84")
# transform into lambert proj
mata_atlantica<-spTransform(mata_atlantica,
                            crs(SA_ecor2))
# crop ecor based on atlantic rainf
MA_ecor <- crop (SA_ecor2, mata_atlantica)
### excluindo aquelas ecoregioes de outros biomas que caem dentro da MA
## (base artigo renan ecography)
MA_ecor <- MA_ecor [1:10,] 
### inserindo no data_frame de variaveis
mata_atl <- ifelse (data_analysis$ecoreg %in% MA_ecor@data$ECO_CODE, 1,0)
data_analysis <- cbind (data_analysis, mata_atl= mata_atl)

### ver se a ecoregiao ocorre nos andes
andes<-readOGR(dsn=here ("Data","Environment",'andes'),
               layer="Lowenberg_Neto_2015")
# tranform into lambert projection
andes<- spTransform (andes,
                     crs(SA_ecor2))

## quais estao nos andes centrais (base artigo renan ecography)
subset_andes <- andes [-c(1,13,14,15,16,8,9,2,3,10,11),]
andes_ecor <- over (subset_andes,
                    SA_ecor2 )

## inserindo no data_frame de variaveis
andes_eco <- ifelse (data_analysis$ecoreg %in% andes_ecor$ECO_CODE, 1,0)
data_analysis <- cbind (data_analysis , andes_eco = andes_eco )

### update ecoregioes de analise
SA_ecor3 <- SA_ecor2[which (SA_ecor2$ECO_CODE %in% data_analysis$ecoreg),]

#______________________________________________
# area of ecoregions
mean(gArea (SA_ecor3,byid=T))/1000
sd(gArea (SA_ecor3,byid=T))/1000

# extent of the area
mean(gArea (SA_ecor[which(SA_ecor$ECO_CODE %in% SA_ecor3$ECO_CODE),],by=T))

# supporting information FigS1
# point position
png(here("Output",
         "Figures",
         "Fig_S1AB.png"), width=20, height=15, units="cm", res=300, family="serif")

par(mar=c(0,0,0,0),mfrow=c(1,2))
plot(SA_ecor3, lwd=1,col="gray95")
lapply (buffer_max, function (i)
  plot(i,pch=19, border = "green",add=T,col="green",cex=1))
lapply (buffer_min, function (i)
  plot(i,add=T,pch=19, border="red", col="red",cex=100))
# map.scale(x=-50,y=-42,cex=0.75,metric=T,relwidth = 0.15)

# bio regions
plot(SA_ecor3, lwd=1,col="gray95")
#map.scale(x=10,y=10,cex=0.65,metric=T,relwidth = 0.2)
plot(subset_andes,col="red",add=T)
plot(MA_ecor,col="green",add=T)

dev.off()


# ----------------------------------------- #
# Obtain species richness
## core
# list of species
list_sp_max <- lapply (filtered_core_species, function (i) 
	lapply (i, function (j) j$binomial))
# length is the richness
riqueza_max <- unlist(lapply (list_sp_max, function (i)## p cada ecoregiao
	lapply (i, function (j)  ## p cada ponto
		length (j)))) # calcule a riqueza
## ecotone
# list of spp
list_sp_min <- lapply (filtered_ecotone_species, function (i) 
	lapply (i, function (j) j$binomial))
# length is the richness
riqueza_min <- unlist(lapply (list_sp_min, function (i)## p cada ecoregiao
	lapply (i, function (j)  ## p cada ponto
		length (j)))) # calcule a riqueza
# total richness 
riqueza_total <- c(riqueza_max, riqueza_min)

### bind richness to the df of analysis
data_analysis <- cbind (data_analysis, 
	                      riqueza = riqueza_total)
# bind richness
# and save after updates
saveRDS (data_analysis, 
         file=here ("Output","res_step2_analysis_DF.rda"))


data_analysis <- readRDS(here ("Output","res_step2_analysis_DF.rda"))

# mean of covariates (Table 1)
# neighbors area
mean(data_analysis[,c("neigh_area")]/1000)
sd(data_analysis[,c("neigh_area")]/1000)

# number neighbors and habitat
round(apply(data_analysis[,c("over_AB","over_FOR")],2,mean),3)
apply(data_analysis[,c("over_AB","over_FOR")],2,sd)

# ------------------------------
# checking spatial autocorreation

# standardizing covarites
data_analysis$over_AB <-decostand(data_analysis$over_AB, "standardize")[,1]
data_analysis$over_FOR <- decostand(data_analysis$over_FOR, "standardize")[,1]
data_analysis$neigh_area <-sqrt(as.numeric(data_analysis$neigh_area))
data_analysis$neigh_area <-decostand(data_analysis$neigh_area, "standardize")[,1]
data_analysis$mata_atl <-as.factor(data_analysis$mata_atl)
data_analysis$andes_eco <-as.factor(data_analysis$andes_eco )
data_analysis$position <-as.factor(data_analysis$position)
data_analysis$tipo_hab <-as.factor(data_analysis$tipo_hab)
data_analysis$int_pos_hab <-as.factor(data_analysis$int_pos_hab)

# pooints to CAR model
points_CAR <-data_analysis[,c("coordinates.V1","coordinates.V2")]
coordinates (points_CAR) <- ~ coordinates.V1 + coordinates.V2
crs (points_CAR) <- crs(SA_ecor3)

points_CAR_pol <-gBuffer (points_CAR,
                          byid=T,
                          width=13200)

neigh <- poly2nb (points_CAR_pol)
winnb <- nb2WB(neigh)

as.matrix(model.matrix(~as.factor (data_analysis$position)-1))
as.matrix(model.matrix(~as.factor (data_analysis$tipo_hab)-1))

# 
with (data_analysis, 
      summary(lm(c(mean_tran_rates_maxDF[[1]],mean_tran_rates_minDF[[1]])
		~ int_pos_hab)))


# neighborhood

filt3<-mst.nb (dist(coordinates(points_CAR_pol)))
filt3

### testando autocorrelacao espacial
data_autoc_test <- data.frame(data_analysis,
			                    TR = c(rowMeans(do.call(cbind,mean_tran_rates_maxDF)),
			                           rowMeans(do.call(cbind,mean_tran_rates_minDF))),
			                    ST = c(rowMeans(do.call(cbind,mean_stasis_time_maxDF)),
			                           rowMeans(do.call(cbind,mean_stasis_time_minDF))),
			                    LT = c(rowMeans(do.call(cbind,mean_last_trans_maxDF)),
			                           rowMeans(do.call(cbind,mean_last_trans_minDF))))

#moran1<-lm.morantest(lme (log(TR) ~ position+tipo_hab, data=dados_teste),nb.filt3)

nb.filt3<-nb2listw (filt3)
nb.filt3 ## spatial weighting values
with(data_autoc_test, moran (TR, nb.filt3, 
                         n = nrow(coordinates(points_CAR_pol)), 
                         S0 = nrow(coordinates(points_CAR_pol))) )
# autocorrelation on TR
with(data_autoc_test, 
	moran.mc (TR, nb.filt3, nsim=999, zero.policy=NULL, alternative="greater", 
	          return_boot=FALSE))
# autocorrelation on ST
with(data_autoc_test, 
     moran.mc (ST, nb.filt3, nsim=999, zero.policy=NULL, alternative="greater", 
               return_boot=FALSE))
# autocorrelation on ST
with(data_autoc_test, 
     moran.mc (LT, nb.filt3, nsim=999, zero.policy=NULL, alternative="greater", 
               return_boot=FALSE))


# finally, bind data analysis dataset to each assemblage level tip based metric data set (produced by each phylogeny) 

# data aTR
dataTR <- lapply (seq (1,length(mean_tran_rates_maxDF)), function (i)
	
	cbind (data_analysis,
		TR = c(mean_tran_rates_maxDF[[i]],
			(mean_tran_rates_minDF[[i]])
	))
)

# data aST
dataST <- lapply (seq (1,length(mean_stasis_time_maxDF)), function (i)
	
	cbind (data_analysis,
		ST = c(mean_stasis_time_maxDF[[i]],
			(mean_stasis_time_minDF[[i]])
	))
)

# data aLT
dataLT <- lapply (seq (1,length(mean_tran_rates_maxDF)), function (i)
	
	cbind (data_analysis,
		LT = c(mean_last_trans_maxDF[[i]],
			(mean_last_trans_minDF[[i]])
	))
)

## complete dataset for models (trying to find the covariate(s) that most explain tip based metric values at assemblage scale)
save (dataTR,
      dataST,
      dataLT,
      file=here("Output", 
                "res_step2_data_for_GLMM.RData"))


