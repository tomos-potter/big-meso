# uploads the mesocosm data.
# for this project we consider one drainage: the Aripo
# this includes the Naranjo stream

upload_data <- function(){
library(data.table)
  
# upload main dataset
data <- fread("all_meso_data.csv")


#Rearrange and make some calculations
#Select fish from Aripo / Naranjo watershed only
data <- data[drainage=='Naranjo' | drainage=="Aripo",]
                   
# remove some experiments
# data <- data[experiment!="PhenoxDensx5050" & experiment!="PhenoxLight" # not using these as 14 day growth period
#              & experiment!= "DensxStructxPhenoxFreq" & experiment!="DensxPhenoxFreq",]

data <- data[experiment=="DensxPheno" | experiment=="PhenoxStruct"
             | experiment=="KilliPhenoxDens" | experiment=="KilliPhenoxSize",]



#Add in immature eggs
data$num.embryos2 <- data$num.embryos
data[which(is.na(data$num.embryos2)==TRUE),'num.embryos2'] <- 0

data$num.imm.eggs2 <- data$num.imm.eggs
data[which(is.na(data$num.imm.eggs2)==TRUE),'num.imm.eggs2'] <- 0

data$num.embryos <- data$num.embryos2 + data$num.imm.eggs2

data[which(data$final.sex=='M'),'num.embryos'] <- NA

data[which(data$survival==0),'num.embryos'] <- NA
data[which(data$survival==0),'embryo.dry.mass'] <- NA

data$num.embryos[which(is.na(data$num.embryos)==TRUE & data$final.sex=='F')] <- 0
# split into guppy and killifish dataframes

g <- data[spp=="guppy" & initial.fish==1,]
k <- data[spp=="killifish" & initial.fish==1,]


# now we need to split the killifish database further, because one experiment measured SL while the others were TL

k1 <- k[experiment=="GupRivCoevo",]
k2 <- k[experiment!="GupRivCoevo",]

k1$init.sl <- k1$initial.length
k1$final.sl <- k1$final.length


# for g2, we need to convert from TL to SL
# use regression coefficients from TL~SL model

# upload 2022 mesocosm data for killifish TL-SL relationship
tlsl_data <- fread("GK_data_meso_dry_2022_pt1.csv")
# subset just killifish data
killi_sl <- tlsl_data[species=="killifish",]



SL_mod <- lm(initial_sl ~ initial_tl*ecotype, data=killi_sl)

int <- SL_mod$coefficients[[1]]
coef_tl <- SL_mod$coefficients[[2]]
coef_eco <- SL_mod$coefficients[[3]]
coef_ecoTL <- SL_mod$coefficients[[4]]

k2$init.sl <- ifelse(k2$riv.pheno=="KO", 
                     int + coef_tl * k2$initial.total.length + coef_eco + coef_ecoTL*k2$initial.total.length,
                     int + coef_tl * k2$initial.total.length)

k2$final.sl <- ifelse(k2$riv.pheno=="KO", 
                      int + coef_tl * k2$final.total.length + coef_eco + coef_ecoTL*k2$final.total.length,
                      int + coef_tl * k2$final.total.length)

# combine them back together

k <- rbind(k1, k2)


#Calculate growth
#only survivors are measured

g$init.sl <- g$initial.length
g$final.sl <- g$final.length

data <- rbind(g, k)

data$growth <- log(data$final.sl / data$init.sl)


#Add in blocks
data[which(data$channel>=1 & data$channel<=4),'block'] <- 1
data[which(data$channel>=5 & data$channel<=8),'block'] <- 2
data[which(data$channel>=9 & data$channel<=12),'block'] <- 3
data[which(data$channel>=13 & data$channel<=16),'block'] <- 4
data$block <- as.factor(data$block)


#Create Channel.id which is unique to the experiment, drainage, replicate, and channel
data$channel.id <- paste(data$experiment, data$drainage,data$replicate,data$channel, sep = ".") 

#Order the data
data <- data[order(data$channel.id,data$init.sl),]



data <- data[which(data$initial.fish==1),c('experiment','drainage','block','channel','channel.id',
                                           'density','riv.density', 'structure', 'riv.size.structure',
                                           'initial.fish','analysis.use','died.early',
                                           'spp', 'riv.pheno','guppy.phenotype','final.sex','color.combo',
                                           'growth','survival','recap.total.length','pregnant','num.embryos',
                                           'stage','off.length','survival','init.sl')]

#Create KO, KG, HP, LP columns. These will get used as dummy variables
data$KO <- 0
data$KG <- 0
data[which(data$riv.pheno=='KO' & data$spp=="killifish"),'KO'] <- 1
data[which(data$riv.pheno=='LP' & data$spp=="killifish"),'KG'] <- 1

data$HP <- 0
data$LP <- 0
data[which(data$guppy.phenotype=='HP' & data$spp=="guppy"),'HP'] <- 1
data[which(data$guppy.phenotype=='LP' & data$spp=="guppy"),'LP'] <- 1

#data$drainage <- droplevels(data$drainage)

#unique(data$drainage)

#data[which(data$drainage==unique(data$drainage)[1]),'drain.num'] <- 1
#data[which(data$drainage==unique(data$drainage)[2]), 'drain.num'] <- 2

#Create gdata, which will be the data set that ultimately gets analyzed
touse <- which( is.na(data$died.early)==TRUE  )
gdata <- data[touse,c('experiment','channel.id','block','spp',
                      'final.sex','growth','num.embryos','pregnant',
                      'off.length','stage','survival','init.sl',
                      'KO','KG', 'LP', 'HP', 
                      'density','riv.density', 'structure','riv.size.structure')]

#Create a unique number for each channel
for (i in unique(gdata$channel.id)){
  gdata[which((gdata$channel.id)==i),'channel.num'] <- which(unique(gdata$channel.id)==i)
}

#Calculate the number of fish in each channel
gdata$count <- 1

densities <- gdata[, .(total.density = sum(count),
                       density.KO = sum(count * KO),
                       density.KG = sum(count*KG),
                       density.LP = sum(count*LP),
                       density.HP = sum(count*HP)), 
                   by=channel.num]

# densities$density.KO[which(density.KO==0)] <- 1
# densities$density.KG[which(density.KG==0)] <- 1
# densities$density.LP[which(density.LP==0)] <- 1
# densities$density.HP[which(density.HP==0)] <- 1

#densities <- ddply(gdata,c('channel.num','KO'),summarise,density=sum(count))
#densities <- reshape2::dcast(densities,channel.num~KO,value.var='density')
#names(densities) <- c('channel.num','density.LP','density.KO')
#densities[which(is.na(densities$density.LP)==T),'density.LP'] <- 0
#densities[which(is.na(densities$density.KO)==T),'density.KO'] <- 0

gdata <- merge(x=gdata,y=densities,by=c('channel.num'),all.x=TRUE)


#To be able to include experiments where there is only one phenotype, need to add the density of the same phenotype.
#This will not affect estimate are these are marked with 'both.pheno' in the model statement. It is just to keep it from throwing errors.
#by including NA's

# need to add in dummy variables for species x ecotype pairings
gdata$KO.HP <- ifelse(gdata$density.KO>0 & gdata$density.HP>0, 1, 0)

gdata$KO.LP <- ifelse(gdata$density.KO>0 & gdata$density.LP>0, 1, 0)

gdata$KG.HP <- ifelse(gdata$density.KG>0 & gdata$density.HP>0, 1, 0)

gdata$KG.LP <- ifelse(gdata$density.KG>0 & gdata$density.LP>0, 1, 0)


# gdata <- gdata[,c('exper','drain.num','block','final.sex','growth',
#                   'num.embryos','pregnant','off.length','stage',
#                   'survival','density.LP','density.KO','initial.total.length',
#                   'channel.num','KO','LP','both.pheno', 'riv.density', 'riv.size.structure')]
# 

#Not used for killifish
#gdata[which(is.na(gdata$off.length)==TRUE),'stage'] <- 0

gdata$growth[which(gdata$spp=="guppy" & gdata$final.sex=="M")]<- NA

gdata$num.embryos[which(is.na(gdata$num.embryos))] <- 0

gdata$num.embryos[which(gdata$final.sex=="M")]<- NA

gdata$growth[which(gdata$init.sl==8)] <- NA

return(gdata)
}

