###############################################################################################################
################################Polyphosphates - Adela Tupa's master thesis data###############################
###############################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#==================================Libraries
library(dplyr)
library(ggplot2)
library(reticulate)
library(readODS)
library(reshape)
library(gridExtra)
library(car)
library(emmeans)
#==================================GGPLOT THEME
#theme_min <- readRDS("/mnt/580CBE2464C5F83D/pracovni/helpfull_R_Matlab_script/ggtheme.rds")
theme_min <- readRDS("ggtheme.rds")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#=============================================================================================================#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Reading all data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#PP <- read_ods("/mnt/580CBE2464C5F83D/pracovni/data_statistika/AdelMGR/MGRAdelData.ods", 2)
PP <- read_ods("MGRAdelData.ods", 2)
PP$muGompertz <- RustParametry$mumax
PP$muL <- coef(RustParametryL)[, 3]
PP$Horizont <- factor(PP$Horizont, levels = c("O", "A"))
PP$Jezero <- factor(PP$Jezero, levels = c("PL", "CT"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#=====================Background data - graphs and statistics
#1. Soil pH
ggplot(PP, aes(factor(Plocha), pH)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  facet_grid(.~Horizont) + theme_min 

anova(glm(pH ~ Jezero, data = PP, family = "poisson"), test = "Chisq")
anova(glm(pH ~ Horizont, data = PP, family = "poisson"), test = "Chisq")
#2. Microbial biomass C, N and P
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(MBC1, na.rm = T),
                                                        ySD = sd(MBC1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min + xlab("Plocha") + ylab(expression(paste("MBC (",mu,"mol ",g^{-1}, ")")))
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(MBN1, na.rm = T),
                                                        ySD = sd(MBN1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(MBP1, na.rm = T),
                                                        ySD = sd(MBP1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 

anova(glm(MBC1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(MBN1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(MBP1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(MBC1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(MBN1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(MBP1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")

ggplot(PP, aes(MBC1, MBN1)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min
ggplot(PP, aes(MBN1, MBP1)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min

#3. Extractable C, N and P
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(TOC1, na.rm = T),
                                                        ySD = sd(TOC1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(TN1, na.rm = T),
                                                        ySD = sd(TN1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(Pi1, na.rm = T),
                                                        ySD = sd(Pi1, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 

anova(glm(TOC1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(TN1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(Pi1 ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(TOC1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(TN1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(Pi1 ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")

ggplot(PP, aes(TOC1, TN1)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min
ggplot(PP, aes(TN1, Pi1)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min


#Kinetic parameters of phosphatase activity
PPKinetics <- melt(PP[, c("Jezero", "Plocha", "Horizont", "VmaxNF", "VmaxF")], id.vars = c("Jezero", "Plocha", "Horizont"))
colnames(PPKinetics) <- c("Jezero", "Plocha", "Horizont", "Treatment", "Vmax")
PPKinetics$Treatment <- c(rep("Fresh", 36), rep("Fumigated", 36))
PPKinetics$Km <- c(as.numeric(PP[, "KmNF"]), as.numeric(PP[, "KmF"]))
PPKinetics$Ki <- c(as.numeric(PP[, "KiNF"]), as.numeric(PP[, "KiF"]))

PPKinetics %>% group_by(Jezero, Plocha, Horizont, Treatment) %>% summarise(y = mean(Vmax, na.rm = T),
                                                                           ySD = sd(Vmax, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, aes(color = Jezero, shape = Treatment)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 
PPKinetics %>% group_by(Jezero, Plocha, Horizont, Treatment) %>% summarise(y = mean(Km, na.rm = T),
                                                                           ySD = sd(Km, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, aes(color = Jezero, shape = Treatment)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 

PPKinetics %>% group_by(Jezero, Plocha, Horizont, Treatment) %>% summarise(y = mean(Ki, na.rm = T),
                                                                           ySD = sd(Ki, na.rm = T)) %>% 
  ggplot(aes(factor(Plocha), y)) + geom_point(cex = 6, aes(color = Jezero, shape = Treatment)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Horizont) + theme_min 

anova(glm(VmaxNF ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KmNF ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KiNF ~ Jezero, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(VmaxNF ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KmNF ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KiNF ~ Horizont, data = PP, family = "Gamma"), test = "Chisq")
anova(glm(VmaxNF ~ factor(Plocha), data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KmNF ~ factor(Plocha), data = PP, family = "Gamma"), test = "Chisq")
anova(glm(KiNF ~ factor(Plocha), data = PP, family = "Gamma"), test = "Chisq")

#=====================Analyzing the presence of PP - growth rates at zero P input
PPGrowth <- PPKinetics[, 1:3]
PPGrowth$Legend <- c(rep("Initial conditions", 36), rep("After incubation", 36))
PPGrowth$CFlush <- c(as.numeric(PP[, "MBC2"])*0.45, as.numeric(PP[, "muMBC"])*0.45)
PPGrowth$NFlush <- c(as.numeric(PP[, "MBN2"])*0.54, as.numeric(PP[, "muMBN"])*0.54)
PPGrowth$PFlush <- c(as.numeric(PP[, "MBP2"])*0.4, as.numeric(PP[, "muMBP"])*0.4)
PPGrowth$CPFlush <- PPGrowth$CFlush/PPGrowth$PFlush
PPGrowth$NPFlush <- PPGrowth$NFlush/PPGrowth$PFlush
PPGrowth$Horizont <- factor(PPGrowth$Horizont, levels = c("O", "A"))
PPGrowth$Legend <- factor(PPGrowth$Legend, levels = c("Initial conditions", "After incubation"))
PPGrowth$Jezero <- factor(PPGrowth$Jezero, levels = c("PL", "CT"))

PPGrowth %>% group_by(Jezero, Horizont, Legend) %>% summarise(y = mean(CFlush/PFlush, na.rm = T),
                                                                      ySD = sd(CFlush/PFlush, na.rm = T)) %>% 
  ggplot(aes(Horizont, y)) + geom_point(cex = 6, aes(color = Jezero, shape = Legend)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_grid(.~Jezero) + theme_min 
facetLabs <- c('A'="Organic Horizon",
               'O'="Litter Layer")
#C to P ratio
ggplot(PPGrowth, aes(Jezero, CFlush/PFlush)) + geom_boxplot(position = "dodge", aes(fill = Legend)) +
  facet_grid(.~Horizont, labeller = as_labeller(facetLabs)) + theme_min + 
  scale_fill_manual(values = c("white", "grey")) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.9), axis.title.x = element_blank()) +
  scale_x_discrete(labels=c('Plešné', 'Čertovo')) + ylab(expression(paste(frac(C, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
  ggtitle("A)")
#N to P ratio
ggplot(PPGrowth, aes(Jezero, NFlush/PFlush)) + geom_boxplot(position = "dodge", aes(fill = Legend)) +
  facet_grid(.~Horizont, labeller = as_labeller(facetLabs)) + theme_min + 
  scale_fill_manual(values = c("white", "grey")) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.9), axis.title.x = element_blank()) +
  scale_x_discrete(labels=c('Plešné', 'Čertovo')) + ylab(expression(paste(frac(N, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
  ggtitle("B)")
#Composite figure - 4 to 12 in landscape
grid.arrange(
  ggplot(PPGrowth, aes(Jezero, CFlush/PFlush)) + geom_boxplot(position = "dodge", aes(fill = Legend)) +
    facet_grid(.~Horizont, labeller = as_labeller(facetLabs)) + theme_min + 
    scale_fill_manual(values = c("white", "grey")) +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.85), axis.title.x = element_blank()) +
    scale_x_discrete(labels=c('Plešné', 'Čertovo')) + ylab(expression(paste(frac(C[Flush], P[Flush])~(mol/mol)))) +
    ggtitle("A)"),
  ggplot(PPGrowth, aes(Jezero, NFlush/PFlush)) + geom_boxplot(position = "dodge", aes(fill = Legend), show.legend = F) +
    facet_grid(.~Horizont, labeller = as_labeller(facetLabs)) + theme_min + 
    scale_fill_manual(values = c("white", "grey")) +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.9), axis.title.x = element_blank()) +
    scale_x_discrete(labels=c('Plešné', 'Čertovo')) + ylab(expression(paste(frac(N[Flush], P[Flush])~(mol/mol)))) +
    ggtitle("B)"), ncol = 2
)
#Statistical test
##Before incubation
anova(glm(CPFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
anova(glm(CPFlush ~ Jezero, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
anova(glm(CPFlush ~ Plocha, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
anova(glm(NPFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
anova(glm(NPFlush ~ Jezero, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
anova(glm(NPFlush ~ Plocha, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
summary(PPGrowth[PPGrowth$Legend == "Initial conditions", "CPFlush"])
summary(PPGrowth[PPGrowth$Legend == "Initial conditions", "NPFlush"])
##===============Calculating average Poly-P contribution to MBP
Pfmean <- coef(glm(PFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
Cfmean <- coef(glm(CFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "Initial conditions"), test = "F")
#Litter horizon
(0.4 - 1/Pfmean[1]/(1/Cfmean[1]/0.24/84))/(0.4 - 0.84 - 1/Pfmean[1]/(1/Cfmean[1]/0.24/84))
#Organic horizon
(0.4 - 1/(Pfmean[1] + Pfmean[2])/(1/(Cfmean[1] + Cfmean[2])/0.24/84))/(0.4 - 0.93 - 1/(Pfmean[1] + Pfmean[2])/(1/(Cfmean[1] + Cfmean[2])/0.24/84))
##=============================================================
##After incubation
anova(glm(CPFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
anova(glm(CPFlush ~ Jezero, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
anova(glm(CPFlush ~ Plocha, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
anova(glm(NPFlush ~ Horizont, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
anova(glm(NPFlush ~ Jezero, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
anova(glm(NPFlush ~ Plocha, data = PPGrowth, family = "Gamma", subset = Legend == "After incubation"), test = "F")
summary(PPGrowth[PPGrowth$Legend == "After incubation", "CPFlush"])
summary(PPGrowth[PPGrowth$Legend == "After incubation", "NPFlush"])
##Before and after
anova(glm(CPFlush ~ Legend, data = PPGrowth, family = "Gamma"), test = "F")
anova(glm(CPFlush ~ Horizont, data = PPGrowth, family = "Gamma"), test = "F")
anova(glm(CPFlush ~ Legend*Horizont, data = PPGrowth, family = "Gamma"), test = "F")
anova(glm(NPFlush ~ Legend, data = PPGrowth, family = "Gamma"), test = "F")
anova(glm(NPFlush ~ Horizont, data = PPGrowth, family = "Gamma"), test = "F")
anova(glm(NPFlush ~ Legend*Horizont, data = PPGrowth, family = "Gamma"), test = "F")

#Growth rate at zero Pi should depend on initial C to P or N to P ratio
PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(muL*24, na.rm = T),
                                                        ySD = sd(muL*24/3, na.rm = T),
                                                        x = mean(log(dMBC), na.rm = T),
                                                        xSD = sd(MBC2*0.45/MBP2/0.4/3, na.rm = T)) %>% 
  ggplot(aes(x, y)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min +
  stat_smooth(method = "lm", se = F, col = "black") + 
  scale_fill_manual(values = c("black", "grey")) + 
  xlab(expression(paste(Initial~frac(C[Flush], P[Flush])~(mol/mol)))) +
  ylab(expression(paste(mu[P0]~(d^{-1})))) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) + ggtitle("A)")

PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(muL*24, na.rm = T),
                                                        ySD = sd(muL*24/3, na.rm = T),
                                                        x = mean(MBN2*0.54/MBP2/0.4, na.rm = T),
                                                        xSD = sd(MBN2*0.54/MBP2/0.4*24/3, na.rm = T)) %>% 
  ggplot(aes(x, y)) + geom_point(cex = 6, pch = 21, aes(fill = Jezero), show.legend = F) + theme_min +
  stat_smooth(method = "lm", se = F, col = "black") + 
  scale_fill_manual(values = c("black", "grey"), labels = c("Litter layer", "Organic horizon")) + 
  xlab(expression(paste(Initial~frac(N[Flush], P[Flush])~(mol/mol)))) +
  ylab(expression(paste(mu[P0]~(d^{-1})))) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) + ggtitle("B)")

PP$CPFlush <- PP$MBC2*0.45/PP$MBP2/0.4
PP$NPFlush <- PP$MBN2*0.54/PP$MBP2/0.4

anova(lm(muL ~ Horizont, PP))


#Composite figure - 8 to 12 inches landscape
grid.arrange(PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(muL*24, na.rm = T),
                                                                     ySD = sd(muL*24, na.rm = T),
                                                                     x = mean(MBC2*0.45/MBP2/0.4, na.rm = T),
                                                                     xSD = sd(MBC2*0.45/MBP2/0.4/3, na.rm = T)) %>% 
               ggplot(aes(x, y)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min +
               stat_smooth(method = "lm", se = F, aes(col = Horizont), show.legend = F) + 
               scale_color_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) + 
               scale_fill_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) + 
               xlab(expression(paste(Initial~frac(C[Flush], P[Flush])~(mol/mol)))) +
               ylab(expression(paste(mu[P0]~(day^{-1})))) +
               geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
               annotate("text", x = 23.5, y = 0.7, label = "***Horizon", cex = 6) +
               #annotate("text", x = 24, y = 0.74, parse = T, label = as.character(expression(paste("*",Initial~frac(C[Flush], P[Flush])))), cex= 6) +
               theme(legend.title = element_blank(), legend.position = c(0.18, 0.1),
                     legend.background=element_rect(fill = alpha("white", 0.01))) + ggtitle("A)"),
             
             PP %>% group_by(Jezero, Plocha, Horizont) %>% summarise(y = mean(muL*24, na.rm = T),
                                                                     ySD = sd(muL*24, na.rm = T),
                                                                     x = mean(MBN2*0.54/MBP2/0.4, na.rm = T),
                                                                     xSD = sd(MBN2*0.54/MBP2/0.4/3, na.rm = T)) %>% 
               ggplot(aes(x, y)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont), show.legend = F) + theme_min +
               stat_smooth(method = "lm", se = F, aes(col = Horizont), show.legend = F) + 
               scale_color_manual(values = c("black", "grey"), labels = c("Litter layer", "Organic horizon")) +
               scale_fill_manual(values = c("black", "grey"), labels = c("Litter layer", "Organic horizon")) + 
               xlab(expression(paste(Initial~frac(N[Flush], P[Flush])~(mol/mol)))) +
               annotate("text", x = 3.35, y = 0.7, label = "***Horizon", cex = 6) +
               #annotate("text", x = 3.4, y = 0.74, parse = T, label = as.character(expression(paste("*",Initial~frac(N[Flush], P[Flush])))), cex= 6) +
               ylab(expression(paste(mu[P0]~(day^{-1})))) +
               geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
               theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) + ggtitle("B)"),
             ncol = 2)

ggplot(PP, aes(MBN2*0.54/MBP2/0.4, muGompertz*24)) + geom_point(aes(color = Jezero))

##Relative change of MBC
PP$wSoil <- c(rep(0.88, 3), rep(1.07, 3), rep(1.31, 3), rep(2.01, 3), rep(0.65, 3), 
              rep(1.13, 3), rep(0.6, 3), rep(0.9, 3), rep(1, 3), rep(1.62, 3), rep(1.3, 3), 
              rep(1.26, 3)) 

PP$dMBC <- with(PP, (muMBC*0.45*(25 + wSoil*DW) - MBC2*0.45*wSoil*DW))
PP$dMBN <- with(PP, (muMBN*0.54*(25 + wSoil*DW) - MBN2*0.54*wSoil*DW))
PP$dMBP <- with(PP, (muMBP*0.4*(25 + wSoil*DW) - MBP2*0.4*wSoil*DW))

#Is the dMBC, dMBN and dMBP statistically significant from zero?
t.test(PP$dMBC)
t.test(PP$dMBN)
t.test(PP$dMBP)
#====================================Confronting data with microscope pictures
##Creating empty data matrix for pictures 
mic_data_matrix <- NULL

##Creating list of paths to representative pictures
pths <- list.files("/mnt/580CBE2464C5F83D/pracovni/data_statistika/AdelMGR/Mikroskop/", full.names = T)

##Extracting pictures data
for(i in pths){
  cimage <- image_read(i)
  cdata <- as.numeric(cimage[[1]][1, , ])
  mic_data_matrix <- rbind(mic_data_matrix, cdata)
}

##Data frame with data (matrix raws) identifications
mic_data_ident <- data.frame(Jezero = c(rep("CT", 12), rep("PL", 12)),
                             Plocha = c(rep(25, 4), rep(265, 4), rep(74, 4), rep(114, 4), rep(147, 4), rep(150, 4)),
                             Horizont = rep(c("A", "A", "O", "O"), 6),
                             Legend = rep(c("Initial conditions", "After incubation"), 12)
)

##Calculating distance matrix (Bray distance is used)
mic_dist <- vegdist(mic_data_matrix)

##Calculating unconstrained principal components
mic_pca <- prcomp(mic_dist, scale. = T)
##Extracting principal components and binding it with mic_data_ident
mic_pca_out <- cbind(mic_data_ident, mic_pca$x)
##Adding medians of measured CPFlush and NPFlush 
stoich <- as.data.frame(PPGrowth %>% group_by(Jezero, Plocha, Horizont, Legend) %>% summarise(CPflush = median(CPFlush, na.rm = T),
                                                                                              NPflush = median(NPFlush, na.rm = T)))
mic_pca_out <- merge(stoich, mic_pca_out, by.x = c("Jezero", "Plocha", "Horizont", "Legend"),
                     by.y = c("Jezero", "Plocha", "Horizont", "Legend"))
##Plotting the results
###C to P flush
ggplot(mic_pca_out[-c(10, 12), ], aes(PC1, CPflush)) + geom_point(cex = 6, pch = 21, aes(fill = Legend)) + theme_min +
  facet_wrap(~Horizont, scales = "free_y", labeller = as_labeller(facetLabs)) + 
  scale_fill_manual(values = c("white", "grey")) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.1)) +
  stat_smooth(method = lm, se = F, formula = y ~ x, color = "black") + 
  ylab(expression(paste(frac(C, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
  xlab("PCA 1 score") + 
  ggtitle("A)") + geom_point(data = mic_pca_out[c(10, 12), ],  aes(PC1, CPflush), cex = 4, pch = 21, fill = "black")
###N to P flush
ggplot(mic_pca_out[-c(10, 12), ], aes(PC1, NPflush)) + geom_point(cex = 6, pch = 21, aes(fill = Legend), show.legend = F) + theme_min +
  facet_wrap(~Horizont, scales = "free_y", labeller = as_labeller(facetLabs)) + 
  scale_fill_manual(values = c("white", "grey")) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.1)) +
  stat_smooth(method = lm, se = F, formula = y ~ x, color = "black") + 
  ylab(expression(paste(frac(N, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
  xlab("PCA 1 score") + 
  ggtitle("B)") + geom_point(data = mic_pca_out[c(10, 12), ],  aes(PC1, NPflush), cex = 4, pch = 21, fill = "black")

##Statistics
summary(lm(CPflush ~ PC1, mic_pca_out[-c(10, 12), ]))
summary(lm(NPflush ~ PC1, mic_pca_out[-c(10, 12), ]))
anova(lm(CPflush ~ PC1, mic_pca_out[-c(10, 12), ]))
anova(lm(NPflush ~ PC1, mic_pca_out[-c(10, 12), ]))

#Composite figure - 4 to 12 inches landscape
grid.arrange(
  ggplot(mic_pca_out[-c(10, 12), ], aes(PC1, CPflush)) + geom_point(cex = 6, pch = 21, aes(fill = Legend)) + theme_min +
    facet_wrap(~Horizont, scales = "free_y", labeller = as_labeller(facetLabs)) + 
    scale_fill_manual(values = c("white", "grey")) +
    theme(legend.title = element_blank(), legend.position = c(0.85, 0.15)) +
    stat_smooth(method = lm, se = F, formula = y ~ x, color = "black") + 
    ylab(expression(paste(frac(C, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
    xlab("PCA 1 score") + 
    ggtitle("A)") + geom_point(data = mic_pca_out[c(10, 12), ],  aes(PC1, CPflush), cex = 4, pch = 21, fill = "black"),
  ggplot(mic_pca_out[-c(10, 12), ], aes(PC1, NPflush)) + geom_point(cex = 6, pch = 21, aes(fill = Legend), show.legend = F) + theme_min +
    facet_wrap(~Horizont, scales = "free_y", labeller = as_labeller(facetLabs)) + 
    scale_fill_manual(values = c("white", "grey")) +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.1)) +
    #stat_smooth(method = lm, se = F, formula = y ~ x, color = "black") + 
    ylab(expression(paste(frac(N, P)~of~CHCl[3], " Flush ", (frac(mol,mol))))) +
    xlab("PCA 1 score") + 
    ggtitle("B)") + geom_point(data = mic_pca_out[c(10, 12), ],  aes(PC1, NPflush), cex = 4, pch = 21, fill = "black"),
  ncol = 1
)
#===============================================================================
#=====================Recovery of PP from extract and the soil
Recoveries <- rbind(PP[, c(1:25, 37:43)], PP[, c(1:25, 37:43)],
                    PP[, c(1:25, 37:43)], PP[, c(1:25, 37:43)])
Recoveries$Legend <- c(rep("Extract - NF", nrow(PP)), rep("Soil - NF", nrow(PP)),
                       rep("Extract - F", nrow(PP)), rep("Soil - F", nrow(PP)))
Recoveries$Recovery <- c(PP$PiExcessNFExtract/PP$PPExtractNF*100, PP$PiExcessNF/PP$Navratnost/PP$PPadd*100,
                         PP$PiExcessFExtract/PP$PPExtractF*100, PP$PiExcessF/PP$Navratnost/PP$PPadd*100)
Recoveries$Soil <- c(rep("No", nrow(PP)), rep("Yes", nrow(PP)),
                     rep("No", nrow(PP)), rep("Yes", nrow(PP)))
Recoveries$Fumigation <- c(rep("No", nrow(PP)), rep("No", nrow(PP)),
                           rep("Yes", nrow(PP)), rep("Yes", nrow(PP)))

glm0 <- glm(Recovery ~ 1, Recoveries, family = Gamma)
add1(glm0, .~.+Jezero+Plocha+Horizont+Soil+Fumigation+Legend+pH, test = "Chisq")

glm1 <- update(glm0, .~.+Legend)
summary(glm1)

em <- emmeans(glm1, "Legend")
contrast(em, "pairwise", adjust = "Tukey")

add1(glm1, .~.*Jezero+Plocha+Horizont+Legend+pH, test = "Chisq")

glm2 <- update(glm1, .~.+Horizont)
summary(glm2)

em2 <- emmeans(glm2, c("Legend", "Horizont"))
summary(em2)
paste0(round(1/summary(em2)$emmean[1], 1),"%")

Recoveries$Legend <- factor(Recoveries$Legend, levels = c("Extract - NF", "Extract - F", "Soil - NF", "Soil - F"))

Recoveries$RecoveryMean <- NA

for(i in 1:nrow(Recoveries)){
  if(Recoveries$Horizont[i] == "O" & Recoveries$Legend[i] == "Extract - F"){
    Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[1], 1),"%")
  }else{
    if(Recoveries$Horizont[i] == "O" & Recoveries$Legend[i] == "Extract - NF"){
      Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[2], 1),"%")
    }else{
      if(Recoveries$Horizont[i] == "O" & Recoveries$Legend[i] == "Soil - F"){
        Recoveries$RecoveryMean[i] <- paste0("***", round(1/summary(em2)$emmean[3], 1),"%")
      }else{
        if(Recoveries$Horizont[i] == "O" & Recoveries$Legend[i] == "Soil - NF"){
          Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[4], 1),"%")
        }else{
          if(Recoveries$Horizont[i] == "A" & Recoveries$Legend[i] == "Extract - F"){
            Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[5], 1),"%")
          }else{
            if(Recoveries$Horizont[i] == "A" & Recoveries$Legend[i] == "Extract - NF"){
              Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[6], 1),"%")
            }else{
              if(Recoveries$Horizont[i] == "A" & Recoveries$Legend[i] == "Soil - F"){
                Recoveries$RecoveryMean[i] <- paste0("***",round(1/summary(em2)$emmean[7], 1),"%")
              }else{
                Recoveries$RecoveryMean[i] <- paste0(round(1/summary(em2)$emmean[8], 1),"%")
              }
            }
          }
        }
      }
    }
  }
}

ggplot(Recoveries, aes(Legend, Recovery)) + geom_boxplot(position = "dodge", aes(fill = Fumigation)) +
  theme_min + facet_grid(.~Horizont, labeller = as_labeller(facetLabs)) +
  scale_fill_manual(values = c("white", "grey"), 
                    labels = c(expression(paste("-",CHCl[3])), expression(paste("+",CHCl[3])))) +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.2), axis.title.x = element_blank()) +
  geom_hline(yintercept = 100, lwd = 1.2) + ylim(0, 110) +
  ylab("Poly-P recovery (%)") + geom_text(aes(x = Legend, y = 110, label = RecoveryMean), cex = 5, check_overlap = T,
                                          fontface = "bold") +
  scale_x_discrete(labels = c("Extract", "Extract", "Soil", "Soil"))
  
#====================================Predict PP depolymerization by phosphatase activity
PP$PiExcessFtoPPadd <- with(PP, PiExcessF/PPadd)

ggplot(PP[PP$PiExcessFtoPPadd<1, ], aes(VmaxNF, PiExcessF/PPadd)) + 
  geom_point(cex = 6, pch = 21, aes(fill = Horizont)) + theme_min +
  scale_fill_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.85)) +
  ylab("Poly-P recovery (%)") +
  xlab(expression(paste(V[MAX], " (", mu,"mol MUB-P  ",g~DW^{-1}~d^{-1}, ")"))) +
  ggtitle("A)")

#===============================================================
PPdynamic<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Deplymerization
    depolymerization <- Vmax*PP/(Km*(1 + Pi/Ki) + PP)
    #Derivatives
    dPP<- - depolymerization
    dPi<- depolymerization 
    
    return(list(c(dPP, dPi)))
  })
}

out <- ode(y=c(PP=PP$PPadd[10], Pi=PP$Pi1[10]), parms=c(Vmax = PP$VmaxF[10], Km = PP$KmF[10], Ki = PP$KiF[10]), 
           PPdynamic, times=seq(0, 24)) #Pocital jsem tu aktivitu na hodinu nebo den?
plot(out)
out[2, 3] - out[1, 3]
PP$PiExcessF[10]

PP$PiExcessPred <- NA
for(i in 1:nrow(PP)){
  outODE <- ode(y=c(PP=PP$PPadd[i], Pi=PP$Pi1[i]), parms=c(Vmax = PP$VmaxNF[i], Km = PP$KmNF[i], Ki = PP$KiF[i]), 
             PPdynamic, times=seq(0, 1))
  PP$PiExcessPred[i] <- outODE[2, 3] - outODE[1, 3]
}

ggplot(PP, aes(PiExcessPred, PiExcessF)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont), show.legend = F) + 
  theme_min + geom_abline(intercept = 0, slope = 1) +
  scale_fill_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.15)) + 
  ylab(expression(paste("Observed ", P[i], " excess (", mu, "mol P g ", DW^{-1}, ")"))) +
  xlab(expression(paste("Predicted ", P[i], " excess (", mu, "mol P g ", DW^{-1}, ")"))) +
  ggtitle("B)")
#====================================================================
#Composite figure
grid.arrange(
  ggplot(PP[PP$PiExcessFtoPPadd<1, ], aes(VmaxNF, PiExcessF/PPadd)) + 
    geom_point(cex = 6, pch = 21, aes(fill = Horizont), show.legend = F) + theme_min +
    scale_fill_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) +
    theme(legend.title = element_blank(), legend.position = c(0.15, 0.85)) +
    ylab("Poly-P recovery (%)") +
    xlab(expression(paste(V[MAX], " (", mu,"mol MUB-P  ",g~DW^{-1}~d^{-1}, ")"))) +
    ggtitle("A)"),
  ggplot(PP, aes(PiExcessPred, PiExcessF)) + geom_point(cex = 6, pch = 21, aes(fill = Horizont), show.legend = T) + 
    theme_min + geom_abline(intercept = 0, slope = 1) +
    scale_fill_manual(values = c("black", "grey"), labels = c("Litter Layer", "Organic Horizon")) +
    theme(legend.title = element_blank(), legend.position = c(0.85, 0.15)) + 
    ylab(expression(paste("Observed ", P[i], " excess (", mu, "mol P g ", DW^{-1}, ")"))) +
    xlab(expression(paste("Predicted ", P[i], " excess (", mu, "mol P g ", DW^{-1}, ")"))) +
    ggtitle("B)"), nrow = 1
)
