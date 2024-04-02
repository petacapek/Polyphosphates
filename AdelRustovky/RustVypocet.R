#==========================================Pocitani rychlosti rustu===========================================#
#Potrebujes balicek "growthrates"
library(growthrates)
library(dplyr)
library(ggplot2)
library(readODS)
#==========================================Priklad na vzorku c. 1=============================================#
#1. Tohle je prepocitavaci funkce
OxiTopF<-function(filepath, headspace, m, BD, Vmedia, DW, Temperature, ID){
  #Reading raw data from file
  rdata<-read.table(file=c(filepath), skip=17, sep=";", header=F)[,c(2,3)]
  colnames(rdata)<-c("Time", "O2")
  #Converting time to hours
  rdata$Time<-rdata$Time/60
  #Converting diff hPa to diff O2 in umols/g
  rdata$O2diff<-with(rdata, abs(O2)*100*(headspace-m/BD-Vmedia)/(Temperature+273.15)/8.314/(m*DW))
  #Add identifiers
  rdata$ID<-ID
  return(rdata)
}

#a tohle musis do ty funkce dodat
##kde je soubor
Fpath <- c("AdelRustovky/MeasData_23020701.csv")
#Headspace v ml (to mam zmereny z drivejska)
headspace = 301.77
#Vaha celeho vzorku v lahvicce - t.j. pisku i pudy
m = 25 + 0.88
#Objemova hmotnost v g/ml (to mam zmereny z drivejska)
BD = 1.4 
#Objem pridaneho media v ml
Vmedia = 5 
#Susina - se zapocitanim sucheho pisku a mnozstvi pridaneho media
DW = (0.88*0.22 + 25)/(0.88 + 25 + 5)
#Teplota ve stupnich celsia
Temperature = 20
#Oznaceni vzorku 
ID = 1

#2. Prepocitas si pomoci funkce tlak na mnozstvi spotrebovaneho kysliku
R1 <- OxiTopF(filepath = Fpath,
              headspace = headspace, m = m, BD = BD, Vmedia = Vmedia, DW = DW, Temperature = Temperature,
              ID = ID)

#3. Vyneses si graf zavislosti spotrebovaneho kysliku (v umol(O2)/g(DW)) na case (hodiny)
ggplot(R1, aes(Time, O2diff)) + geom_point() + theme_min + 
  xlab("Time (hour)") + ylab(expression(paste(O[2], " consumption (", mu, "mol ", g^{-1}, ")"))) +
  stat_function(fun = function(x){8*exp(-exp(0.2*exp(1)*(20 - x)/8 + 1)) + 8*exp(0.01*(x - 10))})

#4. Na zacatku je zmena tlaku velka, protoze se prudce menila teplota tim jak jsme prenaseli vzorky z laborky do inkubatoru
####Pro nasledne fitovani modelu je lepsi vyhodit nulu na zacatku
R1$O2diff[1:3] <- NA
ODataSmooth <- R1 %>% mutate(Tbins = cut(Time, breaks = seq(0, max(R1$Time), length.out = 35))) %>% 
  group_by(Tbins, ID) %>% summarise(Time = median(Time), O2diff = mean(O2diff))
#Estimating growth rate
GParms <- all_easylinear(O2diff ~  Time | ID, data = subset(ODataSmooth, O2diff > 0))
#Visualizing the  fit
plot(GParms)


#Ted se fituje rustova krivka a pomoci ni se odhaduje rustova rychlost
modelfit2 <- fit_growthmodel(FUN = grow_gompertz3, p = c(y0 = 10, mumax = 1, K = 30, lambda = 10), R1$Time, R1$O2diff)
coef(modelfit2) #mumax je ta rustova ryhlost
summary(modelfit2)

#Muzes se podivat jak hezky ti model nafitoval data
ggplot(R1, aes(Time, O2diff)) + geom_point(cex = 6, pch = 21, fill = "grey") + theme_min + 
  xlab("Time (hour)") + ylab(expression(paste(O[2], " consumption (", mu, "mol ", g^{-1}, ")"))) +
  geom_line(data = as.data.frame(grow_gompertz3(R1$Time, coef(modelfit2))), aes(time, y), lwd = 1.2)

#=============================================================================================================#

#==============================================Vsechny vzorky=================================================#
#Pro lepsi porovnatelnost vysledku budeme vsechny data standardizovat na 1 umol MBC na zacatku inkubace
OxiTopFBiomass<-function(filepath, headspace, m, BD, Vmedia, Temperature, ID, MBCInit){
  #Reading raw data from file
  rdata<-read.table(file=c(filepath), skip=17, sep=";", header=F)[,c(2,3)]
  colnames(rdata)<-c("Time", "O2")
  #Converting time to hours
  rdata$Time<-rdata$Time/60
  #Converting diff hPa to diff O2 in umols/g
  rdata$O2diff<-with(rdata, abs(O2)*100*(headspace-m/BD-Vmedia)/(Temperature+273.15)/8.314/(MBCInit))#
  rdata$O2diff[1:5] <- NA
  #Add identifiers
  rdata$ID<-ID
  #Data smoothing
  rdataSmooth <- rdata %>% mutate(Tbins = cut(Time, breaks = seq(0, max(rdata$Time), length.out = 35))) %>% 
    group_by(Tbins, ID) %>% summarise(Time = median(Time), O2diff = mean(O2diff))
  
  return(rdataSmooth)
}
##kde je soubor
Fpaths <- list.files("AdelRustovky/")
#Headspace v ml (to mam zmereny z drivejska)
headspace = 301.77
#Vaha celych vzorku v lahvicce - t.j. pisku i pudy
m = c(rep(25 + 0.88, 3), rep(25 + 1.07, 3), rep(25 + 1.31, 3), rep(25 + 2.01, 3),
      rep(25 + 0.65, 3), rep(25 + 1.13, 3), rep(25 + 0.6, 3), rep(25 + 0.9, 3),
      rep(25 + 1, 3), rep(25 + 1.62, 3), rep(25 + 1.30, 3), rep(25 + 1.26, 3))
#Pocatecni mikrobialni biomasa
MBCi = read_ods("/mnt/580CBE2464C5F83D/pracovni/data_statistika/AdelMGR/MGRAdelData.ods", 2)[, 17]
MBCi <- c(rep(mean(MBCi[1:3]), 3), rep(mean(MBCi[3:6]), 3), rep(mean(MBCi[7:9]), 3), rep(mean(MBCi[10:12]), 3),
          rep(mean(MBCi[13:15]), 3), rep(mean(MBCi[16:18]), 3), rep(mean(MBCi[19:21]), 3), rep(mean(MBCi[22:24]), 3),
          rep(mean(MBCi[25:27]), 3), rep(mean(MBCi[28:30]), 3), rep(mean(MBCi[31:33]), 3), rep(mean(MBCi[34:36]), 3))
DW = read_ods("/mnt/580CBE2464C5F83D/pracovni/data_statistika/AdelMGR/MGRAdelData.ods", 2)[, 7]
DW = ((m-25)*0.22 + 25)/((m-25) + 25 + 5)
MBCInit = MBCi*(m - 25)*DW
#Objemova hmotnost v g/ml (to mam zmereny z drivejska)
BD = 1.4 
#Objem pridaneho media v ml
Vmedia = 5 
#Susina - tu uz ted nepotrebujeme
##DW = (0.88*0.22 + 25)/(0.88 + 25 + 5)
#Teplota ve stupnich celsia
Temperature = 20
#Oznaceni vzorku 
ID = 1:36

###Sem budu ukladat vysledky:
RustData <- data.frame(Time = numeric(), O2diff = numeric(), ID = numeric())

for(i in 1:36){
  RustData <- rbind(RustData,
                    OxiTopFBiomass(filepath = paste0("AdelRustovky/", Fpaths[i]),
                                   headspace = headspace, m = m[i], BD = BD, Vmedia = Vmedia, Temperature = Temperature,
                                   ID = ID[i], MBCInit = MBCInit[i])) #MBCInit = MBCInit[i]), 
}
#Takhle vypadaji vysledky
ggplot(RustData, aes(Time, O2diff)) + geom_point(cex = 6, pch = 21, fill = "grey") + theme_min + facet_wrap(~ID, scales = "free_y") +
  xlab("Time (hour)") + ylab(expression(paste(O[2], " consumption (", mu, "mol ", mu, mol(MBC)^{-1}, ")")))

RustParametryL <- all_easylinear(O2diff ~ Time | ID, data = RustData[!is.na(RustData$O2diff), ])
coef(RustParametryL)
par(mfrow = c(6, 6))
par(mar = c(2.5, 4, 2, 1))
plot(RustParametryL)

#A ted se jen spocitaji rychlosti rustu
RustParametry <- data.frame(y0 = numeric(), mumax = numeric(), K = numeric(), lambda = numeric(), ID = numeric())
RustData <- as.data.frame(RustData)

for(i in 1:36){
  #Ted se fituje rustova krivka a pomoci ni se odhaduje rustova rychlost
  modelfit <- fit_growthmodel(FUN = grow_gompertz3, p = coef(modelfit2), 
                              as.numeric(RustData[RustData$ID == i, "Time"]), 
                              as.numeric(RustData[RustData$ID == i, "O2diff"]))
  RustParametry <- rbind(RustParametry, 
                         data.frame(y0 = as.numeric(coef(modelfit)[1]), 
                                    mumax = as.numeric(coef(modelfit)[2]), 
                                    K = as.numeric(coef(modelfit)[3]), 
                                    lambda = as.numeric(coef(modelfit)[4]), ID = i))
  
}


#Zobrazime si jestli nekde neco nehapruje
RustPredikce <- data.frame(time = numeric(), y = numeric(), ID = numeric())

for(i in 1:36){
  Gomp <- as.data.frame(grow_gompertz3(as.numeric(RustData[RustData$ID == i, "Time"]), RustParametry[i, 1:4]))
  Gomp$ID <- i
  RustPredikce <- rbind(RustPredikce, Gomp)
}

ggplot(RustData, aes(Time, O2diff)) + geom_point(cex = 6, pch = 21, fill = "grey") + theme_min + facet_wrap(~ID, scales = "free_y") +
  xlab("Time (hour)") + ylab(expression(paste(O[2], " consumption (", mu, "mol ", mu, mol(MBC)^{-1}, ")"))) +
  geom_line(data = RustPredikce, aes(time, y), lwd = 1.2, col = "red")


#Exportujeme vysledky
write.csv(RustParametry, "/mnt/580CBE2464C5F83D/pracovni/data_statistika/AdelMGR/AdelRustovky/RustoveParametry.csv")
