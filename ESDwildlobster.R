# Written by: Annika Bell
# last Edited: 2024 August 14

#General directory set up ----

#install packages and library them for all packages needed
library(ggplot2)
install.packages("viridis")
library(viridis)
install.packages("dr4pl")
library(dr4pl)
library(tidyverse)
install.packages("lme4")
library(lme4)
library(tidyverse)
library(tidyr)
# upload data as csv and add to bigelow R folder
e <- read.csv('ESD_trial_data.csv')
e$LobNo <- as.factor(e$LobNo)

# Glucose data manipulation ----
# remove empty ending rows to make a glucose processing data set!
enumb <- subset(e, select = c(LobNo, Location, ESD, Glucose,Brix,Weight,CL,TailL))
gdata <- na.omit(enumb)
#make a molt scale variable from brix based on  
gdata$MoltScale <- ifelse(gdata$Brix<8.5, "soft", 
                          ifelse(gdata$Brix>13, "pre", "inter"))

#general linear model to test ESD on glucose concentration result.
gdata$ESD <- as.factor(gdata$ESD)
glinear_model <- lm(Glucose ~ ESD, data=gdata) 
plot(glinear_model)
summary(glinear_model)

#summarize glucose data frame to get mean glucose amounts by ESD severity category and location
sumgdata <- gdata %>%
  group_by(ESD, Location) %>%
  summarize(mean=mean(Glucose), sd=sd(Glucose), se=sd(Glucose)/sqrt(length(Glucose)), n=length(Glucose))

#these are the graph labels I want to use for my ESD scale
labels <-c("None", "Mild", "Moderate", "Severe")

# graphing Glucose concentration of American lobsters across ESD status
sumgdata$ESD <- as.factor(sumgdata$ESD)
sumgdata$ESD1 <-as.integer(sumgdata$ESD)
gsimple <- ggplot(sumgdata, aes(ESD,mean, fill=ESD1))+geom_col()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(y="Glucose Concentration (mg/mL)", x="ESD Lesions")+
  scale_fill_viridis_c(option="rocket", direction = -1, begin = 0.2, end = 0.8, guide = "none")+scale_x_discrete(labels = c("0" = "None", "1" = "Mild", "2" = "Moderate", "3" ="Severe"))

#glucose brix graphs to see how brix impacts glucose results
sumgbw <- gdata %>%
  group_by(ESD, Brix, Weight, Location) %>%
  summarize(mean=mean(Glucose), sd=sd(Glucose), se=sd(Glucose)/sqrt(length(Glucose)), n=length(Glucose))

sumgbw$ESD <-as.integer(sumgbw$ESD)

gsimpleb <- ggplot(sumgbw, aes(Brix,mean, col=ESD))+geom_point()+ scale_color_viridis(option="rocket", direction = -1, begin = 0.2, end = 0.8)+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(y="Glucose Concentration (mg/mL)", x="Brix (%)")

#graph of weight, esd and glucose readings
gsimplew <- ggplot(sumgbw, aes(Weight,mean, col=ESD))+geom_point()+ scale_color_viridis(option="rocket", direction = -1, begin = 0.2, end = 0.8)+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(y="Glucose Concentration (mg/mL)", x="Weight (g)")

#Glycogen hypothesized graph for poster NOT DATA BASED predicted----
glyhypo <- data.frame(
  disease = c("None", "Mild", "Moderate", "Severe"),
  amount = c(80,75, 50, 20),
  col_D = c(0,1,2,3)
)

Glycogenhypograph <- ggplot(glyhypo, aes(x=fct_inorder(disease),amount, fill=col_D))+geom_col()+
  labs(y="Hypothesized Glycogen Content", x="ESD Lesions")+
  scale_fill_viridis(option="rocket", direction = -1, begin = 0.2, end = 0.8, guide = "none")+
  theme_classic()+theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),axis.ticks.y = element_blank(),axis.text.y = element_blank())

#hepatopancreas data and graphs ----

#upload and subset data frame for HP
# I create separate data frames for each part because of the differing NAs to avoid data loss.
hpd <- read.csv('ESD_trial_data.csv')
hnumber <- subset(hpd, select = c(LobNo, Location, ESD, Hpper,Brix,Weight,CL,TailL))


hdata <- na.omit(hnumber)

#HP and ESD standard graph
graphhp <-  ggplot(hdata, aes(ESD, Hpper, color=ESD))+geom_point(size=2)+theme_classic()+
  scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+
  labs(y="Hepatopancreas (%)", x="ESD Lesions")+
  scale_color_viridis_c(option="rocket", direction = -1, begin = 0.2, end = 0.8, guide="none")+
  scale_x_continuous(label=labels)

#HP graph with ESD and Brix
graphhpb <-  ggplot(hdata, aes(Hpper, Brix, color=ESD))+geom_point(position = "jitter")+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(x="Hepatopancreas (%)", y="Brix (%)")+
  scale_color_viridis_c(option="rocket", direction = -1, begin = 0.2, end = 0.8)

#HP graph with ESD and Molt scale
hdata$MoltScale <- ifelse(hdata$Brix<8.5, "soft", 
                          ifelse(hdata$Brix>13, "pre", "inter"))
graphhpm <-  ggplot(hdata, aes(MoltScale, Hpper, color=ESD))+geom_point(position = "jitter")+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(y="Hepatopancreas (%)", x="Molt Status")+
  scale_color_viridis_c(option="rocket", direction = -1, begin = 0.2, end = 0.8)+
  scale_x_discrete(limits=c("soft", "inter", "pre"), labels=c("soft" ="Newly Molted", "inter" = "Intermolt", "pre"= "About to Molt"))

#hemocyte data and graphs ----
#total hemocyte counts
hemo <- read.csv('ESD_trial_data.csv')

hemot <- subset(hemo, select = c(LobNo, Location, Brix,tcave1, ESD))
tcdata <- na.omit(hemot)

graphtcall <-  ggplot(tcdata, aes(ESD, tcave1))+geom_point()+theme_classic()

tcfilterdata <- tcdata %>%
  group_by(ESD) %>%
  summarize(mean=mean(tcave1), sd=sd(tcave1), se=sd(tcave1)/sqrt(length(tcave1)),n=length(tcave1))

tcfilterdata$ESD <- as.factor(tcfilterdata$ESD)
tcfilterdata$ESD1 <-as.integer(tcfilterdata$ESD)
#total hemocyte counts graph by ESD
graphtc <-  ggplot(tcfilterdata, aes(ESD, mean, fill=ESD1))+geom_col()+theme_classic()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs(y="Total Hemocyte Counts (cells / mL)", x="ESD Lesions")+
  scale_fill_viridis_c(option="rocket", direction = -1, begin = 0.2, end = 0.8, guide = "none")+scale_x_discrete(labels = c("0" = "None", "1" = "Mild", "2" = "Moderate", "3" ="Severe"))

#generalized linear model to see significance of variables on total hemocyte counts
tcdata$ESD <- as.factor(tcdata$ESD) 
tclinear_model <- lm(tcave1 ~ Brix, data=tcdata) 
summary(tclinear_model)


#Differential counts for hemocytes
hemodiff <- read.csv('ESD_trial_data.csv')
hemod <- subset(hemodiff, select = c(LobNo, Location, Brix, ESD,Hdiff,Sdiff,Ldiff))
hddata <- na.omit(hemod)

hddata$MoltScale <- ifelse(hddata$Brix<8.5, "soft", 
                           ifelse(hddata$Brix>13, "pre", "inter"))

#molt stage vs esd graph
graphmesd <-  ggplot(hddata, aes(ESD, MoltScale))+geom_point(position="jitter")+theme_classic()

#generalized linear model to see significance of variables on hyaline hemocytes readings
hddata$ESD <- as.factor(hddata$ESD) 
hlinear_model1 <- lm(Hdiff ~ MoltScale, data=hddata) 
plot(hlinear_model1)
summary(hlinear_model1)

# this makes my variables in the correct order
hddata$aHdiff <- hddata$Hdiff
hddata$bSdiff <- hddata$Sdiff

# pivot longer to make cell types into one variable
hddatasort <-pivot_longer(hddata,
                      cols = c(aHdiff, bSdiff, Ldiff),
                      names_to = "celltype",
                      values_to = "count",
                      names_transform = list(key = forcats::fct_inorder)) 

hddatasort <- hddatasort %>%
  group_by(ESD,celltype) %>%
  summarize(mean=mean(count), sd=sd(count), n=length(count))

#This the differential hemocyte count graph against ESD status
hddatasort$ESD <- as.factor(hddatasort$ESD)
cellsimple <- ggplot(hddatasort, aes(ESD, mean, fill=celltype))+geom_col(position = "dodge")+geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme_classic()+scale_y_continuous(expand = c(0,0,0,0.5),limits=c(0,NA))+labs( x="ESD Lesions", y="Differential Hemocyte Counts (%)", fill= "Cell Type")+
  scale_fill_viridis(option="mako", discrete=TRUE, direction = -1, begin = 0.2, end = 0.8, labels=c("Hyaline", "Small Granulocyte", "Large Granulocyte")) +scale_x_discrete(labels = c("0" = "None", "1" = "Mild", "2" = "Moderate", "3" ="Severe"))

