setwd("~/Activité Professionnelle/CNRS Toulouse (2018-2019)/Publications/Mismatch Energétique")

library(alr3)
library(Deriv) 
library(ggplot2)
library(hier.part)
library(lmodel2)
library(lmtest)
library(nlme)

###############################################################
###############################################################
##### ENERGETIC MISMATCH BETWEEN METABOLISM AND INGESTION #####
###############################################################
###############################################################

########################
### Data preparation ###
########################

# Loading the dataframe
Tab=read.table("Data_RN.txt", h=T, dec=".")
names(Tab)
summary(Tab)
Tab=na.omit(Tab)

# Consider variables as numeric or factor
Tab$Indiv=as.factor(Tab$Indiv)

# Remove negative values
Tab=Tab[Tab$Respi>0,]                              
min(Tab$Respi)
Tab=Tab[Tab$Nutri>0,]              
min(Tab$Nutri)

# Remove outliers
Tab=Tab[!Tab$Indiv=="12",]
Tab=Tab[!Tab$Indiv=="30",]
Tab=Tab[!Tab$Indiv=="76",]
Tab=Tab[!Tab$Indiv=="78",]

# Check data normality
shapiro.test(Tab$Respi)
qqnorm(Tab$Respi)
qqline(Tab$Respi)
hist(Tab$Respi)

shapiro.test(Tab$Nutri)
qqnorm(Tab$Nutri)
qqline(Tab$Nutri)
hist(Tab$Nutri)

# Check data homoscedasticity
bartlett.test(Tab$Respi~Tab$Temp)
bartlett.test(Tab$Nutri~Tab$Temp)

# Center and convert temperature
Tab$InvTem=1/((Tab$Temp+273.15)*8.62*10^-5)
Tab$InvTemC=(Tab$InvTem-mean(Tab$InvTem))


######################
### Metabolic rate ###
######################

# Check mass-temperature interaction
ModelR0=lm(log(Respi)~log(Masse)*InvTemC, data=Tab)
summary(ModelR0)

# Regression models
ModelLM1=lm(log(Respi)~log(Masse)+InvTemC, data=Tab)
ModelPOLY1=glm(log(Respi)~log(Masse)+poly(InvTemC, degree=2, raw=F), data=Tab)

# Model comparison
anova(ModelLM1,ModelPOLY1)
AIC(ModelLM1,ModelPOLY1)
lrtest(ModelLM1,ModelPOLY1)

# Linear model parameters
ModelR1=lm(log(Respi)~log(Masse)+InvTemC, data=Tab)
summary(ModelR1)
confint(ModelR1)
anova(ModelR1)
coef(ModelR1)

# Quadratic model parameters
ModelR2=glm(log(Respi)~log(Masse)+poly(InvTemC, degree=2, raw=T), data=Tab)
summary(ModelR2)
confint(ModelR2)
anova(ModelR2)
coef(ModelR2)

# Check residuals normality
shapiro.test(ModelR2$residuals)
qqnorm(ModelR2$residuals)
qqline(ModelR2$residuals)
hist(ModelR2$residuals)

# Visualize individual residuals
plot(na.omit(Tab$Temp), residuals(ModelR2), ylab="Residuals", xlab="Values", pch=16)     
ResR=data.frame(Tab[,c(1,2,4,5)], Resid=residuals(ModelR2))
ResR[abs(ResR$Resid) > 1,]

# Partition variance
Variables=cbind(Tab$Masse, InvTemC); colnames(Variables)=c("Masse","InvTemC")
hier.part(Tab$Respi, Variables, fam="gaussian", gof="logLik")

# Plot inverse relationship
PlotM=ggplot(data=Tab, aes(x=InvTemC, y=log(Respi))) +
  geom_smooth(method="lm", formula=y~poly(x,2), color="black", linetype="solid", size=1, se=F) +     
  geom_point(stat="identity", color="turquoise3", size=3) + xlim(-1.2,1.2) + ylim(1.0,5.0) +
  ylab(expression('Ln Routine metabolic rate'~'('*µg~C~day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

# Plot raw relationship
ggplot(data=Tab, aes(x=Temp, y=Respi)) +
  geom_smooth(method="lm", formula=exp(log(y))~exp(poly(1/((x+273.15)*8.62*10^(-5)),2)), color="black", linetype="solid", size=1, se=F) +
  geom_point(stat="identity", color="turquoise3", size=3) + xlim(5,21) + ylim(0,80) + 
  ylab(expression('Routine metabolic rate'~'('*µg~C~day^-1*')')) +  xlab(expression('Temperature (°C)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=20)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


######################
### Ingestion rate ###
######################

# Check mass-temperature interaction
ModelI0=lm(log(Nutri)~log(Masse)*InvTemC, data=Tab)
summary(ModelI0)

# Regression models
ModelLM2=lm(log(Nutri)~log(Masse)+InvTemC, data=Tab)
ModelPOLY2=glm(log(Nutri)~log(Masse)+poly(InvTemC, degree=2, raw=F), data=Tab)

# Model comparison
anova(ModelLM2,ModelPOLY2)
AIC(ModelLM2,ModelPOLY2)
lrtest(ModelLM2,ModelPOLY2)

# Linear model parameters
ModelI1=lm(log(Nutri)~log(Masse)+InvTemC, data=Tab)
summary(ModelI1)
confint(ModelI1)
anova(ModelI1)
coef(ModelI1)

# Quadratic model parameters
ModelI2=lm(log(Nutri)~log(Masse)+poly(InvTemC, degree=2, raw=T), data=Tab)
summary(ModelI2)
confint(ModelI2)
anova(ModelI2)
coef(ModelI2)

# Check residuals normality
shapiro.test(ModelI2$residuals)
qqnorm(ModelI2$residuals)
qqline(ModelI2$residuals)
hist(ModelI2$residuals)

# Visualize individual residuals
plot(na.omit(Tab$Temp), residuals(ModelI2), ylab="Residuals", xlab="Values", pch=16)     
ResI=data.frame(Tab[,c(1,2,4,5)], Resid=residuals(ModelI2))
ResI[abs(ResI$Resid) > 1,]

# Partition variance
Variables=cbind(Tab$Masse, InvTemC); colnames(Variables)=c("Masse","InvTemC")
hier.part(Tab$Nutri, Variables, fam="gaussian", gof="logLik")

# Plot inverse relationship
PlotI=ggplot(data=Tab, aes(x=InvTemC, y=log(Nutri))) +
  geom_smooth(method="lm", formula=y~poly(x,2), color="black", linetype="solid", size=1, se=F) + 
  geom_point(stat="identity", color="tomato1", size=3) + xlim(-1.2,1.2) + ylim(3.0,8.0) +
  ylab(expression('Ln Ingestion rate'~'('*µg~C~day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

# Plot raw relationship
ggplot(data=Tab, aes(x=Temp, y=Nutri)) +
  geom_smooth(method="lm", formula=exp(log(y))~exp(poly(1/((x+273.15)*8.62*10^(-5)),2)), color="black", linetype="solid", size=1, se=F) +
  geom_point(stat="identity", color="tomato1", size=3) + xlim(5,21) + ylim(0,1200) +
  ylab(expression('Ingestion rate'~'('*µg~C~day^-1*')')) + xlab(expression('Temperature (°C)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=20)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


###########################################
### Conversion respiration to ingestion ###
###########################################

tiff('Metabolic and Ingestion Rates.tiff', units="in", width=15, height=8, res=300)
Panel=plot_grid(PlotM, PlotI, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression(atop('Normalized inverse temperature to the reference', paste('temperature'~'(1/kT normalized at 12.5 °C in '*eV^-1*')'))), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()

# Extract conversion parameters
ModelRI1=lmodel2(Nutri~Respi, data=Tab, nperm=99)
ModelRI1

ModelRI2=lm(Nutri~Respi, data=Tab)
summary(ModelRI2)
confint(ModelRI2)
abline(ModelRI2)
anova(ModelRI2)

# Plot the raw relationship
ggplot(data=Tab, aes(x=Respi, y=Nutri)) +
  geom_abline(intercept=346.292, slope=9.527, color="black", linetype="dotted", size=1) +
  geom_point(stat="identity", color="chartreuse3", size=3) + xlim(0,80) + ylim(0,1200) + 
  xlab(expression(atop('Routine Metabolic rate', paste('('*µg~C~day^-1*')')))) +
  ylab(expression(atop('Ingestion rate', paste('('*µg~C~day^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=20)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


########################################################
### Energetic efficiency vs temperature relationship ###
########################################################

# Calculate energetic efficiency
Tab$IngE=Tab$Nutri/Tab$Respi

# Convert energetic efficiency
AssE=0.30
Tab$EneE=Tab$IngE*AssE

# Check residuals normality
shapiro.test(Tab$EneE)
qqnorm(Tab$EneE)
qqline(Tab$EneE)
hist(Tab$EneE)

# Linear model parameters
ModelE=lm(EneE~Temp, data=Tab)
summary(ModelE)
confint(ModelE)
anova(ModelE)
coef(ModelE)

# Plot raw relationship
tiff('Energetic Efficiency.tiff', units="in", width=8, height=8, res=300)
ggplot(data=Tab, aes(x=Temp, y=EneE)) +
  geom_smooth(method="lm", formula=y~x, color="black", linetype="solid", size=1, se=F) +
  geom_point(stat="identity", color="chartreuse3", size=3) + xlim(5,22) + ylim(0,30) + 
  ylab(expression('Energetic efficiency')) + xlab(expression('Temperature (°C)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=20)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
dev.off()
