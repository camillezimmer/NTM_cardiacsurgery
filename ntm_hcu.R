setwd("~/Documents/GitHub_local/maddy/QMRAIV")

set.seed(1234)
iters <- 10000

require(EnvStats)

#---- Set parameters as vectors ----

conc_hcu <- vector() #concentration of NTM in HCU
aero_ratio <- vector() #aerosolization ratio
conc_air <- vector() #concentration of NTM in air/aerosols
depo <- vector() #deposition rate
t <- vector() #duration of exposure (hrs)
size_w <- vector() #size of wound (cm^2)
dose <- vector() #dose received (no intervention)
k <- vector() #exponential dose-response parameter
Risk <- vector()

#---- Loop: no intervention ----

for(i in 1:iters)
{
  #conc_hcu[i] <- rlnorm(10000,-2.15526,0.549)
  conc_hcu[i] <- rlnormTrunc(10000,-6.463297,5.918223,max=100)
  #aero_ratio[i] <- rlnorm(100000,16.96646,2.252833)
  #aero_ratio[i] <- rnorm(10000,7.368443,0.978393)
  conc_air[i] <- 0.0001*conc_hcu[i]+14.556
  #conc_air[i] <- conc_hcu[i]*aero_ratio[i]
  depo[i] <- 0.0058*conc_air[i]
  size_w[i] <- runif(10000,7.5,12.5)
  t[i] <- runif(10000,46/60,294/60)
  dose[i] <- size_w[i]*depo[i]*t[i]
  #k[i] <- 0.00000000312 #endpoint: death
  #k[i] <- 0.0000011 #endpoint: lung lesions
  k[i] <- rlnorm(10000,-13.742,0.208) #endpoint: subclinical infection
  Risk[i] <- 1-exp(-k[i]*dose[i])
}

###interventions

dis_lr <- vector() #log-reduction from disinfecting HCU water
conc_hcu_dis <- vector() #concentration in HCU after disinfection
dose_dis <- vector() #dose received (disinfection)
Risk_dis <- vector()

filt_lr <- vector() #log-reduction from filtration
dose_filt <- vector() #dose received (filtration)
dose_flow <- vector() #dose received (laminar flow)
Risk_filt <- vector()
Risk_flow <- vector()

#---- Loop: intervention 1 (disinfection) ----
for(i in 1:iters)
{
  conc_hcu[i]
  dis_lr[i] <-
  conc_hcu_dis[i] <- conc_hcu*dis_lr[i]
  conc_air[i] <- 2005*log(conc_hcu_dis[i])-19822
  depo[i] <- 0.0058*conc_air[i]
  size_w[i] <- 
  t[i] <-
  dose[i] <- size_w[i]*dep[i]*t[i]
  k[i] <- TriRand(0.225,0.491,1.39) 
  
}

#---- Loop: intervention 2 (filter) ----

#Put code here :)

#---- Loop: intervention 3 (laminar flow) ----

#Put code here :)

#---- DOSE RESPONSE ----
# exp.dr <- function(k,dose) 1-exp(-k*dose)
# Risk <- exp.dr(k[i],dose[i]) 
# 
# Risk <- 1-exp(-k[i]*dose[i])

Risk_dis <- exp.dr(k[i],dose_dis)
Risk_filt <- exp.dr(k[i],dose_filt)
Risk_flow <- exp.dr(k[i],dose_flow)

#---- RISK CHARACTERIZATION ----

#INFERENCES 

Risks <- cbind(Risk,Risk_dis,Risk_filt,Risk_flow); colnames(Risks) <- c("None","Disinfection","Filtration","Laminar Flow")
outs <- data.frame(mean(Risk), median(Risk), sd(Risk), min(Risk), quantile(Risk,probs=0.05), quantile(Risk, probs=0.95), max(Risk))
outs_dis <- data.frame(mean(Risk_dis), median(Risk_dis), sd(Risk_dis), min(Risk_dis), quantile(Risk_dis,probs=0.05), quantile(Risk_dis, probs=0.95), max(Risk_dis))
outs_filt <- data.frame(mean(Risk_filt), median(Risk_filt), sd(Risk_filt), min(Risk_filt), quantile(Risk_filt,probs=0.05), quantile(Risk_filt, probs=0.95), max(Risk_filt))
outs_flow <- data.frame(mean(Risk_flow), median(Risk_flow), sd(Risk_flow), min(Risk_flow), quantile(Risk_flow,probs=0.05), quantile(Risk_flow, probs=0.95), max(Risk_flow))
colnames(outs) <- c("Mean", "Median", "Std Dev", "Min", "Lower 95th", "Upper 95th", "Max")
colnames(outs_dis) <- c("Mean", "Median" ,"Std Dev", "Min" ,"Lower 95th", "Upper 95th", "Max")
colnames(outs_filt) <- c("Mean", "Median" ,"Std Dev", "Min" ,"Lower 95th", "Upper 95th", "Max")
colnames(outs_flow) <- c("Mean", "Median" ,"Std Dev", "Min" ,"Lower 95th", "Upper 95th", "Max")
outputs <- rbind(outs,outs_dis,outs_filt,outs_flow)
rownames(outputs) <- c("None","Disinfection","Filtration","Laminar Flow")
outputs

#VISUALIZATION

#Distribution of input parameters
#hist(conc_hcu, main = "", xlab = "Concentration (CFU/mL)", xlim=c(0,1)); rug(conc_hcu, col = "blue")
hist(conc_hcu)
hist(aero_ratio)
hist(conc_air)
hist(depo)
hist(t)
hist(size_w)
hist(dose)
hist(k)

hist(Risk)
hist(Risk_dis)
hist(Risk_filt)
hist(Risk_flow)


RiskNT_hist <- hist(RiskNT, plot=FALSE); RiskNT_P = RiskNT_hist$counts/iters
RiskT_hist <- hist(RiskT, plot=FALSE); RiskT_P = RiskT_hist$counts/iters

png("LR_lg.png", res=300, width = 6, height = 7, units = "in") #change plot name for pop size (sm, med, lg)
par(mfrow=c(3,2))
hist(RiskNT, main="Risk: Raw", xlab = "Probability of infection"); rug(RiskNT, col="#22269E")
hist(RiskT, main= "Risk: Treated", xlab = "Probability of infection"); rug(RiskT, col="#F18A00")
plot(RiskNT_P, main="Density Raw", ylab = "Probability of Risk Occurrence")
plot(RiskT_P, main="Density Treated", ylab = "Probability of Risk Occurrence")
boxplot(Risks, main= "Boxplot by Exposure", ylab= "Risk", col=c("#ebeafe","#fee7fb"))
boxplot(Risks, main= "Boxplot by Exposure - Log", ylab= "Log Risk", log="y", col=c("#ebeafe","#fee7fb"))
dev.off()

#---- SENSITIVITY ANALYSES----
# ---- LOAD PACKAGES ----

if("reshape" %in% rownames(installed.packages())==FALSE){install.packages("reshape", 
                                                                          dependencies = TRUE, repos = "http://cran.us.r-project.org"); require("reshape")}else{require("reshape")}
if("ggplot2" %in% rownames(installed.packages())==FALSE){install.packages("ggplot2", 
                                                                          dependencies = TRUE, repos = "http://cran.us.r-project.org"); require("ggplot2")}else{require("ggplot2")}

# ---- GENERALLY USEFUL VARIABLES ----
old.par <- par( no.readonly = TRUE)

# ---- SCATTERS TO EVALUATE USE OF SPEARMAN CORRELATION METHOD ----  ##NOT WORKING

png("ScatterRiskNT_LR_lg.png", width = 1000, height = 300) #not treated #change plot name for pop size (sm, med, lg)
par(mfrow = c(1,5), cex=1, oma=c(0,0,0,0))
plot(a0,RiskNT, log="y", xlab="Intercept", ylab="Risk Untreated")
lines(lowess(a0,RiskNT), lwd = 4, col = "#de712b")
plot(a1,RiskNT, log="y", xlab="Slope", ylab="")
lines(lowess(a1,RiskNT), lwd = 4, col = "#de712b")
plot(res_err,RiskNT, log="y", xlab="Residual error", ylab="")
lines(lowess(res_err,RiskNT), lwd = 4, col = "#de712b")
plot(decay,RiskNT, log="y", xlab="Decay", ylab="")
lines(lowess(decay, RiskNT), lwd = 4, col = "#de712b")
plot(pop_size,RiskNT, log="y", xlab="Population size", ylab="")
lines(lowess(pop_size,RiskNT), lwd = 4, col = "#de712b")
dev.off()

png("ScatterRiskTreat_LR_lg.png", width = 1000, height = 300) #not treated #change plot name for pop size (sm, med, lg)
par(mfrow = c(1,6), cex=1, oma=c(0,0,0,0))
plot(a0,RiskT, log="y", xlab="Intercept", ylab="Risk Untreated")
lines(lowess(a0,RiskT), lwd = 4, col = "#de712b")
plot(a1,RiskT, log="y", xlab="Slope", ylab="")
lines(lowess(a1,RiskT), lwd = 4, col = "#de712b")
plot(res_err,RiskT, log="y", xlab="Residual error", ylab="")
lines(lowess(res_err,RiskT), lwd = 4, col = "#de712b")
plot(decay,RiskT, log="y", xlab="Decay", ylab="")
lines(lowess(decay, RiskT), lwd = 4, col = "#de712b")
plot(treat,RiskT, log="y", xlab="Water treatment", ylab="")
lines(lowess(treat, RiskT), lwd = 4, col = "#de712b")
plot(pop_size,RiskT, log="y", xlab="Population size", ylab="")
lines(lowess(pop_size,RiskT), lwd = 4, col = "#de712b")
dev.off()

Risk.Vars <- data.frame(pop_size, a0, a1, res_err, decay, pred_excrtn, wp, depth, dist, IR, k) 
colnames(Risk.Vars) <- c("P-size", "Intercept", "Slope", "Residual error", "Decay", "Virus excreted", "WP", "Depth", "Dist", "Ingestion rate", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars) #not treated
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(RiskNT,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
RiskLR.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme(legend.position = "none") + 
  labs(y = expression(bold(~Spearman~rho)), x = expression(bold(~Model~Variable)), 
       title = expression(" "))

Risk.Vars.T <- data.frame(pop_size, a0, a1, res_err, decay, pred_excrtn, wp, depth, dist, IR, k, treat) 
colnames(Risk.Vars.T) <- c("P-size", "Intercept", "Slope", "Residual error", "Decay", "Virus excreted", "WP", "Depth", "Dist", "Ingestion rate", "k", "Treatment")
Sens.Risk.T <- matrix(NA, 1, ncol(Risk.Vars.T)); colnames(Sens.Risk.T) <- colnames(Risk.Vars.T) #treated
for(j in 1:ncol(Sens.Risk.T))
{
  Sens.Risk.T[,j] <- cor(log(Risk.Vars.T[,j],10),log(RiskT,10), method = "spearman", use="complete.obs")
}
Risk.Sens.T <- melt(Sens.Risk.T); colnames(Risk.Sens.T) <- c("dummy","Variable","Rho")
RiskLR.Sens.T.PLOT <- ggplot(Risk.Sens.T,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme(legend.position = "none") + 
  labs(y = expression(bold(~Spearman~rho)), x = expression(bold(~Model~Variable)), 
       title = expression(" "))

ggsave("Sensitivity_NT_LR_lg.png", RiskLR.Sens.PLOT) #change plot name for pop size (sm, med, lg)
ggsave("Sensitivity_T_LR_lg.png", RiskLR.Sens.T.PLOT) #change plot name for pop size (sm, med, lg)

Risk.sens <- ggarrange(RiskLR.Sens.PLOT, RiskLR.Sens.T.PLOT,
                       labels = c("Untreated", "Treated"), ncol = 2, nrow = 1)
png("Sensitivity_LR_lg.png", res=300, width = 9, height = 7, units = "in") #change plot name for pop size (sm, med, lg)
print(Risk.sens)
#dev.off()


#---- Daily and Annualized Risks ----
# DailyRisks <- Risks
# 
# d2A_n=2*365 #2L/day 
# 
# daily2annualRisks <- function(DailyRisks, numPerDay=1, numPerYear=d2A_n){
#   calcAnnualRisk <- function(DailyRiskCol){
#     sampledRisks <- rep(sample(DailyRiskCol,numPerYear,replace=TRUE),numPerDay)
#     AnnualRisk <- 1-prod(1-sampledRisks)
#     return(AnnualRisk)
#   }
#   AnnualRisks <- c()
#   for(i in 1:length(DailyRisks))
#     AnnualRisks[i] <- calcAnnualRisk(DailyRisks) 
#   return(AnnualRisks)
# }
# 
# AnnualRisksNT <- daily2annualRisks(RiskNT)
# AnnualRisksT <-  daily2annualRisks(RiskT)