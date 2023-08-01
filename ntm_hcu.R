setwd("~/Documents/GitHub_local/maddy/QMRAIV")

if("triangle" %in% rownames(installed.packages())==FALSE){install.packages("triangle");require(triangle)}else(require(triangle))
source("TriRand.r")

require(stats)
TriRand <- function(minValue, likeValue, maxValue)
{
  z = runif(1); .Random.seed[1:1];
  t = sqrt(z*(maxValue-minValue)*(likeValue-minValue))+minValue
  tt = maxValue-sqrt((1-z)*(maxValue-minValue)*(maxValue-likeValue))
  if (tt < likeValue) {return(t)} else{return(tt)}
}

set.seed(1234)
iters <- 10000

#set parameters as vectors

conc_hcu <- vector() #concentration of NTM in HCU

conc_aer <- vector() #concentration of NTM in aerosols
t <- vector() #duration of exposure
size_w <- vector() #size of wound
dis_lr <- vector() #log-reduction from disinfection
filt_lr <- vector() #log-reduction from filtration
k <- vector() #exponential dose-response parameter
dose <- vector() #dose received (no intervention)
dose_dis <- vector() #dose received (disinfection)
dose_filt <- vector() #dose received (filtration)
dose_flow <- vector() #dose received (laminar flow)
Risk <- vector()
Risk_dis <- vector()
Risk_filt <- vector()
Risk_flow <- vector()

#Loop: no intervention

for(i in 1:iters)
{
  pop_size[i] <- runif(1,50,100) #size of population vaccinated; VARIABLE (small, med, large populations)
  a0[i] <- TriRand(5.56,5.98,6.39) #intercept
  a1[i] <- TriRand(0.83,0.96,1.14) #slope
  res_err[i] <- TriRand(0.33,0.39,0.46) #random error
  pred_excrtn[i] <- ((a1[i]*log(pop_size[i],10))+a0[i]+res_err[i]) #predicted excretion via linear regression
  pred_excrtn[i] <- 10^(pred_excrtn[i]) #log10 scale 
  wp[i] <- TriRand(59.3,65.8,73.0) #wetted perimeter
  depth[i] <- TriRand(8,13,18) #depth  
  dist[i] <- TriRand(100,367010,733920) #distance between community A and B
  river[i] <- wp[i]*depth[i]*dist[i] #volume of river
  Conc_vax[i] <- pred_excrtn[i]/river[i] #concentration of total poliovirus in river
  decay[i] <- TriRand(0.0855, 0.1411, 0.2581) #decay rate
  Conc_cap[i] <- Conc_vax[i]*10^(-decay[i]) #concentration of untreated captured water
  IR[i] <- TriRand(0.5,0.75,1) #ingestion rate
  dose_notreat[i] <- Conc_cap[i]*IR[i] 
  k[i] <- TriRand(0.225,0.491,1.39) 
  
}

#Loop: intervention 1 (disinfection)
for(i in 1:iters)
{
  pop_size[i] <- runif(1,50,100) #VARIABLE (small, med, large populations)
  a0[i] <- TriRand(5.56,5.98,6.39) 
  a1[i] <- TriRand(0.83,0.96,1.14) 
  res_err[i] <- TriRand(0.33,0.39,0.46) 
  pred_excrtn[i] <- ((a1[i]*log(pop_size[i],10))+a0[i]+res_err[i])
  pred_excrtn[i] <- 10^(pred_excrtn[i])
  wp[i] <- TriRand(59.3,65.8,73.0)
  depth[i] <- TriRand(8,13,18)   
  dist[i] <- TriRand(100,367010,733920)  
  river[i] <- wp[i]*depth[i]*dist[i]
  Conc_vax[i] <- pred_excrtn[i]/river[i]
  decay[i] <- TriRand(0.0855, 0.1411, 0.2581)
  Conc_cap[i] <- Conc_vax[i]*10^(-decay[i])
  treat[i] <- runif(1,1,2.4) #log reduction from treatment        
  Conc_treat[i] <- Conc_cap[i]*10^(-treat[i]) #concentration of treated captured water
  IR[i] <- TriRand(0.5,0.75,1)  
  dose_treat[i] <- Conc_treat[i]*IR[i]
  k[i] <- TriRand(0.225,0.491,1.39)
}

#Loop: intervention 2 (filter)
#Loop: intervention 1 (laminar flow)

#Dose response 
exp.dr <- function(k,dose) 1-exp(-k*dose)
Risk <- exp.dr(k[i],dose) 
Risk_dis <- exp.dr(k[i],dose_dis)
Risk_filt <- exp.dr(k[i],dose_filt)
Risk_flow <- exp.dr(k[i],dose_flow)

#RISK CHARACTERIZATION

#Inferences

Risks <- cbind(RiskNT,RiskT); colnames(Risks) <- c("Raw", "Treated")
outsNT <- data.frame(mean(RiskNT), median(RiskNT), sd(RiskNT), min(RiskNT), quantile(RiskNT,probs=0.05), quantile(RiskNT, probs=0.95), max(RiskNT))
outsT <- data.frame(mean(RiskT), median(RiskT), sd(RiskT), min(RiskT), quantile(RiskT,probs=0.05), quantile(RiskT, probs=0.95), max(RiskT))
colnames(outsNT) <- c("Mean", "Median", "Std Dev", "Min", "Lower 95th", "Upper 95th", "Max")
colnames(outsT) <- c("Mean", "Median" ,"Std Dev", "Min" ,"Lower 95th", "Upper 95th", "Max")
outputs <- rbind(outsNT,outsT)
rownames(outputs) <- c("Untreated", "Treated")
outputs

#Visualization

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

#Sensitivity analysis
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



DailyRisks <- Risks

d2A_n=2*365 #2L/day 

daily2annualRisks <- function(DailyRisks, numPerDay=1, numPerYear=d2A_n){
  calcAnnualRisk <- function(DailyRiskCol){
    sampledRisks <- rep(sample(DailyRiskCol,numPerYear,replace=TRUE),numPerDay)
    AnnualRisk <- 1-prod(1-sampledRisks)
    return(AnnualRisk)
  }
  AnnualRisks <- c()
  for(i in 1:length(DailyRisks))
    AnnualRisks[i] <- calcAnnualRisk(DailyRisks) 
  return(AnnualRisks)
}

AnnualRisksNT <- daily2annualRisks(RiskNT)
AnnualRisksT <-  daily2annualRisks(RiskT)