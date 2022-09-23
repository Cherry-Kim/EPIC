#OxBS Array Processing
library(sesame)
library(colorRamps)
library(ggplot2)
library(data.table)
library(dplyr)
library(gplots)
library(grDevices)
library(reshape2)
library(tidyverse)
library(minfi)
library(RColorBrewer)
library(ENmix)

#STEP2. Make a signal summary dataset for all the IDAT files and run SeSAMe to generate and normalize ?values for each sample
ssets <- lapply(searchIDATprefixes("./"),readIDATpair)
OxBSbetas <- openSesame(ssets)
#To make a bean plot of the β-value distributions from minfi:
ox <- OxBSbetas[,c(1,9,3,11,5,13,7,15)]
bs <- OxBSbetas[,c(2,10,4,12,6,14,8,16)]
oxbs <- cbind(ox,bs)
groups <- c("groupa","groupa","groupb","groupb","groupc","groupc","groupd","groupd","groupa","groupa","groupb","groupb","groupc","groupc","groupd","groupd")
par(mar=c(8,9,7,5))
densityBeanPlot(oxbs, main = "", sampGroups=groups)

#STEP3. Transform the β-value data matrix into a data frame and remove all probes that do not have a β-value for all samples queried.
OxBSbetas_df <- data.frame(OxBSbetas)
OxBSbetas_df <- OxBSbetas_df[complete.cases(OxBSbetas_df),]

#STEP4. Calculate true 5hmCβ-values by subtracting the OxBS array β-values from the BS array β-values for each individual sample.
OxBSbetas_df <- mutate(OxBSbetas_df, ctrl_5hmC= ctrl_BS - ctrl_OX)
OxBSbetas_df <- mutate(OxBSbetas_df, ctrl2_5hmC= ctrl_BS2 - ctrl_OX2)
OxBSbetas_df <- mutate(OxBSbetas_df, Enza_5hmC= Enza_BS - Enza_OX)
OxBSbetas_df <- mutate(OxBSbetas_df, Enza2_5hmC= BS2 - Enza_OX2)
plot(density(OxBSbetas_df$ctrl_5hmC), col = "#190B28", lty =2, lwd = 2,main="", xlim=c(-1,1),ylim=c(0,55))
lines(density(OxBSbetas_df$ctrl_BS), col = "#190B28", lty =1, lwd = 2)
lines(density(OxBSbetas_df$ctrl_OX), col = "#190B28", lty =3, lwd = 2)
lines(density(OxBSbetas_df$ctrl2_5hmC), col = "#2c0e5e", lty =2, lwd = )
lines(density(OxBSbetas_df$ctrl_BS2), col = "#2c0e5e", lty =1, lwd = 2)
lines(density(OxBSbetas_df$ctrl_OX2), col = "#2c0e5e", lty =3, lwd = 2)
lines(density(OxBSbetas_df$Enza_5hmC), col = "#0073C2FF", lty =2, lwd = 2)
lines(density(OxBSbetas_df$Enza_BS), col = "#0073C2FF", lty =1, lwd = 2)
lines(density(OxBSbetas_df$Enza_OX), col = "#0073C2FF", lty =3, lwd = 2)
lines(density(OxBSbetas_df$Enza2_5hmC), col = "#0066FF", lty =2, lwd = 2)
lines(density(OxBSbetas_df$Enza_BS2), col = "#0066FF", lty =1, lwd = 2)
lines(density(OxBSbetas_df$Enza_OX2), col = "#0066FF", lty =3, lwd = 2)
legend("topright",legend=c("5mC+5hmC","5mC","5hmC"),lty=c(1,3,2), lwd = 2,box.lty=0)
text(-0.8,52, labels="Control",col="#190B28")
text(-0.8,50, labels="Control2",col="#2c0e5e")
text(-0.8,48, labels="Enza",col="#0073C2FF")
text(-0.8,46, labels="Enza_2",col="#0066FF")

#STEP5. Determine the number of CpG probes with a 5hmCβ-values above and below 0.
ctrl_5hmC_above <- subset(OxBSbetas_df, ctrl_5hmC> 0) 
ctrl_5hmC_below <- subset(OxBSbetas_df, ctrl_5hmC< 0)   
Enza_5hmC_above <- subset(OxBSbetas_df, Enza_5hmC> 0) 
Enza_5hmC_below <- subset(OxBSbetas_df, Enza_5hmC< 0)   
data1 <- read.table(text="              ctrl  Enza  
 β-val>0 1000 1000  
 β-val<0 1000 1000", header=TRUE)
data2=as.matrix(data1)
b<-barplot(data2, beside= TRUE, col=c("black","cornsilk4"), ylim=c(0,170000), ylab="#CpGs",cex.axis=0.7,cex.names=0.7)
legend("topleft",legend=rownames(data2),box.lty=0)
tx2 <- data2
text(b,tx2+10, as.character(tx2),pos = 3, cex = 0.5, col = "black")

#STEP7. To correct for the number of 5hmC?values below 0, use the OxBS-MLE command from ENmix [53]. First, isolate the ?values for BS array and then isolate the ?values for OxBS array.
OxBSbetas2 <- OxBSbetas
colnames(OxBSbetas2) <- c("ctrl","ctrl","Enza","Enza","ctrl2","ctrl2","Enza2","Enza2")
beta.BS <- OxBSbetas2[,c(2,10,4,12,6,14,8,16)]
beta.oxBS <- OxBSbetas2[,c(1,9,3,11,5,13,7,15)]
#STEP8. Next, isolate the intensity values independently for both BS array and the OxBS array. A critical note is that all samples must remain in the same order and be named the same thing between BS array and OxBS array
ctrl <- totalIntensities(ssets$'ctrl')
ctrl2 <- totalIntensities(ssets$'ctrl2')
Enza <- totalIntensities(ssets$'Enza')
Enza2 <- totalIntensities(ssets$'Enza2')
N.BS <- cbind(ctrl,ctrl2,Enza,Enza2)
N.BS <- N.BS[order(row.names(N.BS)),]
ctrl <- totalIntensities(ssets$'ctrl')
ctrl2 <- totalIntensities(ssets$'ctrl2')
Enza <- totalIntensities(ssets$'Enza')
Enza2 <- totalIntensities(ssets$'Enza2')
N.oxBS <- cbind(ctrl,ctrl2,Enza,Enza2)
N.oxBS <- N.oxBS[order(row.names(N.oxBS)),]
#STEP9. Using the isolated β-values and intensity values from above, run OxBS-MLE to recalculate 5mC and 5hmCβ-values.
OxBS.EN <- oxBS.MLE(beta.BS, beta.oxBS, N.BS, N.oxBS)
OxBS.df <- data.frame(OxBS.EN)
colnames(OxBS.df) <- c("5mC.ctrl","5mC.ctrl2","5mC.Enza","5mC.Enza2","5mC.Enza2.CB839","5mC.Enza2.CB8392","5mC.Enza2.TRC105","5mC.Enza2.TRC1052","5hmC.ctrl","5hmC.ctrl2","5hmC.Enza","5hmC.Enza2","5hmC.Enza2.CB839","5hmC.Enza2.CB8392","5hmC.Enza2.TRC105","5hmC.Enza2.TRC1052")
OxBS.df <- data.matrix(OxBS.df)
hmC <- OxBS.df[,c(9:16)]
mC <- OxBS.df[,c(1:8)]
hmcmc <- cbind(hmC,mC)
groups <- c("groupa","groupa","groupb","groupb","groupc","groupc","groupd","groupd","groupa","groupa","groupb","groupb","groupc","groupc","groupd","groupd")
par(mar=c(8,9,7,5))
densityBeanPlot(hmcmc, main = "", sampGroups=groups)

OxBS.df2 <- OxBS.df[complete.cases(OxBS.df),]
OxBS.df2 <- cbind(rownames(OxBS.df2), data.frame(OxBS.df2), row.names = NULL)
#STEP8. Calculate β-values for each drug treatment relative to each other for both 5mC and 5hmC
OxBS.df2 <- mutate(OxBS.df2, Enza.5hmC.db = (X5hmC.Enza+X5hmC.Enza2) - (X5hmC.ctrl+X5hmC.ctrl2))
OxBS.df2 <- mutate(OxBS.df2, Enza.5mC.db = (X5mC.Enza+X5mC.Enza2) - (X5mC.ctrl+X5mC.ctrl2))

#STEP9. Using the calculated βvalues, define the direction of change for each modification or note if the change is not significant using the following criteria
OxBS.df2 <- mutate(OxBS.df2, EnzavsCtrl_5hmC_direction =
ifelse(Enza.5hmC.db >= 0.1, "Up",
ifelse(Enza.5hmC.db <=-0.1, "Down", "NotSig")))
OxBS.df2 <- mutate(OxBS.df2, EnzavsCtrl_5mC_direction =
ifelse(Enza.5mC.db >= 0.2, "Up",
ifelse(Enza.5mC.db <= -0.2, "Down", "NotSig")))

OxBS.df2 <- mutate(OxBS.df2, EnzavsCtrl.states =
ifelse(EnzavsCtrl_5hmC_direction == "Up" & EnzavsCtrl_5mC_direction=="NotSig", "State1",
ifelse(EnzavsCtrl_5hmC_direction == "NotSig" & EnzavsCtrl_5mC_direction== "Up", "State2",
ifelse(EnzavsCtrl_5hmC_direction == "Up" & EnzavsCtrl_5mC_direction == "Up","State3",
ifelse(EnzavsCtrl_5hmC_direction == "Up" & EnzavsCtrl_5mC_direction == "Down","State4",
ifelse(EnzavsCtrl_5hmC_direction == "Down" & EnzavsCtrl_5mC_direction == "Up", "State5",
ifelse(EnzavsCtrl_5hmC_direction == "Down" & EnzavsCtrl_5mC_direction == "Down","State6",
ifelse(EnzavsCtrl_5hmC_direction == "NotSig" & EnzavsCtrl_5mC_direction =="Down","State7",
ifelse(EnzavsCtrl_5hmC_direction == "Down" & EnzavsCtrl_5mC_direction == "NotSig", "State8",
ifelse(EnzavsCtrl_5hmC_direction == "NotSig" & EnzavsCtrl_5mC_direction =="NotSig", "State9", "else"))))))))))
table(OxBS.df2$EnzavsCtrl.states)
data <- data.frame( category=c("state1", "state2", "state3","state4", "state5", "state6","state7", "state8", "state9"),
   count=c(1000,1000,1000,1000,1000, 1000,1000,1000,1000))
data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))
mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF","#800000","red","#228b22","#0000cd","#006400")
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
       geom_rect() + scale_fill_manual(values = mycols)+  theme_void()+
       coord_polar(theta="y") +  xlim(c(2, 4) )

=
