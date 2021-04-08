# Load required packages
library(readxl)  # read excel spreadsheets
library(pastecs) # for statdesc
library(ggplot2) # figures
library(car) # figures
library(RColorBrewer)

# Set working directory
setwd("E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials")

#show working directory
print(getwd())

# Load spreadsheet

spreadsheet<-"E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/TH_WSB_bNHE trials_blood parameters_02042021.xlsx"

BOBS<-read_excel(
    path = spreadsheet, 
    sheet = "BOBS")
BOBS<-as.data.frame(BOBS) # make it a data frame

BOBS<-BOBS[!(is.na(BOBS$Hct)),] # remove empty rows

rm(spreadsheet)
##------------------------------------------------------------------------------------------------
#clean up data for graphs

BOBS<-BOBS[,-c(2,3,4,5,6,7,8,9, 16,17,18,19)]
BOBS$Date<-as.character(BOBS$Date)
BOBS<-subset(BOBS, !Date=="2021-01-27") #remove in vivo stress data
BOBS<-subset(BOBS, !Date=="2021-02-10") #remove fixation trial data

BOBS$drug<-as.factor(BOBS$drug)

str(BOBS)
sum(is.na(BOBS))

#relevel drug
BOBS$drug<-factor(BOBS$drug, levels=c("Ctrl", "ISO", "ISO+Am", "ISO+CA"),
                          labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))

##------------------------------------------------------------------------------------------------
# Hct
##------------------------------------------------------------------------------------------------

#overall means and errorbars
mean<-mean(BOBS$Hct, na.rm = TRUE)
SD<-sd(BOBS$Hct, na.rm = TRUE)
n<-as.numeric(sum(!is.na(BOBS$Hct)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=2)
SE<-formatC(round(SE, 3 ), format='f', digits=2)

Hct_lab<-bquote("Average Hct ="~.(mean)*"±"*.(SE)*"%")

rm(SD, n, SE, mean)

#drug means and errorbars

mean<-aggregate(BOBS$Hct, by=list(BOBS$drug), FUN=mean)
SD<-aggregate(BOBS$Hct, by=list(BOBS$drug), FUN=sd)
n<-aggregate(BOBS$Hct, by=list(BOBS$drug), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

#check parametric assumptions

anova<-lm(Hct~drug, data= BOBS)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.5, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(BOBS$Hct, BOBS$drug, center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]

rm(anova, anova_res, summary)

#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
rm(small_p)

#axis labels

ylab<-(expression(paste("Haematocrit (%)")))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


Fig_Hct<-ggplot(data=mean, aes(x=drug, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete()+
  scale_y_continuous(name=ylab, breaks = seq(0, 8, 2), limits=c(0, 8))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=0.2)+
  geom_jitter(data=BOBS, aes(x=drug, y=Hct, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=0.5),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(Hct_lab), x=0.5 , y=0.2, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.88, 0.2),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_Hct

#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_Hct_02162021.jpeg", unit="cm", width=9, height=7, plot=Fig_Hct)
#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_Hct_02162021.pdf", unit="cm", width=9, height=7, plot=Fig_Hct, useDingbats=FALSE)

rm(mean, pval_drug, ylab, Hct_lab)

##------------------------------------------------------------------------------------------------
# Hb
##------------------------------------------------------------------------------------------------
#overall means and errorbars
mean<-mean(BOBS$Hb, na.rm = TRUE)
SD<-sd(BOBS$Hb, na.rm = TRUE)
n<-as.numeric(sum(!is.na(BOBS$Hb)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=3)
SE<-formatC(round(SE, 3 ), format='f', digits=3)

Hb_lab<-bquote("Average [Hb] ="~.(mean)*"±"*.(SE)~"(mM)")

rm(SD, n, SE, mean)

#drug means and errorbars
mean<-aggregate(BOBS$Hb, by=list(BOBS$drug), FUN=mean)
SD<-aggregate(BOBS$Hb, by=list(BOBS$drug), FUN=sd)
n<-aggregate(BOBS$Hb, by=list(BOBS$drug), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

#check parametric assumptions

anova<-lm(Hb~drug, data= BOBS)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.05, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(BOBS$Hb, BOBS$drug, center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]

rm(anova, anova_res, summary)

#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
rm(small_p)

#axis labels

ylab<-(expression(paste("Haemoglobin (mM)")))

Fig_Hb<-ggplot(data=mean, aes(x=drug, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete()+
  scale_y_continuous(name=ylab, breaks = seq(0, 0.3, 0.1), limits=c(0, 0.3))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=0.2)+
  geom_jitter(data=BOBS, aes(x=drug, y=Hb, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=0.5),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=0.3, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(Hb_lab), x=0.5 , y=0.01, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.88, 0.2),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_Hb

#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_Hb_02162021.jpeg", unit="cm", width=9, height=7, plot=Fig_Hb)
#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_Hb_02162021.pdf", unit="cm", width=9, height=7, plot=Fig_Hb, useDingbats=FALSE)

rm(mean, pval_drug, ylab, Hb_lab)

##------------------------------------------------------------------------------------------------
# MCHC
##------------------------------------------------------------------------------------------------
#overall means and errorbars
mean<-mean(BOBS$MCHC, na.rm = TRUE)
SD<-sd(BOBS$MCHC, na.rm = TRUE)
n<-as.numeric(sum(!is.na(BOBS$MCHC)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=2)
SE<-formatC(round(SE, 3 ), format='f', digits=2)

MCHC_lab<-bquote("Average MCHC (mM L"^-1~"RBC) ="~.(mean)*"±"*.(SE))

rm(SD, n, SE, mean)

#drug means and errorbars
mean<-aggregate(BOBS$MCHC, by=list(BOBS$drug), FUN=mean)
SD<-aggregate(BOBS$MCHC, by=list(BOBS$drug), FUN=sd)
n<-aggregate(BOBS$MCHC, by=list(BOBS$drug), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

#check parametric assumptions

anova<-lm(MCHC~drug, data= BOBS)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(BOBS$MCHC, BOBS$drug, center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]

rm(anova, anova_res, summary)

#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
rm(small_p)

#axis labels

ylab<-(expression(paste("MCHC (mM L RBC"^"-1",")")))

Fig_MCHC<-ggplot(data=mean, aes(x=drug, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete()+
  scale_y_continuous(name=ylab, breaks = seq(0, 5, 1), limits=c(0, 5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=0.2)+
  geom_jitter(data=BOBS, aes(x=drug, y=MCHC, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=0.5),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=5, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(MCHC_lab), x=0.5 , y=0.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.88, 0.3),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_MCHC

#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_MCHC_02162021.jpeg", unit="cm", width=9, height=7, plot=Fig_MCHC)
#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_MCHC_02162021.pdf", unit="cm", width=9, height=7, plot=Fig_MCHC, useDingbats=FALSE)

rm(mean, pval_drug, ylab, MCHC_lab)

##------------------------------------------------------------------------------------------------
#pHe
##------------------------------------------------------------------------------------------------
#overall means and errorbars
mean<-mean(BOBS$pHe, na.rm = TRUE)
SD<-sd(BOBS$pHe, na.rm = TRUE)
n<-as.numeric(sum(!is.na(BOBS$pHe)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=3)
SE<-formatC(round(SE, 3 ), format='f', digits=3)

pHe_lab<-bquote("Average pH"[e]~"="~.(mean)*"±"*.(SE))

rm(SD, n, SE, mean)

#drug means and errorbars
mean<-aggregate(BOBS$pHe, by=list(BOBS$drug), FUN=mean)
SD<-aggregate(BOBS$pHe, by=list(BOBS$drug), FUN=sd)
n<-aggregate(BOBS$pHe, by=list(BOBS$drug), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

#check parametric assumptions

anova<-lm(pHe~drug, data= BOBS)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(BOBS$pHe, BOBS$drug, center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]

rm(anova, anova_res, summary)

#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
rm(small_p)

#axis labels

ylab<-(expression(paste("pH"[e])))

Fig_pHe<-ggplot(data=mean, aes(x=drug, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete()+
  scale_y_continuous(name=ylab, breaks = seq(7.2, 8.2, 0.2), limits=c(7.2, 8.2))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=0.2)+
  geom_jitter(data=BOBS, aes(x=drug, y=pHe, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=0.5),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=8.2, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(pHe_lab), x=0.5 , y=7.25, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.88, 0.2),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_pHe

#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_pHe_02162021.jpeg", unit="cm", width=9, height=7, plot=Fig_pHe)
#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_pHe_02162021.pdf", unit="cm", width=9, height=7, plot=Fig_pHe, useDingbats=FALSE)

rm(mean, pval_drug, ylab, pHe_lab)


##------------------------------------------------------------------------------------------------
# pHi
##------------------------------------------------------------------------------------------------
#overall means and errorbars
mean<-mean(BOBS$pHi, na.rm = TRUE)
SD<-sd(BOBS$pHi, na.rm = TRUE)
n<-as.numeric(sum(!is.na(BOBS$pHi)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=3)
SE<-formatC(round(SE, 3 ), format='f', digits=3)

pHi_lab<-bquote("Average pH"[i]~"="~.(mean)*"±"*.(SE))

rm(SD, n, SE, mean)

#drug means and errorbars
mean<-aggregate(BOBS$pHi, by=list(BOBS$drug), FUN=mean)
SD<-aggregate(BOBS$pHi, by=list(BOBS$drug), FUN=sd)
n<-aggregate(BOBS$pHi, by=list(BOBS$drug), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

#check parametric assumptions

anova<-lm(pHi~drug, data= BOBS)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(BOBS$pHi, BOBS$drug, center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]

rm(anova, anova_res, summary)

#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
rm(small_p)

#axis labels

ylab<-(expression(paste("pH"[i])))

Fig_pHi<-ggplot(data=mean, aes(x=drug, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete()+
  scale_y_continuous(name=ylab, breaks = seq(7, 8, 0.2), limits=c(7, 8))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=0.2)+
  geom_jitter(data=BOBS, aes(x=drug, y=pHi, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=0.5),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(pHi_lab), x=0.5 , y=7.05, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.88, 0.2),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_pHi

#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_pHi_02162021.jpeg", unit="cm", width=9, height=7, plot=Fig_pHi)
#ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC_BOBS_pHi_02162021.pdf", unit="cm", width=9, height=7, plot=Fig_pHi, useDingbats=FALSE)

rm(mean, pval_drug, ylab, pHi_lab)


##------------------------------------------------------------------------------------------------
#Combined graph
##------------------------------------------------------------------------------------------------
library(cowplot)

RBC_comb<-ggdraw() +
  draw_plot(Fig_Hct, x = 0.031, y = 2/3, width = (1-0.031), height = 1/3)+
  draw_plot(Fig_Hb, x = 0, y = 1/3, width = 1, height = 1/3)+
  draw_plot(Fig_MCHC, x = 0.017, y = 0, width = (1-0.017), height = 1/3)+
  
  draw_plot_label(label = c("A", "B", "C"), size = 20,
                  x = c(-0.01, -0.01, -0.01), y = c(1, 2/3, 1/3))
RBC_comb

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_RBC params_combined_02162021.jpeg", unit="cm", width = 9, height = 21, plot=RBC_comb)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_RBC params_combined_02162021.pdf", unit="cm", width = 9, height = 21, plot=RBC_comb, useDingbats=FALSE)

rm(RBC_comb)

pH_comb<-ggdraw() +
draw_plot(Fig_pHe, x = 0, y = 1/2, width = 1, height = 1/2)+
  draw_plot(Fig_pHi, x = 0, y = 0, width = 1, height = 1/2)+
  
  draw_plot_label(label = c("A", "B"), size = 20,
                  x = c(-0.01, -0.01), y = c(1, 1/2))
pH_comb

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_pH_combined_02162021.jpeg", unit="cm", width = 9, height = 14, plot=pH_comb)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_pH_combined_02162021.pdf", unit="cm", width = 9, height = 14, plot=pH_comb, useDingbats=FALSE)

rm(pH_comb)


Comb<-ggdraw() +
  draw_plot(Fig_Hct, x = 0.018, y = 2/3, width = (0.5-0.018), height = 1/3)+
  draw_plot(Fig_MCHC, x = 0.01, y = 1/3, width = (0.5-0.01), height = 1/3)+
  draw_plot(Fig_pHe, x = 0, y = 0, width = 0.5, height = 1/3)+
  draw_plot(Fig_Hb, x = 0.5, y = 2/3, width = 0.5, height = 1/3)+
  draw_plot(Fig_pHi, x = 0.5, y = 0, width = 0.5, height = 1/3)+
  
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 20,
                  x = c(-0.01, 0.49, -0.01, -0.01, 0.49), y = c(1, 1, 2/3, 1/3, 1/3))
Comb

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_blood params_combined_02162021.jpeg", unit="cm", width = 18, height = 21, plot=Comb)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_BOBS_FigS3_blood params_combined_02162021.pdf", unit="cm", width = 18, height = 21, plot=Comb, useDingbats=FALSE)

rm(Comb)


rm(Fig_Hb, Fig_Hct, Fig_MCHC, Fig_pHe, Fig_pHi, col_pal)
