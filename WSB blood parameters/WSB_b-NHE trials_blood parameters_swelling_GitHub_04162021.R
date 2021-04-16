# Load required packages
library(readxl)  # read excel spreadsheets
library(pastecs) # for statdesc
library(ggplot2) # figures
library(car) # figures
library(RColorBrewer)
citation()
RStudio.Version()

# Set working directory
setwd("select path/WSB blood parameters")

#show working directory
print(getwd())

# Load spreadsheet

spreadsheet<-"select path/WSB blood parameters/TH_WSB_bNHE trials_blood parameters_02042021.xlsx"

swelling<-read_excel(
    path = spreadsheet, 
    sheet = "swelling")
swelling<-as.data.frame(swelling) # make it a data frame

swelling<-swelling[!(is.na(swelling$Hct)),] # remove empty rows

rm(spreadsheet)
##------------------------------------------------------------------------------------------------
#clean up data for graphs

swelling<-swelling[,-c(6,7,8)]
swelling$drug<-as.factor(swelling$drug)
swelling$Date<-as.factor(swelling$Date)

str(swelling)

##------------------------------------------------------------------------------------------------
#average mass and Hct

mean_mass<-mean(swelling$mass)
SD_mass<-sd(swelling$mass)
n_mass<-as.numeric(length(unique(swelling$mass)))
SE_mass<-SD_mass/sqrt(n_mass)

mean_Hct<-mean(swelling$Hct)
SD_Hct<-sd(swelling$Hct)
n_Hct<-as.numeric(length(swelling$Hct))
SE_Hct<-SD_Hct/sqrt(n_Hct)

rm(mean_mass, SD_mass, SE_mass, n_mass, mean_Hct, SD_Hct, n_Hct, SE_Hct)

##------------------------------------------------------------------------------------------------

##------------------------------------------------------------------------------------------------
# Hct
##------------------------------------------------------------------------------------------------

#means and errorbars
mean<-aggregate(swelling$Hct, by=list(swelling$drug, swelling$time), FUN=mean)
SD<-aggregate(swelling$Hct, by=list(swelling$drug, swelling$time), FUN=sd)
n<-aggregate(swelling$Hct, by=list(swelling$drug, swelling$time), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, time=mean$Group.2, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

(mean$mean[mean$drug=="ISO"&mean$time==60]-mean$mean[mean$drug=="ini"])/
  mean$mean[mean$drug=="ini"]
  

mean_ini<-mean[1,]
mean<-mean[-1,]
mean_ini$drug[1]<-"DMSO"
mean_DMSO<-mean_ini
mean_ini$drug[1]<-"ISO"
mean_ISO<-mean_ini
mean_ini$drug[1]<-"ISO+Am"
mean_ISOAm<-mean_ini
mean_ini$drug[1]<-"ini"
mean<-rbind(mean, mean_DMSO, mean_ISO, mean_ISOAm)
rm(mean_DMSO, mean_ISO, mean_ISOAm)

swelling_ini<-subset(swelling, drug=="ini")
swelling_run<-subset(swelling, !drug=="ini")

#check parametric assumptions

anova<-lm(Hct~drug*time, data= swelling_run)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(swelling_run$Hct, swelling_run$drug, center=median)
leveneTest(swelling_run$Hct, as.factor(swelling_run$time), center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]
pval_time<-summary$`Pr(>F)`[2]
pval_inter<-summary$`Pr(>F)`[3]

rm(anova, anova_res, summary)
#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_time<-formatC(round(pval_time, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_time<-ifelse(pval_time < small_p, paste("time: italic(P) < ", small_p),
                  sprintf("time: italic(P)=='%s'", pval_time))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxtime: italic(P) < ", small_p),
                   sprintf("drugxtime: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc analysis
library(multcompView)

posthoc_Hct_60<-pairwise.t.test(subset(swelling$Hct, swelling$time==60), 
                                       subset(swelling$drug, swelling$time==60), 
                                       p.adjust.method="BH")

posthoc_Hct_60


siglab_Hct<-c("a",  "c", "b")
#                 DMSO, ISO, ISO+Am, 

#axis labels
xlab<-(expression(paste("Time (min)")))
ylab<-(expression(paste("Haematocrit (%)")))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")

Fig_Hct<-ggplot(data=mean, aes(x=time, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_colour_manual(name="Drug", values=c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_x_continuous(name=xlab, breaks = c(0, 10, 30, 60), labels=c("ini", 10, 30, 60), 
                     limits=c(-2, 65))+
  scale_y_continuous(name=ylab, breaks = seq(20, 35, 5), limits=c(20, 37))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=10, position=position_dodge(2))+
  geom_label(x=mean_ini$time, y=mean_ini$mean, colour="white", label="o", size=5,
             show.legend = F, inherit.aes = F)+
  geom_hline(yintercept=mean_ini$mean, colour="black", linetype="dashed", alpha=0.3)+
  
  geom_point(data=mean_ini, aes(x=time, y=mean),colour="black", fill="black",
             position=position_dodge(2), shape=21, size=1, 
             show.legend = F, inherit.aes = F)+
  geom_errorbar(data=mean_ini, aes(x=time, ymin=mean-se, ymax=mean+se),colour="black", 
                width=3.5, position=position_dodge(2), show.legend = FALSE, inherit.aes = F)+
  geom_jitter(data=swelling_ini, aes(x=time, y=Hct), colour="black", fill="black",
              size=1, alpha=0.3, show.legend = FALSE, inherit.aes = F)+
  geom_line()+
  geom_jitter(data=swelling_run, aes(x=time, y=Hct, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=2, dodge.width=2),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(data=subset(mean, time==60),
            y=c(mean$mean[mean$drug=="DMSO"&mean$time==60],
                mean$mean[mean$drug=="ISO"&mean$time==60],
                mean$mean[mean$drug=="ISO+Am"&mean$time==60]+0.6),
            x=64, colour="black",
            label= siglab_Hct,
            size=3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=65 , y=37, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_time, x=65 , y=35.8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=65 , y=34.6, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.15, 0.85),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_Hct

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_Hct_02122021.tiff", unit="cm", width=9, height=7, plot=Fig_Hct)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_Hct_02122021.pdf", unit="cm", width=9, height=7, plot=Fig_Hct, useDingbats=FALSE)

rm(swelling_ini, swelling_run, mean, mean_ini, pval_drug, 
   pval_inter, pval_time, ylab, xlab, col_pal, siglab_Hct, posthoc_Hct_60)


##------------------------------------------------------------------------------------------------
# Hb
##------------------------------------------------------------------------------------------------

#means and errorbars
mean<-aggregate(swelling$Hb, by=list(swelling$drug, swelling$time), FUN=mean, na.rm = TRUE)
SD<-aggregate(swelling$Hb, by=list(swelling$drug, swelling$time), FUN=sd, na.rm = TRUE)
n<-aggregate(swelling$Hb, by=list(swelling$drug, swelling$time), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, time=mean$Group.2, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

mean_ini<-mean[1,]
mean<-mean[-1,]
mean_ini$drug[1]<-"DMSO"
mean_DMSO<-mean_ini
mean_ini$drug[1]<-"ISO"
mean_ISO<-mean_ini
mean_ini$drug[1]<-"ISO+Am"
mean_ISOAm<-mean_ini
mean_ini$drug[1]<-"ini"
mean<-rbind(mean, mean_DMSO, mean_ISO, mean_ISOAm)
rm(mean_DMSO, mean_ISO, mean_ISOAm)

swelling_ini<-subset(swelling, drug=="ini")
swelling_run<-subset(swelling, !drug=="ini")

#check parametric assumptions

anova<-lm(Hb~drug*time, data= swelling_run)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.05, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(swelling_run$Hb, swelling_run$drug, center=median)
leveneTest(swelling_run$Hb, as.factor(swelling_run$time), center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]
pval_time<-summary$`Pr(>F)`[2]
pval_inter<-summary$`Pr(>F)`[3]

rm(anova, anova_res, summary)
#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_time<-formatC(round(pval_time, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_time<-ifelse(pval_time < small_p, paste("time: italic(P) < ", small_p),
                  sprintf("time: italic(P)=='%s'", pval_time))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxtime: italic(P) < ", small_p),
                   sprintf("drugxtime: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc analysis

posthoc_Hb_60<-pairwise.t.test(subset(swelling$Hb, swelling$time==60), 
                                subset(swelling$drug, swelling$time==60), 
                                p.adjust.method="BH")

posthoc_Hb_60


#axis labels
xlab<-(expression(paste("Time (min)")))
ylab<-(expression(paste("Haemoglobin (mM)")))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


Fig_Hb<-ggplot(data=mean, aes(x=time, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_colour_manual(name="Drug", values=c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_x_continuous(name=xlab, breaks = c(0, 10, 30, 60), labels=c("ini", 10, 30, 60), 
                     limits=c(-2, 65))+
  scale_y_continuous(name=ylab, breaks = seq(0.2, 1, 0.2), limits=c(0.2, 1))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=10, position=position_dodge(2))+
  geom_label(x=mean_ini$time, y=mean_ini$mean, colour="white", label="o", size=5,
             show.legend = F, inherit.aes = F)+
  geom_hline(yintercept=mean_ini$mean, colour="black", linetype="dashed", alpha=0.3)+
  
  geom_point(data=mean_ini, aes(x=time, y=mean),colour="black", fill="black",
             position=position_dodge(2), shape=21, size=1, 
             show.legend = F, inherit.aes = F)+
  geom_errorbar(data=mean_ini, aes(x=time, ymin=mean-se, ymax=mean+se),colour="black", 
                width=3.5, position=position_dodge(2), show.legend = FALSE, inherit.aes = F)+
  geom_jitter(data=na.omit(swelling_ini), aes(x=time, y=Hb), colour="black", fill="black",
              size=1, alpha=0.3, show.legend = FALSE, inherit.aes = F)+
  geom_line()+
  geom_jitter(data=na.omit(swelling_run), aes(x=time, y=Hb, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=2, dodge.width=2),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=65 , y=1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_time, x=65 , y=0.95, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=65 , y=0.9, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.15, 0.85),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_Hb

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_Hb_02122021.tiff", unit="cm", width=9, height=7, plot=Fig_Hb)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_Hb_02122021.pdf", unit="cm", width=9, height=7, plot=Fig_Hb, useDingbats=FALSE)

rm(swelling_ini, swelling_run, mean, mean_ini, pval_drug, 
   pval_inter, pval_time, ylab, xlab, col_pal, posthoc_Hb_60)

#overall means and errorbars
mean<-mean(swelling$Hb, na.rm = TRUE)
SD<-sd(swelling$Hb, na.rm = TRUE)
n<-as.numeric(sum(!is.na(swelling$Hb)))
SE<-SD/sqrt(n)
mean<-formatC(round(mean, 3 ), format='f', digits=3)
SE<-formatC(round(SE, 3 ), format='f', digits=3)

Hb_lab<-bquote("[Hb] (mM) ="~.(mean)*"±"*.(SE))

rm(SD, n, SE, mean, Fig_Hb)

##------------------------------------------------------------------------------------------------
# MCHC
##------------------------------------------------------------------------------------------------

#means and errorbars
mean<-aggregate(swelling$MCHC, by=list(swelling$drug, swelling$time), FUN=mean)
SD<-aggregate(swelling$MCHC, by=list(swelling$drug, swelling$time), FUN=sd)
n<-aggregate(swelling$MCHC, by=list(swelling$drug, swelling$time), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, time=mean$Group.2, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

mean_ini<-mean[1,]
mean<-mean[-1,]
mean_ini$drug[1]<-"DMSO"
mean_DMSO<-mean_ini
mean_ini$drug[1]<-"ISO"
mean_ISO<-mean_ini
mean_ini$drug[1]<-"ISO+Am"
mean_ISOAm<-mean_ini
mean_ini$drug[1]<-"ini"
mean<-rbind(mean, mean_DMSO, mean_ISO, mean_ISOAm)
rm(mean_DMSO, mean_ISO, mean_ISOAm)

swelling_ini<-subset(swelling, drug=="ini")
swelling_run<-subset(swelling, !drug=="ini")

#check parametric assumptions

anova<-lm(MCHC~drug*time, data= swelling_run)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(swelling_run$MCHC, swelling_run$drug, center=median)
leveneTest(swelling_run$MCHC, as.factor(swelling_run$time), center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]
pval_time<-summary$`Pr(>F)`[2]
pval_inter<-summary$`Pr(>F)`[3]

rm(anova, anova_res, summary)
#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_time<-formatC(round(pval_time, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_time<-ifelse(pval_time < small_p, paste("time: italic(P) < ", small_p),
                  sprintf("time: italic(P)=='%s'", pval_time))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxtime: italic(P) < ", small_p),
                   sprintf("drugxtime: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc analysis

posthoc_MCHC_60<-pairwise.t.test(subset(swelling$MCHC, swelling$time==60), 
                                subset(swelling$drug, swelling$time==60), 
                                p.adjust.method="BH")

posthoc_MCHC_60


siglab_MCHC<-c("a",  "b", "a")
#                 DMSO, ISO, ISO+Am, 


#axis labels
xlab<-(expression(paste("Time (min)")))
ylab<-(expression(paste("MCHC (mM L RBC"^"-1",")")))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


Fig_MCHC<-ggplot(data=mean, aes(x=time, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_colour_manual(name="Drug", values=c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_x_continuous(name=xlab, breaks = c(0, 10, 30, 60), labels=c("ini", 10, 30, 60), 
                     limits=c(-2, 65))+
  scale_y_continuous(name=ylab, breaks = seq(1, 4, 1), limits=c(1, 4))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=10, position=position_dodge(2))+
  geom_label(x=mean_ini$time, y=mean_ini$mean, colour="white", label="o", size=5,
             show.legend = F, inherit.aes = F)+
  geom_hline(yintercept=mean_ini$mean, colour="black", linetype="dashed", alpha=0.3)+
  
  geom_point(data=mean_ini, aes(x=time, y=mean),colour="black", fill="black",
             position=position_dodge(2), shape=21, size=1, 
             show.legend = F, inherit.aes = F)+
  geom_errorbar(data=mean_ini, aes(x=time, ymin=mean-se, ymax=mean+se),colour="black", 
                width=3.5, position=position_dodge(2), show.legend = FALSE, inherit.aes = F)+
  geom_jitter(data=swelling_ini, aes(x=time, y=MCHC), colour="black", fill="black",
              size=1, alpha=0.3, show.legend = FALSE, inherit.aes = F)+
  geom_line()+
  geom_jitter(data=swelling_run, aes(x=time, y=MCHC, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=2, dodge.width=2),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(data=subset(mean, time==60),
            y=c(mean$mean[mean$drug=="DMSO"&mean$time==60]+0.05,
                mean$mean[mean$drug=="ISO"&mean$time==60],
                mean$mean[mean$drug=="ISO+Am"&mean$time==60]-0.15),
            x=64, colour="black",
            label= siglab_MCHC,
            size=3, show.legend = FALSE)+
  
  geom_text(label=deparse(Hb_lab), x=65 , y=1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  
  geom_text(label=pval_drug, x=65 , y=4, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_time, x=65 , y=3.8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=65 , y=3.6, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.15, 0.85),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_MCHC

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_MCHC_02122021.tiff", unit="cm", width=9, height=7, plot=Fig_MCHC)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_MCHC_02122021.pdf", unit="cm", width=9, height=7, plot=Fig_MCHC, useDingbats=FALSE)

rm(swelling_ini, swelling_run, mean, mean_ini, pval_drug, 
   pval_inter, pval_time, ylab, xlab, col_pal, Hb_lab)

##------------------------------------------------------------------------------------------------
# pHe
##------------------------------------------------------------------------------------------------

#means and errorbars
mean<-aggregate(swelling$pHe, by=list(swelling$drug, swelling$time), FUN=mean, na.rm = TRUE)
SD<-aggregate(swelling$pHe, by=list(swelling$drug, swelling$time), FUN=sd, na.rm = TRUE)
n<-aggregate(swelling$pHe, by=list(swelling$drug, swelling$time), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, time=mean$Group.2, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

mean_ini<-mean[1,]
mean<-mean[-1,]
mean_ini$drug[1]<-"DMSO"
mean_DMSO<-mean_ini
mean_ini$drug[1]<-"ISO"
mean_ISO<-mean_ini
mean_ini$drug[1]<-"ISO+Am"
mean_ISOAm<-mean_ini
mean_ini$drug[1]<-"ini"
mean<-rbind(mean, mean_DMSO, mean_ISO, mean_ISOAm)
rm(mean_DMSO, mean_ISO, mean_ISOAm)

swelling_ini<-subset(swelling, drug=="ini")
swelling_run<-subset(swelling, !drug=="ini")

#check parametric assumptions

anova<-lm(pHe~drug*time, data= swelling_run)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.05, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(swelling_run$pHe, swelling_run$drug, center=median)
leveneTest(swelling_run$pHe, as.factor(swelling_run$time), center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]
pval_time<-summary$`Pr(>F)`[2]
pval_inter<-summary$`Pr(>F)`[3]

rm(anova, anova_res, summary)
#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_time<-formatC(round(pval_time, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_time<-ifelse(pval_time < small_p, paste("time: italic(P) < ", small_p),
                  sprintf("time: italic(P)=='%s'", pval_time))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxtime: italic(P) < ", small_p),
                   sprintf("drugxtime: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc analysis

posthoc_pHe_60<-pairwise.t.test(subset(swelling$pHe, swelling$time==60), 
                                 subset(swelling$drug, swelling$time==60), 
                                 p.adjust.method="BH")

posthoc_pHe_60


siglab_pHe<-c("a",  "b", "a")
#                 DMSO, ISO, ISO+Am, 

#axis labels
xlab<-(expression(paste("Time (min)")))
ylab<-(expression(paste("pH"[e])))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


Fig_pHe<-ggplot(data=mean, aes(x=time, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_colour_manual(name="Drug", values=c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_x_continuous(name=xlab, breaks = c(0, 10, 30, 60), labels=c("ini", 10, 30, 60), 
                     limits=c(-2, 65))+
  scale_y_continuous(name=ylab, breaks = seq(7.4, 8, 0.2), limits=c(7.4, 8.1))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=10, position=position_dodge(2))+
  geom_label(x=mean_ini$time, y=mean_ini$mean, colour="white", label="o", size=5,
             show.legend = F, inherit.aes = F)+
  geom_hline(yintercept=mean_ini$mean, colour="black", linetype="dashed", alpha=0.3)+
  
  geom_point(data=mean_ini, aes(x=time, y=mean),colour="black", fill="black",
             position=position_dodge(2), shape=21, size=1, 
             show.legend = F, inherit.aes = F)+
  geom_errorbar(data=mean_ini, aes(x=time, ymin=mean-se, ymax=mean+se),colour="black", 
                width=3.5, position=position_dodge(2), show.legend = FALSE, inherit.aes = F)+
  geom_jitter(data=na.omit(swelling_ini), aes(x=time, y=pHe), colour="black", fill="black",
              size=1, alpha=0.3, show.legend = FALSE, inherit.aes = F)+
  geom_line()+
  geom_jitter(data=na.omit(swelling_run), aes(x=time, y=pHe, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=2, dodge.width=2),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(data=subset(mean, time==60),
            y=c(mean$mean[mean$drug=="DMSO"&mean$time==60]+0.02,
                mean$mean[mean$drug=="ISO"&mean$time==60],
                mean$mean[mean$drug=="ISO+Am"&mean$time==60]-0.02),
            x=64, colour="black",
            label= siglab_pHe,
            size=3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=65 , y=8.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_time, x=65 , y=8.05, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=65 , y=8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.15, 0.85),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_pHe

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_pHe_02122021.tiff", unit="cm", width=9, height=7, plot=Fig_pHe)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_pHe_02122021.pdf", unit="cm", width=9, height=7, plot=Fig_pHe, useDingbats=FALSE)

rm(swelling_ini, swelling_run, mean, mean_ini, pval_drug, 
   pval_inter, pval_time, ylab, xlab, col_pal)

##------------------------------------------------------------------------------------------------
# pHi
##------------------------------------------------------------------------------------------------

#means and errorbars
mean<-aggregate(swelling$pHi, by=list(swelling$drug, swelling$time), FUN=mean, na.rm = TRUE)
SD<-aggregate(swelling$pHi, by=list(swelling$drug, swelling$time), FUN=sd, na.rm = TRUE)
n<-aggregate(swelling$pHi, by=list(swelling$drug, swelling$time), FUN=function(x)sum(!is.na(x)))
SE<-SD$x/sqrt(n$x)
mean<-data.frame(drug=mean$Group.1, time=mean$Group.2, mean=mean$x, se=SE, n=n$x)
rm(SD, n, SE)

mean_ini<-mean[1,]
mean<-mean[-1,]
mean_ini$drug[1]<-"DMSO"
mean_DMSO<-mean_ini
mean_ini$drug[1]<-"ISO"
mean_ISO<-mean_ini
mean_ini$drug[1]<-"ISO+Am"
mean_ISOAm<-mean_ini
mean_ini$drug[1]<-"ini"
mean<-rbind(mean, mean_DMSO, mean_ISO, mean_ISOAm)
rm(mean_DMSO, mean_ISO, mean_ISOAm)

swelling_ini<-subset(swelling, drug=="ini")
swelling_run<-subset(swelling, !drug=="ini")

#check parametric assumptions

anova<-lm(pHi~drug*time, data= swelling_run)

anova_res<-data.frame(anova$residuals)
stat.desc(anova_res, basic=FALSE, norm=TRUE)
plot(fitted(anova), residuals(anova))
qplot(sample=anova_res$anova.residuals)
ggplot(anova_res, aes(anova.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.05, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(anova_res$anova.residuals, na.rm=TRUE),
                                     sd=sd(anova_res$anova.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(swelling_run$pHi, swelling_run$drug, center=median)
leveneTest(swelling_run$pHi, as.factor(swelling_run$time), center=median)

#run analysis
summary<-Anova(anova, type="II")
summary

#safe p-values
pval_drug<-summary$`Pr(>F)`[1]
pval_time<-summary$`Pr(>F)`[2]
pval_inter<-summary$`Pr(>F)`[3]

rm(anova, anova_res, summary)
#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_time<-formatC(round(pval_time, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_time<-ifelse(pval_time < small_p, paste("time: italic(P) < ", small_p),
                  sprintf("time: italic(P)=='%s'", pval_time))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxtime: italic(P) < ", small_p),
                   sprintf("drugxtime: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc


posthoc_pHi_60<-pairwise.t.test(subset(swelling$pHi, swelling$time==60), 
                                subset(swelling$drug, swelling$time==60), 
                                p.adjust.method="BH")

posthoc_pHi_60

siglab_pHi<-c("a",  "a", "a")
#                 DMSO, ISO, ISO+Am, 

#axis labels
xlab<-(expression(paste("Time (min)")))
ylab<-(expression(paste("pH"[i])))

#display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


Fig_pHi<-ggplot(data=mean, aes(x=time, y=mean, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_colour_manual(name="Drug", values=c(col_pal[3], col_pal[4], col_pal[5]))+
  scale_x_continuous(name=xlab, breaks = c(0, 10, 30, 60), labels=c("ini", 10, 30, 60), 
                     limits=c(-2, 65))+
  scale_y_continuous(name=ylab, breaks = seq(7.3, 7.8, 0.1), limits=c(7.3, 7.8))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=10, position=position_dodge(2))+
  geom_label(x=mean_ini$time, y=mean_ini$mean, colour="white", label="o", size=5,
             show.legend = F, inherit.aes = F)+
  geom_hline(yintercept=mean_ini$mean, colour="black", linetype="dashed", alpha=0.3)+
  
  geom_point(data=mean_ini, aes(x=time, y=mean),colour="black", fill="black",
             position=position_dodge(2), shape=21, size=1, 
             show.legend = F, inherit.aes = F)+
  geom_errorbar(data=mean_ini, aes(x=time, ymin=mean-se, ymax=mean+se),colour="black", 
                width=3.5, position=position_dodge(2), show.legend = FALSE, inherit.aes = F)+
  geom_jitter(data=na.omit(swelling_ini), aes(x=time, y=pHi), colour="black", fill="black",
              size=1, alpha=0.3, show.legend = FALSE, inherit.aes = F)+
  geom_line()+
  geom_jitter(data=na.omit(swelling_run), aes(x=time, y=pHi, fill=drug, colour=drug),
              position=position_jitterdodge(jitter.width=2, dodge.width=2),
              size=1, alpha=0.3, show.legend = FALSE)+
  
  geom_text(data=subset(mean, time==60),
            y=c(mean$mean[mean$drug=="DMSO"&mean$time==60],
                mean$mean[mean$drug=="ISO"&mean$time==60]+0.02,
                mean$mean[mean$drug=="ISO+Am"&mean$time==60]),
            x=64, colour="black",
            label= siglab_pHi,
            size=3, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=65 , y=7.8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_time, x=65 , y=7.76, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=65 , y=7.72, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = c(0.15, 0.85),
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
Fig_pHi

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_pHi_02122021.tiff", unit="cm", width=9, height=7, plot=Fig_pHi)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_RBC swelling_pHi_02122021.pdf", unit="cm", width=9, height=7, plot=Fig_pHi, useDingbats=FALSE)

rm(swelling_ini, swelling_run, mean, mean_ini, pval_drug, 
   pval_inter, pval_time, ylab, xlab, col_pal)

##------------------------------------------------------------------------------------------------
#Combined graph
##------------------------------------------------------------------------------------------------
library(cowplot)

g_comb<-ggdraw() +
  draw_plot(Fig_Hct, x = 0, y = 1/2, width = 1/2, height = 1/2)+
  draw_plot(Fig_MCHC, x = 1/2, y = 1/2, width = 1/2, height = 1/2)+
  draw_plot(Fig_pHe, x = 0, y = 0, width = 1/2, height = 1/2)+
  draw_plot(Fig_pHi, x = 1/2, y = 0, width = 1/2, height = 1/2)+
  
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20,
                  x = c(0, 0.49, 0, 0.49), y = c(1, 1, 0.5, 0.5))
g_comb

ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_swelling_Fig2_swelling combined_02122021.jpeg", unit="cm", width = 16.5, height = 14, plot=g_comb)
ggsave(filename="E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis/WSB_bNHE_swelling_Fig2_swelling combined_02122021.pdf", unit="cm", width = 16.5, height = 14, plot=g_comb, useDingbats=FALSE)

rm(g_comb)