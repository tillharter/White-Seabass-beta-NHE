# Load required packages
library(readxl)  # read excel spreadsheets
library(pastecs) # for statdesc
library(ggplot2) # figures
library(car) # figures
library(RColorBrewer)


# Set working directory
setwd("select path/WSB OECs")

#show working directory
print(getwd())

# Load spreadsheet

spreadsheet<-"select path/WSB OECs/TH_WSB_bNHE trials_OECs_02052021.xlsx"

OEC_dat<-read_excel(
    path = spreadsheet, 
    sheet = "Respiratory")
OEC_dat<-as.data.frame(OEC_dat) # make it a data frame

OEC_dat<-OEC_dat[!(is.na(OEC_dat$date)),] # remove empty rows

rm(spreadsheet)
##------------------------------------------------------------------------------------------------
#clean up data for graph

OEC_dat<-OEC_dat[,-c(9, 11,13, 18, 19, 20, 21,22,23,24,25)]
#OEC_dat<-subset(OEC_dat, drug=="DMSO")
sum(is.na(OEC_dat))

##------------------------------------------------------------------------------------------------
#average mass and Hct

mean_mass<-mean(OEC_dat$mass)
SD_mass<-sd(OEC_dat$mass)
n_mass<-as.numeric(length(unique(OEC_dat$mass)))
SE_mass<-SD_mass/sqrt(n_mass)

mean_Hct<-mean(OEC_dat$Hct)
SD_Hct<-sd(OEC_dat$Hct)
n_Hct<-as.numeric(length(OEC_dat$Hct))
SE_Hct<-SD_Hct/sqrt(n_Hct)

rm(mean_mass, SD_mass, SE_mass, n_mass, mean_Hct, SD_Hct, n_Hct, SE_Hct)


##------------------------------------------------------------------------------------------------
#Bohr effect
OEC_dat$logp50<-log10(OEC_dat$p50)

#regression through all data

bohr_pHe_reg<-lm(logp50~pHe, data= subset(OEC_dat, !N==8))
bohr_pHe_sum<-summary(bohr_pHe_reg)

bohr_pHe_coef<-as.character(formatC(round(bohr_pHe_sum$coefficients[2,1],2), format='f', digits=2))
bohr_pHe_inter<-as.character(formatC(round(bohr_pHe_sum$coefficients[1,1],2), format='f', digits=2))

bohr_pHe_reg_lab<-bquote("log(P"[50]*") ="~.(bohr_pHe_coef)~"x pH"[e]~"+"~.(bohr_pHe_inter))
bohr_pHe_R2<-as.character(formatC(round(bohr_pHe_sum$adj.r.squared, 2), format='f', digits=2))
bohr_pHe_R2<-bquote(italic(R)[adj.]^"2"==~.(bohr_pHe_R2))
  
bohr_pHe_pval<-as.character(formatC(round(anova(bohr_pHe_reg)$'Pr(>F)'[1],3)), format='f', digits=2)
small_p=0.001
bohr_pHe_pval<-ifelse(bohr_pHe_pval < small_p, paste("italic(P) < ", small_p),
                  sprintf("italic(P)=='%s'", bohr_pHe_pval))
rm(small_p,bohr_pHe_coef, bohr_pHe_inter, bohr_pHe_sum)


bohr_pHi_reg<-lm(logp50~pHi, data= subset(OEC_dat, !N==8))
bohr_pHi_sum<-summary(bohr_pHi_reg)

bohr_pHi_coef<-as.character(formatC(round(bohr_pHi_sum$coefficients[2,1],2), format='f', digits=2))
bohr_pHi_inter<-as.character(formatC(round(bohr_pHi_sum$coefficients[1,1],2), format='f', digits=2))

bohr_pHi_reg_lab<-bquote("log(P"[50]*") ="~.(bohr_pHi_coef)~"x pH"[e]~"+"~.(bohr_pHi_inter))
bohr_pHi_R2<-as.character(formatC(round(bohr_pHi_sum$adj.r.squared, 2), format='f', digits=2))
bohr_pHi_R2<-bquote(italic(R)[adj.]^"2"==~.(bohr_pHi_R2))

bohr_pHi_pval<-as.character(formatC(round(anova(bohr_pHi_reg)$'Pr(>F)'[1],3)), format='f', digits=2)
small_p=0.001
bohr_pHi_pval<-ifelse(bohr_pHi_pval < small_p, paste("italic(P) < ", small_p),
                      sprintf("italic(P)=='%s'", bohr_pHi_pval))
rm(small_p,bohr_pHi_coef, bohr_pHi_inter, bohr_pHi_sum)

#average slope and se by individual
bohr<-data.frame()

for (i in 1:length(unique(OEC_dat$N))){
  dat<-subset(OEC_dat, N==i)
  bohr_pHe_reg<-lm(logp50~pHe, data= dat)
  bohr_pHe_sum<-summary(bohr_pHe_reg)
  bohr_pHe_coef<-bohr_pHe_sum$coefficients[2,1]
  
  bohr_pHi_reg<-lm(logp50~pHi, data= dat)
  bohr_pHi_sum<-summary(bohr_pHi_reg)
  bohr_pHi_coef<-bohr_pHi_sum$coefficients[2,1]
  values<-data.frame(bohr_pHe=bohr_pHe_coef, bohr_pHi=bohr_pHi_coef)
  
  bohr<-rbind(bohr, values)
}

bohr<-bohr[-8,]
bohr_mean_pHe<-mean(bohr$bohr_pHe)
bohr_SD_pHe<-sd(bohr$bohr_pHe)
bohr_n_pHe<-as.numeric(length(bohr$bohr_pHe))
bohr_SE_pHe<-bohr_SD_pHe/sqrt(bohr_n_pHe)

bohr_mean_pHi<-mean(bohr$bohr_pHi)
bohr_SD_pHi<-sd(bohr$bohr_pHi)
bohr_n_pHi<-as.numeric(length(bohr$bohr_pHi))
bohr_SE_pHi<-bohr_SD_pHi/sqrt(bohr_n_pHi)

rm(bohr_n_pHe, bohr_n_pHi, bohr_SD_pHe, bohr_SD_pHi, dat, i, 
   values, bohr, bohr_pHe_coef, bohr_pHi_coef)

#supplemental Figure S1 A and B

xlab<-(expression(paste("pH"[e])))
ylab<-(expression(paste("log(P"[50],")")))

FigS1A<-ggplot(data=OEC_dat, aes(x=pHe, y=logp50))+
  scale_x_continuous(name=xlab, breaks = seq(7.2, 8.2, 0.2),limits=c(7.2, 8.2))+
  scale_y_continuous(name=ylab, breaks = seq(1, 2, 0.2),limits=c(1, 2.1))+
  geom_point(colour="grey", show.legend = F)+
  geom_smooth(aes(group=N), colour="grey", method='lm', se=F, formula= y~x, show.legend = F)+
  geom_smooth(method='lm', formula= y~x, show.legend = F)+
  geom_text(label=bohr_pHe_pval, x=7.2 , y=1.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(bohr_pHe_R2), x=7.45 , y=1.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(bohr_pHe_reg_lab), x=7.2 , y=1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
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
FigS1A
#ggsave(filename="WSB_b-NHE_FigS1B_Bohr pHe_02122021.pdf", device = "pdf", unit="cm", width=12, height=8, plot=FigS1A)


xlab<-(expression(paste("pH"[i])))
ylab<-(expression(paste("log(P"[50],")")))

FigS1B<-ggplot(data=subset(OEC_dat, !N==8), aes(x=pHi, y=logp50))+
  scale_x_continuous(name=xlab, breaks = seq(7, 7.8, 0.2),limits=c(7, 7.8))+
  scale_y_continuous(name=ylab, breaks = seq(1, 2, 0.2),limits=c(1, 2.1))+
  geom_point(colour="grey", show.legend = F)+
  geom_smooth(aes(group=N), colour="grey", method='lm', se=F, formula= y~x, show.legend = F)+
  geom_smooth(method='lm', formula= y~x, show.legend = F)+
  geom_text(label=bohr_pHi_pval, x=7.0 , y=1.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(bohr_pHi_R2), x=7.2 , y=1.1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(bohr_pHi_reg_lab), x=7.0 , y=1, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
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
FigS1B
#ggsave(filename="WSB_b-NHE_FigS1A_Bohr pHi_02122021.tiff", unit="cm", width=12, height=8, plot=FigS1B)

#Root effect from other dataset
R_coef<--52.4
R_se<-1.8

rm(bohr_pHe_reg, bohr_pHi_reg, bohr_pHe_sum, bohr_pHi_sum, xlab, ylab)

#labels for main Figure

bohr_mean_pHe<-as.character(formatC(round(bohr_mean_pHe, 2 ), format='f', digits=2))
bohr_SE_pHe<-as.character(formatC(round(bohr_SE_pHe, 2 ), format='f', digits=2))
bohr_mean_pHi<-as.character(formatC(round(bohr_mean_pHi, 2 ), format='f', digits=2))
bohr_SE_pHi<-as.character(formatC(round(bohr_SE_pHi, 2 ), format='f', digits=2))
bohr_pHe_lab<-bquote("B (pH"[e]~") ="~.(bohr_mean_pHe)*"±"*.(bohr_SE_pHe))
bohr_pHi_lab<-bquote("B (pH"[i]~") ="~.(bohr_mean_pHi)*"±"*.(bohr_SE_pHi))

R_coef<-as.character(formatC(round(R_coef, 1 ), format='f', digits=1))
R_se<-as.character(formatC(round(R_se, 1 ), format='f', digits=1))
R_lab<-bquote("R"~"="~.(R_coef)*"±"*.(R_se)*"%")

rm(bohr_mean_pHe, bohr_mean_pHi, bohr_SE_pHe, bohr_SE_pHi, R_coef, R_se,
   bohr_pHe_pval, bohr_pHi_pval, bohr_pHe_R2, bohr_pHi_R2, bohr_pHe_reg_lab, bohr_pHi_reg_lab)

##------------------------------------------------------------------------------------------------
#non-bicarb buffer capacity

nb_buffer<-data.frame()

for (i in 1:length(unique(OEC_dat$N))){
  dat<-subset(OEC_dat, N==i)
  HCO3_pHe_reg<-lm(HCO3~pHe, data= dat)
  HCO3_pHe_sum<-summary(HCO3_pHe_reg)
  HCO3_pHe_coef<-HCO3_pHe_sum$coefficients[2,1]
  
  nb_buffer<-rbind(nb_buffer, HCO3_pHe_coef)
}
#regression through all data
HCO3_reg<-lm(HCO3~pHe, data= subset(OEC_dat, !N==7))

HCO3_sum<-summary(HCO3_reg)

HCO3_coef<-as.character(formatC(round(HCO3_sum$coefficients[2,1],2), format='f', digits=2))
HCO3_inter<-as.character(formatC(round(HCO3_sum$coefficients[1,1],2), format='f', digits=2))

HCO3_reg_lab<-bquote("[HCO"[3]^"-"~"] (mM) ="~.(HCO3_coef)~"x pH"[e]~"+"~.(HCO3_inter))
HCO3_R2<-as.character(formatC(round(HCO3_sum$adj.r.squared, 2), format='f', digits=2))
HCO3_R2<-bquote(italic(R)[adj.]^"2"==~.(HCO3_R2))

HCO3_pval<-as.character(formatC(round(anova(HCO3_reg)$'Pr(>F)'[1],3)), format='f', digits=2)
small_p=0.001
HCO3_pval<-ifelse(HCO3_pval < small_p, paste("italic(P) < ", small_p),
                      sprintf("italic(P)=='%s'", HCO3_pval))
rm(small_p,HCO3_coef, HCO3_inter, HCO3_sum)


#average slope and se by individual
colnames(nb_buffer)<-"slopes"
nb_buffer<-subset(nb_buffer, slopes<0)
nb_buffer_mean<-mean(nb_buffer$slopes)
nb_buffer_SD<-sd(nb_buffer$slopes)
nb_buffer_n<-as.numeric(length(nb_buffer$slopes))
nb_buffer_SE<-nb_buffer_SD/sqrt(nb_buffer_n)

#Hct correction after Wood and McDonald 1982

nb_buffer_25<--28.35*0.25-2.59
nb_buffer_5<--28.35*0.05-2.59
dif<-nb_buffer_25-nb_buffer_5

nb_buffer_wb<-nb_buffer_mean+dif

rm(nb_buffer_25, nb_buffer_5, dif, dat, i)


#supplemental Figure
xlab<-(expression(paste("pH"[e])))
ylab<-(expression(paste("plasma [HCO"[3]^"-", "] (mM)")))

FigS1C<-ggplot(data=subset(OEC_dat, !N==7), aes(x=pHe, y=HCO3))+
  scale_x_continuous(name=xlab, breaks = seq(7.2, 8.2, 0.2),limits=c(7.2, 8.2))+
  scale_y_continuous(name=ylab, breaks = seq(5, 12, 2),limits=c(5, 12))+
  geom_point(colour="grey", show.legend = F)+
  geom_smooth(aes(group=N), colour="grey", method='lm', se=F, formula= y~x, show.legend = F)+
  geom_smooth(method='lm', formula= y~x, show.legend = F)+
  geom_text(label=HCO3_pval, x=7.93 , y=12, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(HCO3_R2), x=8.2 , y=12, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(HCO3_reg_lab), x=8.2 , y=11.4, colour="black", size=3,
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
FigS1C
#ggsave(filename="WSB_b-NHE_FigS2_nonbirab buffer_02122021.tiff", unit="cm", width=12, height=8, plot=FigS1C)

#labels for main Figure
nb_buffer_mean<-as.character(formatC(round(nb_buffer_mean, 2 ), format='f', digits=2))
nb_buffer_SE<-as.character(formatC(round(nb_buffer_SE, 2 ), format='f', digits=2))
nb_buffer_lab<-bquote("Non-bicarb buffer (Hct 5%) ="~.(nb_buffer_mean)*"±"*.(nb_buffer_SE))

rm(nb_buffer, nb_buffer_SD, HCO3_pHe_coef, HCO3_reg, HCO3_pHe_sum, ylab, xlab, HCO3_pHe_reg,
   nb_buffer_mean, nb_buffer_SE, nb_buffer_n, nb_buffer_wb, HCO3_R2, HCO3_pval, HCO3_reg_lab)

##------------------------------------------------------------------------------------------------
#pHi vs pHe

#regression through all data
pHi_reg<-lm(pHi~pHe, data= subset(OEC_dat, !N==8))

pHi_sum<-summary(pHi_reg)

pHi_coef<-as.character(formatC(round(pHi_sum$coefficients[2,1],2), format='f', digits=2))
pHi_inter<-as.character(formatC(round(pHi_sum$coefficients[1,1],2), format='f', digits=2))

pHi_reg_lab<-bquote("pH"[e]~"="~.(pHi_coef)~"x pH"[e]~"+"~.(pHi_inter))
pHi_R2<-as.character(formatC(round(pHi_sum$adj.r.squared, 2), format='f', digits=2))
pHi_R2<-bquote(italic(R)[adj.]^"2"==~.(pHi_R2))

pHi_pval<-as.character(formatC(round(anova(pHi_reg)$'Pr(>F)'[1],3)), format='f', digits=2)
small_p=0.001
pHi_pval<-ifelse(pHi_pval < small_p, paste("italic(P) < ", small_p),
                  sprintf("italic(P)=='%s'", pHi_pval))
rm(small_p, pHi_sum)

#average slope and se by individual

pHi_vs_pHe<-data.frame()

for (i in 1:length(unique(OEC_dat$N))){
  dat<-subset(OEC_dat, N==i)
  pHi_pHe_reg<-lm(pHi~pHe, data= dat)
  pHi_pHe_sum<-summary(pHi_pHe_reg)
  pHi_pHe_coef<-pHi_pHe_sum$coefficients[2,1]
  
  pHi_vs_pHe<-rbind(pHi_vs_pHe, pHi_pHe_coef)
}

rm(pHi_pHe_coef, pHi_pHe_reg, pHi_pHe_sum, i, dat)

colnames(pHi_vs_pHe)<-"slopes"
pHi_vs_pHe<-subset(pHi_vs_pHe, !slopes<0.2)
pHi_vs_pHe_mean<-mean(pHi_vs_pHe$slopes)
pHi_vs_pHe_SD<-sd(pHi_vs_pHe$slopes)
pHi_vs_pHe_n<-as.numeric(length(pHi_vs_pHe$slopes))
pHi_vs_pHe_SE<-pHi_vs_pHe_SD/sqrt(pHi_vs_pHe_n)

#supplemental Figure
xlab<-(expression(paste("pH"[e])))
ylab<-(expression(paste("pH"[i])))

FigS1D<-ggplot(data=subset(OEC_dat, !N==8), aes(x=pHe, y=pHi))+
  scale_x_continuous(name=xlab, breaks = seq(7.2, 8.2, 0.2),limits=c(7.2, 8.2))+
  scale_y_continuous(name=ylab, breaks = seq(7.0, 7.8, 0.2),limits=c(7.0, 7.8))+
  geom_point(colour="grey", show.legend = F)+
  geom_smooth(aes(group=N), colour="grey", method='lm', se=F, formula= y~x, show.legend = F)+
  geom_smooth(method='lm', formula= y~x, show.legend = F)+
  geom_text(label=pHi_pval, x=7.2 , y=7.8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(pHi_R2), x=7.45 , y=7.8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(pHi_reg_lab), x=7.2 , y=7.73, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
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
FigS1D

#ggsave(filename="WSB_b-NHE_FigS3_phi_phe_02112021.tiff", unit="cm", width=12, height=8, plot=FigS1D)

#labels for main Figure
pHi_vs_pHe_mean<-as.character(formatC(round(pHi_vs_pHe_mean, 2 ), format='f', digits=2))
pHi_vs_pHe_SE<-as.character(formatC(round(pHi_vs_pHe_SE, 2 ), format='f', digits=2))
pHi_vs_pHe_lab<-bquote("pH"[i]~"="~.(pHi_vs_pHe_mean)*"±"*.(pHi_vs_pHe_SE)*" x pH"[e])

rm(pHi_vs_pHe, pHi_vs_pHe_mean, pHi_vs_pHe_SE, pHi_vs_pHe_n, pHi_vs_pHe_SD, xlab, ylab,
   pHi_inter, pHi_coef, pHi_R2, pHi_pval, pHi_reg_lab, pHi_reg)

#Combined graph
##------------------------------------------------------------------------------------------------
library(cowplot)

g_comb<-ggdraw() +
  draw_plot(FigS1A, x = 0, y = 1/2, width = 1/2, height = 1/2)+
  draw_plot(FigS1B, x = 1/2, y = 1/2, width = 1/2, height = 1/2)+
  draw_plot(FigS1D, x = 0, y = 0, width = 1/2, height = 1/2)+
  draw_plot(FigS1C, x = 1/2, y = 0, width = 1/2, height = 1/2)+
  
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20,
                  x = c(-0.01, 1/2, -0.01, 1/2), y = c(1, 1, 1/2, 1/2))
g_comb

ggsave(filename="WSB_bNHE_Figure S1_combined.jpeg", unit="cm", width = 18, height = 14, plot=g_comb)
ggsave(filename="WSB_bNHE_Figure S1_combined.pdf", unit="cm", width = 18, height = 14, plot=g_comb, useDingbats=FALSE)

rm(g_comb, FigS1A, FigS1B, FigS1C, FigS1D)

##------------------------------------------------------------------------------------------------
#OECs

mean_p50<-aggregate(OEC_dat$p50, by=list(OEC_dat$PCO2_BOBS), FUN=mean, na.rm=TRUE)
se_p50<-aggregate(OEC_dat$p50, by=list(OEC_dat$PCO2_BOBS), FUN=sd, na.rm=TRUE)
n_p50<-aggregate(OEC_dat$p50, by=list(OEC_dat$PCO2_BOBS), FUN=length)
se_p50<-se_p50$x/(sqrt(n_p50$x))

mean_nH<-aggregate(OEC_dat$nH, by=list(OEC_dat$PCO2_BOBS), FUN=mean, na.rm=TRUE)
se_nH<-aggregate(OEC_dat$nH, by=list(OEC_dat$PCO2_BOBS), FUN=sd, na.rm=TRUE)
n_nH<-aggregate(OEC_dat$nH, by=list(OEC_dat$PCO2_BOBS), FUN=length)
se_nH<-se_nH$x/(sqrt(n_nH$x))

mean<-data.frame(CO2=mean_p50$Group.1, p50=mean_p50$x, se_p50, n_p50=n_p50$x,
                 nH=mean_nH$x, se_nH, n_nH=n_nH$x)
rm(mean_nH, mean_p50, n_nH, n_p50, se_nH, se_p50)

#Stats
reg_p50<-lm(p50~PCO2_BOBS, data=OEC_dat)

#check parametric assumptions
reg_p50_res<-data.frame(reg_p50$residuals)
stat.desc(reg_p50_res, basic=FALSE, norm=TRUE)
plot(fitted(reg_p50), residuals(reg_p50))
qplot(sample=reg_p50_res$reg_p50.residuals)
ggplot(reg_p50_res, aes(reg_p50.residuals))+
  geom_histogram(aes(y=..density..), binwidth=5, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(reg_p50_res$reg_p50.residuals, na.rm=TRUE),
                                     sd=sd(reg_p50_res$reg_p50.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(OEC_dat$p50, OEC_dat$PCO2_BOBS, center=median)

anova_reg_p50<-Anova(reg_p50, type="II")
anova_reg_p50


reg_nH<-lm(nH~PCO2_BOBS, data=OEC_dat)

#check parametric assumptions
reg_nH_res<-data.frame(reg_nH$residuals)
stat.desc(reg_nH_res, basic=FALSE, norm=TRUE)
plot(fitted(reg_nH), residuals(reg_nH))
qplot(sample=reg_nH_res$reg_nH.residuals)
ggplot(reg_nH_res, aes(reg_nH.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(reg_nH_res$reg_nH.residuals, na.rm=TRUE),
                                     sd=sd(reg_nH_res$reg_nH.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(OEC_dat$nH, OEC_dat$PCO2_BOBS, center=median)

anova_reg_nH<-Anova(reg_nH, type="II")
anova_reg_nH
#safe p-values

pval_p50<-anova_reg_p50$`Pr(>F)`[1]
pval_nH<-anova_reg_nH$`Pr(>F)`[1]


rm(anova_reg_p50, reg_p50, reg_p50_res,
   anova_reg_nH, reg_nH, reg_nH_res)

#significance labels

small_p=0.001

pval_p50<-formatC(round(pval_p50, 3 ), format='f', digits=3)
pval_nH<-formatC(round(pval_nH, 3 ), format='f', digits=3)

pval_p50<-ifelse(pval_p50 < small_p, pval_p50<-small_p,
                 pval_p50<-pval_p50)
pval_nH<-ifelse(pval_nH < small_p, pval_p50<-small_p,
                 pval_nH<-pval_nH)
rm(small_p)

pval_p50<-bquote("P"[50]*": "*italic(P)==.(pval_p50))
pval_nH<-bquote("n"[H]*": "*italic(P)==.(pval_nH))

#simulate curves

oec_0.3<-function(x){(x^mean$nH[1])/(mean$p50[1]^mean$nH[1]+x^mean$nH[1])*100}
oec_0.6<-function(x){(x^mean$nH[2])/(mean$p50[2]^mean$nH[2]+x^mean$nH[2])*100}
oec_1<-function(x){(x^mean$nH[3])/(mean$p50[3]^mean$nH[3]+x^mean$nH[3])*100}
oec_1.5<-function(x){(x^mean$nH[4])/(mean$p50[4]^mean$nH[4]+x^mean$nH[4])*100}
oec_2.5<-function(x){(x^mean$nH[5])/(mean$p50[5]^mean$nH[5]+x^mean$nH[5])*100}

  
##------------------------------------------------------------------------------------------------
#figure labels
xlab<-(expression(paste("PO"[2], " (mmHg)")))
ylab<-(expression(paste("Hb-O"[2], " sat. (%)")))
display.brewer.pal(8,"Reds")
col_pal<-brewer.pal(6, "Reds")
col_pal<-c("#FCBBA1", "#FC9272", "#FB6A4A","#DE2D26", "#A50F15")


#Figure 1
Fig1<-ggplot(data.frame(x = c(0, 160)), aes(x))+
  scale_x_continuous(name=xlab, breaks = seq(0, 160, 20),limits=c(0, 203))+
  scale_y_continuous(name=ylab, breaks = seq(0, 100, 20),limits=c(0, 100))+
  stat_function(fun=oec_0.3, xlim = c(0, 160), colour=col_pal[5])+
  stat_function(fun=oec_0.6, xlim = c(0, 160), colour=col_pal[4])+
  stat_function(fun=oec_1, xlim = c(0, 160), colour=col_pal[3])+
  stat_function(fun=oec_1.5, xlim = c(0, 160), colour=col_pal[2])+
  stat_function(fun=oec_2.5, xlim = c(0, 160), colour=col_pal[1])+
  geom_text(label=deparse(bohr_pHe_lab), x=210 , y=32, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(bohr_pHi_lab), x=210 , y=24, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(R_lab), x=210 , y=16, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(pHi_vs_pHe_lab), x=210 , y=8, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=deparse(nb_buffer_lab), x=210 , y=0, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  
  geom_text(label=deparse(bquote("PCO"[2])), x=165 , y=101, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$CO2[1],1), x=165 , y=95, colour=col_pal[5], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$CO2[2],1), x=165 , y=85, colour=col_pal[4], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(sprintf("%.1f", mean$CO2[3])), x=165 , y=76, colour=col_pal[3], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$CO2[4],2), x=165 , y=67, colour=col_pal[2], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$CO2[5],2), x=165 , y=60, colour=col_pal[1], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  
  geom_text(label=deparse(bquote("P"[50])), x=183 , y=101, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$p50[1],3), x=183 , y=95, colour=col_pal[5], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$p50[2],3), x=183 , y=85, colour=col_pal[4], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(sprintf("%.1f", mean$p50[3],3)), x=183 , y=76, colour=col_pal[3], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$p50[4],3), x=183 , y=67, colour=col_pal[2], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$p50[5],3), x=183 , y=60, colour=col_pal[1], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  
  geom_text(label=deparse(bquote("n"[H])), x=200 , y=101, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$nH[1],3), x=200 , y=95, colour=col_pal[5], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$nH[2],3), x=200 , y=85, colour=col_pal[4], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$nH[3],2), x=200 , y=76, colour=col_pal[3], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$nH[4],2), x=200 , y=67, colour=col_pal[2], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=signif(mean$nH[5],2), x=200 , y=60, colour=col_pal[1], size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(pval_p50), x=0 , y=100, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
  geom_text(label=deparse(pval_nH), x=0 , y=93, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=0)+
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
Fig1

ggsave(filename="WSB_b-NHE_Fig1 OECs_02122021.jpeg", unit="cm", width=12, height=8, plot=Fig1)
ggsave(filename="WSB_b-NHE_Fig1 OECs_02122021.pdf", unit="cm", width=12, height=8, plot=Fig1, useDingbats=FALSE)

rm(Fig1, mean, xlab, ylab,
   R_lab, nb_buffer_lab, pHi_vs_pHe_lab, oec_0.3, oec_0.6, oec_1, oec_1.5, oec_2.5,
   bohr_pHe_lab, bohr_pHi_lab, col_pal, pval_nH, pval_p50)
