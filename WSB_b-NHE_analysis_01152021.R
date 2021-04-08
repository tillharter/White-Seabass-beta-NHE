library(ggplot2)
library(dplyr)
library(stringr)
library(pastecs) # for statdesc
library(car) # for leveneTest
require(minpack.lm) #for nlsLM
library(RColorBrewer) #for colour palett


#--------------------------------------------------------------------------------------------------
# Load BOBS text files
#--------------------------------------------------------------------------------------------------

# Set working directory
setwd("E:/Uni/Postdoc/Data/RBC_b-NHE_project/WSB_b-NHE trials/Analysis")

#show working directory
print(getwd())

#list all files
myfiles <- list.files(path="E:/Uni/Postdoc/data/RBC_b-NHE_project/WSB_b-NHE trials/",
                         recursive=T, pattern=".txt", include.dirs = T, full.names=T)
print(myfiles)

#exclude the gas protocol text files
myfiles<-subset(myfiles, str_detect(myfiles, "step")==FALSE)

#exclude the gin vivo stress protocol
myfiles<-subset(myfiles, str_detect(myfiles, "01272021")==FALSE)

#read files
sum<-data.frame()

for (i in 1:(NROW(myfiles))) {

# Load file
data <- read.delim(myfiles[i], header = TRUE, sep = "\t", dec = ".", skip=21)

#select colums for analysis
data<-select(data, 1, 10, 14, "Absorption..390.018.", "Absorption..430.103.")
colnames(data)<-c("date", "O2", "CO2", "A390", "A430")

#delete last row which is always off
data<-head(data,-1)

data$date<-as.character(data$date)

#--------------------------------------------------------------------------------------------------
# Calculate Hb-O2 extraction
#--------------------------------------------------------------------------------------------------

#calculate isosbesthic correction
data$iso<-data$A430/data$A390

#calibration values
high_1_iso<-data[(1:500),]
high_1_iso<-na.omit(subset(high_1_iso, data$O2==99.71&data$CO2==0.29))
high_1_time<-rownames(high_1_iso)
high_1_time<-tail(high_1_time, n=10)
high_1_time<-mean(as.numeric(high_1_time))
high_1_iso<-tail(high_1_iso$iso, n=10)
high_1_iso<-mean(high_1_iso)

high_2_iso<-na.omit(subset(data, data$O2==99.71&data$CO2==0.29))
high_2_time<-rownames(high_2_iso)
high_2_time<-tail(high_2_time, n=10)
high_2_time<-mean(as.numeric(high_2_time))
high_2_iso<-tail(high_2_iso$iso, n=10)
high_2_iso<-mean(high_2_iso)

low_1_iso<-na.omit(subset(data, data$O2==0&data$CO2==0.29))
low_1_time<-rownames(low_1_iso)
low_1_time<-tail(low_1_time, n=10)
low_1_time<-mean(as.numeric(low_1_time))
low_1_iso<-tail(low_1_iso$iso, n=10)
low_1_iso<-mean(low_1_iso)

low_2_iso<-na.omit(subset(data, data$O2==0&data$CO2==CO2))
low_2_time<-rownames(low_2_iso)
low_2_time<-tail(low_2_time, n=10)
low_2_time<-mean(as.numeric(low_2_time))
low_2_iso<-tail(low_2_iso$iso, n=10)
low_2_iso<-mean(low_2_iso)


cal_dat<-data.frame(order=c(1,1, 2,2), time=c(high_1_time, high_2_time, low_1_time, low_2_time),
                    abs=c(high_1_iso, high_2_iso, low_1_iso, low_2_iso))
rm(high_1_iso, high_1_time, high_2_iso, high_2_time, low_1_iso, low_1_time, low_2_iso, low_2_time)

high_lm<-lm(abs~time, data=subset(cal_dat, order==1))
summary<-summary(high_lm)
high_inter<-summary$coefficients[1,1]
high_slope<-summary$coefficients[2,1]

low_lm<-lm(abs~time, data=subset(cal_dat, order==2))
summary<-summary(low_lm)
low_inter<-summary$coefficients[1,1]
low_slope<-summary$coefficients[2,1]

#save i if NAs are produced
if(any(is.na(cal_dat$abs)))corrupted_file<-i
if(any(is.nan(cal_dat$abs)))corrupted_file<-i

rm(summary, cal_dat, high_lm, low_lm)


#data points

ifelse(data$O2[1500]==20.4,
       set_O2<-20.4,
       ifelse(data$O2[1500]==3.15,
              set_O2<-3.15,
              set_O2<-3.00))

CO2_0.3_iso<-(subset(data, data$O2==set_O2&data$CO2==0.29))
CO2_0.3_time<-rownames(CO2_0.3_iso)
CO2_0.3_time<-tail(CO2_0.3_time, n=10)
CO2_0.3_time<-mean(as.numeric(CO2_0.3_time))
CO2_0.3_iso<-tail(CO2_0.3_iso$iso, n=10)
CO2_0.3_iso<-mean(CO2_0.3_iso)

CO2_0.5_iso<-(subset(data, data$O2==set_O2&data$CO2==0.48))
CO2_0.5_time<-rownames(CO2_0.5_iso)
CO2_0.5_time<-tail(CO2_0.5_time, n=10)
CO2_0.5_time<-mean(as.numeric(CO2_0.5_time))
CO2_0.5_iso<-tail(CO2_0.5_iso$iso, n=10)
CO2_0.5_iso<-mean(CO2_0.5_iso)

CO2_1_iso<-(subset(data, data$O2==set_O2&data$CO2==0.98))
CO2_1_time<-rownames(CO2_1_iso)
CO2_1_time<-tail(CO2_1_time, n=10)
CO2_1_time<-mean(as.numeric(CO2_1_time))
CO2_1_iso<-tail(CO2_1_iso$iso, n=10)
CO2_1_iso<-mean(CO2_1_iso)

CO2_1.5_iso<-(subset(data, data$O2==set_O2&data$CO2==1.47))
CO2_1.5_time<-rownames(CO2_1.5_iso)
CO2_1.5_time<-tail(CO2_1.5_time, n=10)
CO2_1.5_time<-mean(as.numeric(CO2_1.5_time))
CO2_1.5_iso<-tail(CO2_1.5_iso$iso, n=10)
CO2_1.5_iso<-mean(CO2_1.5_iso)

CO2_2_iso<-(subset(data, data$O2==set_O2&data$CO2==1.95))
CO2_2_time<-rownames(CO2_2_iso)
CO2_2_time<-tail(CO2_2_time, n=10)
CO2_2_time<-mean(as.numeric(CO2_2_time))
CO2_2_iso<-tail(CO2_2_iso$iso, n=10)
CO2_2_iso<-mean(CO2_2_iso)

CO2_3_iso<-(subset(data, data$O2==set_O2&data$CO2==2.9))
CO2_3_time<-rownames(CO2_3_iso)
CO2_3_time<-tail(CO2_3_time, n=10)
CO2_3_time<-mean(as.numeric(CO2_3_time))
CO2_3_iso<-tail(CO2_3_iso$iso, n=10)
CO2_3_iso<-mean(CO2_3_iso)

data<-data.frame(date=data$date[i], CO2=c(0.3, 0.5, 1, 1.5,2,3), set_O2=set_O2,
                 time=c(CO2_0.3_time, CO2_0.5_time, CO2_1_time, CO2_1.5_time, CO2_2_time, CO2_3_time),
                 abs=c(CO2_0.3_iso, CO2_0.5_iso, CO2_1_iso, CO2_1.5_iso, CO2_2_iso, CO2_3_iso))

rm(CO2_0.3_time, CO2_0.5_time, CO2_1_time, CO2_1.5_time, CO2_2_time, CO2_3_time,
   CO2_0.3_iso, CO2_0.5_iso, CO2_1_iso, CO2_1.5_iso, CO2_2_iso, CO2_3_iso, set_O2)

#corrected values

data$oxy<-data$time*high_slope+high_inter
data$deoxy<-data$time*low_slope+low_inter
data$sat<-(data$deoxy-data$abs)/(data$deoxy-data$oxy)

rm(high_inter, high_slope, low_inter, low_slope)

#name treatments
ifelse(str_detect(myfiles[i], "CA1"), 
               data$drug<-"CA1",
       ifelse(str_detect(myfiles[i], "Am1000"), 
              data$drug<-"Am1000",
       ifelse(str_detect(myfiles[i], "ISO10"), 
              data$drug<-"ISO10",
       ifelse(str_detect(myfiles[i], "DMSO"), 
              data$drug<-"DMSO",
              data$drug<-"unknown"))))
#name O2 levels
ifelse(str_detect(myfiles[i], "normoxia"), 
              data$O2<-"normoxia",
              ifelse(str_detect(myfiles[i], "hypoxia"), 
                     data$O2<-"hypoxia",
                     data$O2<-"unknown"))
#name runs
data$run<-i

#name days
ifelse(data$date=="1/13/2021",
       data$day<-1,
       ifelse(data$date=="1/14/2021",
              data$day<-2,
              ifelse(data$date=="1/19/2021",
                     data$day<-3,
                     ifelse(data$date=="1/21/2021",
                            data$day<-4,
                            ifelse(data$date=="1/26/2021",
                                   data$day<-5,
                                   ifelse(data$date=="1/27/2021",
                                          data$day<-7,
                                          ifelse(data$date=="1/28/2021",
                                                 data$day<-6,
                                                 data$day<-NA)))))))

sum<-rbind(sum, data)
}

rm(data, i, myfiles)

#make separate file with in vivo stress data
#sum_all<-sum
#sum<-subset(sum, !(sum$date=="1/27/2021"))

#check for any NAs
sum(is.na(sum))
#sum<-na.omit(sum)

#---------------------------------------------------------------------------------------------
#Figure 3
#---------------------------------------------------------------------------------------------

# make groups by O2
sum_normoxia<-subset(sum, O2=="normoxia")
sum_hypoxia<-subset(sum, O2=="hypoxia")

#relevel drug
sum_normoxia$drug<-factor(sum_normoxia$drug, levels=c("DMSO", "ISO10", "Am1000", "CA1"),
                 labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))
sum_hypoxia$drug<-factor(sum_hypoxia$drug, levels=c("DMSO", "ISO10", "Am1000", "CA1"),
                          labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))
str(sum_hypoxia$drug)

#for results text
mean_DMSO<-aggregate(sum$sat[sum$drug=="DMSO"], by=list(sum$O2[sum$drug=="DMSO"], sum$CO2[sum$drug=="DMSO"]), FUN=mean, na.rm=T)
se_DMSO<-aggregate(sum$sat[sum$drug=="DMSO"], by=list(sum$O2[sum$drug=="DMSO"], sum$CO2[sum$drug=="DMSO"]), FUN=sd, na.rm=T)
n_DMSO<-aggregate(sum$sat[sum$drug=="DMSO"], by=list(sum$O2[sum$drug=="DMSO"], sum$CO2[sum$drug=="DMSO"]), FUN=length)
se_DMSO<-se_DMSO$x/sqrt(n_DMSO$x)
mean_DMSO<-data.frame(CO2= mean_DMSO$Group.1, O2=mean_DMSO$Group.2, mean= mean_DMSO$x,
                         SE=se_DMSO, n=n_DMSO$x)
mean_DMSO

rm(mean_DMSO,se_DMSO,n_DMSO)

mean_DMSO<-aggregate(sum$sat[sum$drug=="DMSO"],
                     by=list(sum$O2[sum$drug=="DMSO"], sum$CO2[sum$drug=="DMSO"]), FUN=mean)

#calculate relative change in Hb-O2 sat
sum_normoxia$sat<-sum_normoxia$sat-1
sum_hypoxia$sat<-sum_hypoxia$sat-1


# fit summary models normoxia
                        
normox_hill_DMSO <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                 data=subset(sum_normoxia, drug=="DMSO"), 
                 start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                 control=nls.control(maxiter=1024))

normox_MM_DMSO <- nlsLM(sat ~ (a*CO2)/(b+CO2), 
               data=subset(sum_normoxia, drug=="DMSO"), 
               start = list(a=-1, b=1), trace=TRUE, 
               control=nls.control(maxiter=1024))

normox_exp_DMSO <- nlsLM(sat ~ (a*exp(b*CO2)), 
                 data=subset(sum_normoxia, drug=="DMSO"), 
                 start = list(a=-1, b=1), trace=TRUE, 
                 control=nls.control(maxiter=1024))

AIC(normox_hill_DMSO, normox_MM_DMSO, normox_exp_DMSO)

rm(normox_MM_DMSO, normox_exp_DMSO)

modsum_normox_DMSO=summary(normox_hill_DMSO)
summary(normox_hill_DMSO)
normox_DMSO<-function(x){(((modsum_normox_DMSO$coefficients[1,1]*x^modsum_normox_DMSO$coefficients[2,1])/
                       ((modsum_normox_DMSO$coefficients[3,1]^modsum_normox_DMSO$coefficients[2,1])+
                          x^modsum_normox_DMSO$coefficients[2,1]))+1)*100}

normox_ISO <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                 data=subset(sum_normoxia, drug=="ISO"), 
                 start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                 control=nls.control(maxiter=1024))

modsum_normox_ISO=summary(normox_ISO)
summary(normox_ISO)
normox_ISO<-function(x){(((modsum_normox_ISO$coefficients[1,1]*x^modsum_normox_ISO$coefficients[2,1])/
                          ((modsum_normox_ISO$coefficients[3,1]^modsum_normox_ISO$coefficients[2,1])+
                             x^modsum_normox_ISO$coefficients[2,1]))+1)*100}

normox_Am <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                  data=subset(sum_normoxia, drug=="ISO+Am"), 
                  start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                  control=nls.control(maxiter=1024))

modsum_normox_Am=summary(normox_Am)
summary(normox_Am)
normox_Am<-function(x){(((modsum_normox_Am$coefficients[1,1]*x^modsum_normox_Am$coefficients[2,1])/
                         ((modsum_normox_Am$coefficients[3,1]^modsum_normox_Am$coefficients[2,1])+
                            x^modsum_normox_Am$coefficients[2,1]))+1)*100}

normox_CA <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                  data=subset(sum_normoxia, drug=="ISO+CA"), 
                  start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                  control=nls.control(maxiter=1024))

modsum_normox_CA=summary(normox_CA)
summary(normox_CA)
normox_CA<-function(x){(((modsum_normox_CA$coefficients[1,1]*x^modsum_normox_CA$coefficients[2,1])/
                         ((modsum_normox_CA$coefficients[3,1]^modsum_normox_CA$coefficients[2,1])+
                            x^modsum_normox_CA$coefficients[2,1]))+1)*100}

normox_halfsat<-data.frame(drug=c("DMSO","ISO", "Am","CA"), 
                     CO2_50=c(modsum_normox_DMSO$coefficients[3,1],
                              modsum_normox_ISO$coefficients[3,1],
                              modsum_normox_Am$coefficients[3,1],
                              modsum_normox_CA$coefficients[3,1]),
                     SEM=c(modsum_normox_DMSO$coefficients[3,2],
                           modsum_normox_ISO$coefficients[3,2],
                           modsum_normox_Am$coefficients[3,2],
                           modsum_normox_CA$coefficients[3,2]))

rm(normox_hill_DMSO)

# fit summary models hypoxia

hypox_hill_DMSO <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                 data=subset(sum_hypoxia, drug=="DMSO"), 
                 start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                 control=nls.control(maxiter=1024))

hypox_MM_DMSO <- nlsLM(sat ~ (a*CO2)/(b+CO2), 
               data=subset(sum_hypoxia, drug=="DMSO"), 
               start = list(a=1, b=1), trace=TRUE, 
               control=nls.control(maxiter=1024))

hypox_exp_DMSO <- nlsLM(sat ~ (a*exp(b*CO2)), 
                data=subset(sum_hypoxia, drug=="DMSO"), 
                start = list(a=1, b=1), trace=TRUE, 
                control=nls.control(maxiter=1024))

AIC(hypox_hill_DMSO, hypox_MM_DMSO, hypox_exp_DMSO)

rm(hypox_MM_DMSO, hypox_exp_DMSO)

modsum_hypox_DMSO=summary(hypox_hill_DMSO)
summary(hypox_hill_DMSO)
hypox_DMSO<-function(x){(((modsum_hypox_DMSO$coefficients[1,1]*x^modsum_hypox_DMSO$coefficients[2,1])/
                           ((modsum_hypox_DMSO$coefficients[3,1]^modsum_hypox_DMSO$coefficients[2,1])+
                              x^modsum_hypox_DMSO$coefficients[2,1]))+1)*100}

hypox_ISO <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                  data=subset(sum_hypoxia, drug=="ISO"), 
                  start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                  control=nls.control(maxiter=1024))

modsum_hypox_ISO=summary(hypox_ISO)
summary(hypox_ISO)
hypox_ISO<-function(x){(((modsum_hypox_ISO$coefficients[1,1]*x^modsum_hypox_ISO$coefficients[2,1])/
                          ((modsum_hypox_ISO$coefficients[3,1]^modsum_hypox_ISO$coefficients[2,1])+
                             x^modsum_hypox_ISO$coefficients[2,1]))+1)*100}

hypox_Am <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                 data=subset(sum_hypoxia, drug=="ISO+Am"), 
                 start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                 control=nls.control(maxiter=1024))

modsum_hypox_Am=summary(hypox_Am)
summary(hypox_Am)
hypox_Am<-function(x){(((modsum_hypox_Am$coefficients[1,1]*x^modsum_hypox_Am$coefficients[2,1])/
                         ((modsum_hypox_Am$coefficients[3,1]^modsum_hypox_Am$coefficients[2,1])+
                            x^modsum_hypox_Am$coefficients[2,1]))+1)*100}

hypox_CA <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                 data=subset(sum_hypoxia, drug=="ISO+CA"), 
                 start = list(a=-0.5, b=2, c=0.8), trace=TRUE, 
                 control=nls.control(maxiter=1024))

modsum_hypox_CA=summary(hypox_CA)
summary(hypox_CA)
hypox_CA<-function(x){(((modsum_hypox_CA$coefficients[1,1]*x^modsum_hypox_CA$coefficients[2,1])/
                         ((modsum_hypox_CA$coefficients[3,1]^modsum_hypox_CA$coefficients[2,1])+
                            x^modsum_hypox_CA$coefficients[2,1]))+1)*100}

hypox_halfsat<-data.frame(drug=c("DMSO","ISO", "Am","CA"), 
                     CO2_50=c(modsum_hypox_DMSO$coefficients[3,1],
                              modsum_hypox_ISO$coefficients[3,1],
                              modsum_hypox_Am$coefficients[3,1],
                              modsum_hypox_CA$coefficients[3,1]),
                     SEM=c(modsum_hypox_DMSO$coefficients[3,2],
                           modsum_hypox_ISO$coefficients[3,2],
                           modsum_hypox_Am$coefficients[3,2],
                           modsum_hypox_CA$coefficients[3,2]))

rm(hypox_hill_DMSO)


#adjust scales for sat

sum_normoxia$sat<-(sum_normoxia$sat+1)*100

sum_hypoxia$sat<-(sum_hypoxia$sat+1)*100

# axis labels

xlab<-(expression(paste("PCO"[2], " (%)")))
ylab<-(expression(paste("Hb-O"[2], " sat. (%)")))

display.brewer.all(colorblindFriendly = TRUE)
col_pal<-brewer.pal(8, "Dark2")


F3<-ggplot(data=sum_normoxia, aes(x=CO2, y=sat, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug",  values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug",  values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1.5, alpha=0.7, shape=21, show.legend = FALSE)+
  geom_jitter(data=sum_hypoxia, position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1.5, alpha=0.7, shape=21, fill="white", show.legend = FALSE)+
  stat_summary(fun=mean, geom="point", shape=21, size=3, show.legend = T)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1, size=0.6, show.legend = T)+
  
  stat_summary(data=sum_hypoxia, fun.data=mean_se, geom="errorbar", width=0.1, size=0.6,
               show.legend = T)+
  stat_summary(data=sum_hypoxia, fun=mean, geom="point", shape=21, size=3,  
               fill="white", show.legend = F)+
  
   stat_function(fun=normox_DMSO, colour=col_pal[3],xlim = c(0.29, 3))+
   stat_function(fun=normox_ISO, colour=col_pal[4], xlim = c(0.29, 3))+
   stat_function(fun=normox_Am, colour=col_pal[5], xlim = c(0.29, 3))+
   stat_function(fun=normox_CA, colour=col_pal[6], xlim = c(0.29, 3))+
  
   stat_function(fun=hypox_DMSO, colour=col_pal[3], xlim = c(0.29, 3))+
   stat_function(fun=hypox_ISO, colour=col_pal[4], xlim = c(0.29, 3))+
   stat_function(fun=hypox_Am, colour=col_pal[5], xlim = c(0.29, 3))+
   stat_function(fun=hypox_CA, colour=col_pal[6], xlim = c(0.29, 3))+
  
   scale_x_continuous(name = xlab, labels=c("0", "", 0.5, 1, 1.5, 2, 3),
                      breaks = c(0, 0.3, 0.5, 1, 1.5, 2, 3), limits=c(0, 3.2))+
   scale_y_continuous(name = ylab, breaks=seq(20, 100, 10), limits=c(20, 100.5))+
  #geom_text(label=pval_drug, x=2.4 , y=100, colour="black", size=3,
  #           show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  # geom_text(label=pval_O2, x=2.4 , y=95, colour="black", size=3,
  #           show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  # geom_text(label=pval_inter, x=2.4 , y=90, colour="black", size=3,
  #           show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 12),
        axis.title.x = element_text(face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 10),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 10),
        strip.text.x = element_text(face="bold", colour="black", size = 10),
        legend.text = element_text(face="italic", colour="black", size = 8),
        legend.title = element_text(colour="black", size = 10),
        legend.position = c(0.9, 0.9),
        legend.margin = NULL,
        legend.key.height=unit(0.5,"line"),
        legend.key.width=unit(0.5,"line"))
F3

ggsave(filename="Figure 3_WSB bNHE_Sat_02042021.jpeg", unit="cm", width=12, height=12, plot=F3)
ggsave(filename="Figure 3_WSB bNHE_Sat_02042021.pdf", unit="cm", width=12, height=12, plot=F3, useDingbats=FALSE)


rm(modsum_hypox_DMSO, modsum_hypox_ISO, modsum_hypox_Am,modsum_hypox_CA, 
   modsum_normox_DMSO, modsum_normox_ISO, modsum_normox_Am,modsum_normox_CA, 
   sum_hypoxia, sum_normoxia, hypox_Am, hypox_CA, hypox_DMSO, hypox_ISO,
   normox_Am, normox_CA, normox_DMSO, normox_ISO, F3, xlab, ylab)

#---------------------------------------------------------------------------------------------
#Figure 4
#---------------------------------------------------------------------------------------------
#Fig 4A Halfmax

#halftime statistical analysis

model_dat<-sum
model_dat$sat<-(model_dat$sat-1)*100


model_params<-data.frame()

for (i in 1:(dim(model_dat)[1]/6)){
  
  data<-subset(model_dat, run==i)
  
  hill <- nlsLM(sat ~ (a*CO2^b)/(c^b+CO2^b), 
                data=data, 
                start = list(a=-50, b=0.5, c=2), trace=TRUE, 
                control=nls.control(maxiter=1024))
  
  modsum_hill=summary(hill)
  
  data<-data.frame(date=data$date[1], drug=data$drug[1], O2=data$O2[1], run=data$run[1], a=modsum_hill$coefficients[1,1],
                   b=modsum_hill$coefficients[2,1], c=modsum_hill$coefficients[3,1])
  
  model_params<-rbind(model_params, data)
  
  rm(data, hill, modsum_hill, i)
  
}

rm(model_dat)

#for results text
mean_halfsat<-aggregate(model_params$c, by=list(model_params$drug, model_params$O2), FUN=mean, na.rm=T)
se_halfsat<-aggregate(model_params$c, by=list(model_params$drug, model_params$O2), FUN=sd, na.rm=T)
n_halfsat<-aggregate(model_params$c, by=list(model_params$drug, model_params$O2), FUN=length)
se_halfsat<-se_halfsat$x/sqrt(n_halfsat$x)
mean_halfsat<-data.frame(CO2= mean_halfsat$Group.1, O2=mean_halfsat$Group.2, mean= mean_halfsat$x,
                        SE=se_halfsat, n=n_halfsat$x)
mean_halfsat

rm(mean_halfsat,se_halfsat,n_halfsat)

model_params$drug<-factor(model_params$drug, levels=c("DMSO", "ISO10", "Am1000", "CA1"),
                             labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))
str(model_params)

#2-way ANOVA

reg_half_t<-lm(c~drug*O2, data=model_params)

#transform to meet homogeneity of variances

model_params$trans_c<-log10(model_params$c)
reg_half_t<-lm(trans_c~drug*O2, data=model_params)

#check parametric assumptions
reg_half_t_res<-data.frame(reg_half_t$residuals)
stat.desc(reg_half_t_res, basic=FALSE, norm=TRUE)
plot(fitted(reg_half_t), residuals(reg_half_t))
qplot(sample=reg_half_t_res$reg_half_t.residuals)
ggplot(reg_half_t_res, aes(reg_half_t.residuals))+
  geom_histogram(aes(y=..density..), binwidth=0.05, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(reg_half_t_res$reg_half_t.residuals, na.rm=TRUE),
                                     sd=sd(reg_half_t_res$reg_half_t.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(model_params$trans_c, model_params$drug, center=median)
leveneTest(model_params$trans_c, model_params$O2, center=median)

anova_reg_half_t<-Anova(reg_half_t, type="II")
anova_reg_half_t

#safe p-values

pval_drug<-anova_reg_half_t$`Pr(>F)`[1]
pval_O2<-anova_reg_half_t$`Pr(>F)`[2]
pval_inter<-anova_reg_half_t$`Pr(>F)`[3]

rm(anova_reg_half_t, reg_half_t, reg_half_t_res)


#significance labels

small_p=0.001

pval_drug<-formatC(round(pval_drug, 3 ), format='f', digits=3)
pval_O2<-formatC(round(pval_O2, 3 ), format='f', digits=3)
pval_inter<-formatC(round(pval_inter, 3 ), format='f', digits=3)

pval_drug<-ifelse(pval_drug < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug))
pval_O2<-ifelse(pval_O2 < small_p, paste("O2: italic(P) < ", small_p),
                sprintf("O2: italic(P)=='%s'", pval_O2))
pval_inter<-ifelse(pval_inter < small_p, paste("drugxO2: italic(P) < ", small_p),
                   sprintf("drugxO2: italic(P)=='%s'", pval_inter))
rm(small_p)

#posthoc analysis
library(multcompView)

posthoc_half_t_normox<-pairwise.t.test(subset(model_params$c, model_params$O2=="normoxia"& 
                                         !model_params$date=="1/13/2021"), 
                                    subset(model_params$drug, model_params$O2=="normoxia"& 
                                             !model_params$date=="1/13/2021"), 
                                    p.adjust.method="BH", paired=T)

posthoc_half_t_normox
#make_cld(posthoc_half_t_normox)

siglab_normoxia<-c("a",  "c", "b",  "c")
#                 DMSO, ISO, ISO+Am, ISO+CA  

posthoc_half_t_hypox<-pairwise.t.test(subset(model_params$c, model_params$O2=="hypoxia"& 
                                         !model_params$date=="1/13/2021"), 
                                subset(model_params$drug, model_params$O2=="hypoxia"& 
                                         !model_params$date=="1/13/2021"), 
                                p.adjust.method="BH", paired=T)

posthoc_half_t_hypox
#make_cld(posthoc_half_t_hypox)

siglab_hypoxia<-c("x",  "z", "y",  "z")
#                DMSO,  ISO,  ISO+Am, ISO+CA  

rm(posthoc_half_t_hypox, posthoc_half_t_normox)

aggregate(model_params$c, by=list(model_params$drug, model_params$O2), FUN=mean)


#axis labels
spacer=0.1
ylab<-(expression(paste("EC"[50],"PCO"[2], " (%)")))

F4A<-ggplot(data=subset(model_params, O2=="normoxia"), aes(x=drug, y=c, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", limits = c("DMSO", "ISO", "ISO+Am", "ISO+CA"),
                    values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", limits = c("DMSO", "ISO", "ISO+Am", "ISO+CA"),
                      values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete(labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))+
  scale_y_continuous(name = ylab, breaks = seq(0, 1.4, 0.2), limits=c(0, 1.5))+
  
  geom_jitter(position=position_jitterdodge(jitter.width=0.5, dodge.width=0.2),
              size=2, alpha=0.5, shape=21, show.legend = FALSE)+
  geom_jitter(data=subset(model_params, O2=="hypoxia"), 
              position=position_jitterdodge(jitter.width=0.5, dodge.width=0.2),
              size=2, alpha=0.5, shape=1, show.legend = FALSE)+
  
  stat_summary(fun=mean, geom="point", shape=21, size=3, show.legend = T)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, size=0.8, show.legend = T)+
  
  stat_summary(data=subset(model_params, O2=="hypoxia"), fun.data=mean_se, geom="errorbar",
               width=0.2, size=0.8, show.legend = T)+
  stat_summary(data=subset(model_params, O2=="hypoxia"), fun=mean, geom="point", shape=21, 
               size=3, fill="white", show.legend = F)+
  
  geom_text(data=normox_halfsat,
            y=c(max(model_params$c[model_params$O2=="normoxia"&model_params$drug=="DMSO"])+spacer,
                max(model_params$c[model_params$O2=="normoxia"&model_params$drug=="ISO"])+spacer,
                max(model_params$c[model_params$O2=="normoxia"&model_params$drug=="ISO+Am"])+spacer,
                max(model_params$c[model_params$O2=="normoxia"&model_params$drug=="ISO+CA"])+spacer),
            x=c(1,2,3,4), colour="black",
            label=siglab_normoxia,
            size=2.5, show.legend = FALSE)+
  
  geom_text(data=hypox_halfsat,
            y=c(max(model_params$c[model_params$O2=="hypoxia"&model_params$drug=="DMSO"])+spacer,
                max(model_params$c[model_params$O2=="hypoxia"&model_params$drug=="ISO"])+spacer,
                max(model_params$c[model_params$O2=="hypoxia"&model_params$drug=="ISO+Am"])+spacer,
                max(model_params$c[model_params$O2=="hypoxia"&model_params$drug=="ISO+CA"])+spacer),
            x=c(1,2,3,4), colour="black",
            label=siglab_hypoxia,
            size=2.5, show.legend = FALSE)+
  
  geom_text(label=pval_drug, x=4.5 , y=0.2, colour="black", size=3,
             show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_O2, x=4.5 , y=0.1, colour="black", size=3,
             show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter, x=4.5 , y=0, colour="black", size=3,
             show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x =  element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = "none",
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
F4A

ggsave(filename="Figure 4A_WSB bNHE_halfsat_02162021.jpeg", unit="cm", width=8, height=12, plot=F4A)
ggsave(filename="Figure 4A_WSB bNHE_halfsat_02162021.pdf", unit="cm", width=8, height=12, plot=F4A, useDingbats=FALSE)

rm(pval_drug, pval_inter, pval_O2, ylab, siglab_normoxia, spacer, normox_halfsat, hypox_halfsat)

#---------------------------------------------------------------------------------------------
#Figure 4B Vmax

#for results text
mean_vmax<-aggregate(model_params$a, by=list(model_params$O2), FUN=mean, na.rm=T)
se_vmax<-aggregate(model_params$a, by=list(model_params$O2), FUN=sd, na.rm=T)
n_vmax<-aggregate(model_params$a, by=list(model_params$O2), FUN=length)
se_vmax<-se_vmax$x/sqrt(n_vmax$x)
mean_vmax<-data.frame(O2=mean_vmax$Group.1, mean= mean_vmax$x,
                         SE=se_vmax, n=n_vmax$x)
mean_vmax

rm(mean_vmax,se_vmax,n_vmax)

#2-way ANOVA

reg_vmax<-lm(a~drug*O2, data=model_params)

#check parametric assumptions
reg_vmax_res<-data.frame(reg_vmax$residuals)
stat.desc(reg_vmax_res, basic=FALSE, norm=TRUE)
plot(fitted(reg_vmax), residuals(reg_vmax))
qplot(sample=reg_vmax_res$reg_vmax.residuals)
ggplot(reg_vmax_res, aes(reg_vmax.residuals))+
  geom_histogram(aes(y=..density..), binwidth=1, colour="Black", fill="White")+
  labs(x="Activity", y="Individual")+
  stat_function(fun=dnorm, args=list(mean=mean(reg_vmax_res$reg_vmax.residuals, na.rm=TRUE),
                                     sd=sd(reg_vmax_res$reg_vmax.residuals, na.rm=TRUE)), colour="Black", size=1)

leveneTest(model_params$a, model_params$drug, center=median)
leveneTest(model_params$a, model_params$O2, center=median)

anova_reg_vmax<-Anova(reg_vmax, type="II")
anova_reg_vmax

#safe p-values

pval_drug_vmax<-anova_reg_vmax$`Pr(>F)`[1]
pval_O2_vmax<-anova_reg_vmax$`Pr(>F)`[2]
pval_inter_vmax<-anova_reg_vmax$`Pr(>F)`[3]

rm(anova_reg_vmax, reg_vmax, reg_vmax_res)

#significance labels

small_p=0.001

pval_drug_vmax<-formatC(round(pval_drug_vmax, 3 ), format='f', digits=3)
pval_O2_vmax<-formatC(round(pval_O2_vmax, 3 ), format='f', digits=3)
pval_inter_vmax<-formatC(round(pval_inter_vmax, 3 ), format='f', digits=3)

pval_drug_vmax<-ifelse(pval_drug_vmax < small_p, paste("drug: italic(P) < ", small_p),
                  sprintf("drug: italic(P)=='%s'", pval_drug_vmax))
pval_O2_vmax<-ifelse(pval_O2_vmax < small_p, paste("O2: italic(P) < ", small_p),
                sprintf("O2: italic(P)=='%s'", pval_O2_vmax))
pval_inter_vmax<-ifelse(pval_inter_vmax < small_p, paste("drugxO2: italic(P) < ", small_p),
                   sprintf("drugxO2: italic(P)=='%s'", pval_inter_vmax))
rm(small_p)

#get Root effect for control treatment

as.character(formatC(round(mean(model_params$a[model_params$drug=="DMSO"&model_params$O2=="normoxia"]), 2 ), format='f', digits=2))
as.character(formatC(round(sd(model_params$a[model_params$drug=="DMSO"&model_params$O2=="normoxia"])/
  sqrt(as.numeric(length(model_params$a[model_params$drug=="DMSO"&model_params$O2=="normoxia"]))), 
  2 ), format='f', digits=2))

#axis labels
ylab<-(expression(paste("Max. ", Delta, "Hb-O"[2], " sat. (%)")))


F4B<-ggplot(data=subset(model_params, O2=="normoxia"), aes(x=drug, y=a, fill=drug, colour=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug",  values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  scale_x_discrete(labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))+
  scale_y_continuous(name = ylab, breaks = seq(-100, 0,  20), limits=c(-100, 0))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.5, dodge.width=0.2),
              size=2, alpha=0.5, shape=21, show.legend = FALSE)+
  geom_jitter(data=subset(model_params, O2=="hypoxia"), 
              position=position_jitterdodge(jitter.width=0.5, dodge.width=0.2),
              size=2, alpha=0.5, shape=1, show.legend = FALSE)+
  
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, size=0.8, show.legend = T)+
  stat_summary(fun=mean, geom="point", shape=21, size=3, show.legend = F)+
  
  stat_summary(data=subset(model_params, O2=="hypoxia"), fun.data=mean_se, geom="errorbar", 
               width=0.2, size=0.8, show.legend = T)+
  stat_summary(data=subset(model_params, O2=="hypoxia"), fun=mean, geom="point", 
               shape=21, fill="white", size=3, show.legend = T)+
  
  geom_text(label=pval_drug_vmax, x=4.5 , y=0, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_O2_vmax, x=4.5 , y=-6, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  geom_text(label=pval_inter_vmax, x=4.5 , y=-12, colour="black", size=3,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
        axis.title.x =  element_blank(),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
        strip.text.x = element_blank(),
        legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.title = element_text(family="sans", face="plain", colour="black", size = 8),
        legend.key.height=unit(0.6,"line"),
        legend.key.width=unit(1,"line"),
        legend.position = "none",
        legend.margin = NULL,
        axis.line = element_line(colour = "black", size=0.5),
        panel.border = element_blank(),
        panel.background = element_blank())
F4B

ggsave(filename="Figure 4B_WSB bNHE_Vmax_02162021.jpeg", unit="cm", width=8, height=12, plot=F4B)
ggsave(filename="Figure 4B_WSB bNHE_Vmax_02162021.pdf", unit="cm", width=8, height=12, plot=F4B, useDingbats=FALSE)

rm(model_params, model_params_hypox, model_params_normox,
   pval_drug_vmax, pval_inter_vmax, pval_O2_vmax, ylab)

#Combined graph
##------------------------------------------------------------------------------------------------
library(cowplot)

g_comb<-ggdraw() +
  draw_plot(F4A, x = 0, y = 1/2, width = 1, height = 1/2)+
  draw_plot(F4B, x = 0, y = 0, width = 1, height = 1/2)+
  
  draw_plot_label(label = c("A", "B"), size = 20,
                  x = c(-0.01, -0.01), y = c(1, 1/2))
g_comb

ggsave(filename="WSB_bNHE_Figure 4_WSB bNHE_kinetics_combined.jpeg", unit="cm", width = 9, height = 14, plot=g_comb)
ggsave(filename="WSB_bNHE_Figure 4_WSB bNHE_kinetics_combined.pdf", unit="cm", width = 9, height = 14, plot=g_comb, useDingbats=FALSE)

rm(F4A, F4B, g_comb, siglab_hypoxia)

#---------------------------------------------------------------------------------------------
#Figure 5
#---------------------------------------------------------------------------------------------

#changes relative to DMSO by day
sum$day<-factor(sum$day)
delta_sat<-data.frame()

for (i in 1:as.numeric(tail(levels(sum$day),n=1))){
  
  rel_hypox<-subset(sum, sum$day==i& sum$O2=="hypoxia")
  rel_hypox$rel_sat<-rel_hypox$sat-
    rel_hypox$sat[rel_hypox$drug=="DMSO"]
  rel_normox<-subset(sum, sum$day==i& sum$O2=="normoxia")
  rel_normox$rel_sat<-rel_normox$sat-
    rel_normox$sat[rel_normox$drug=="DMSO"]  
 
  delta_sat<-rbind(delta_sat, rel_hypox, rel_normox)

}

rm(rel_hypox, rel_normox, i)

delta_sat$rel_sat<-delta_sat$rel_sat*100

#for results text
mean_delta<-aggregate(delta_sat$rel_sat, by=list(delta_sat$O2, delta_sat$CO2, delta_sat$drug), FUN=mean, na.rm=T)
se_delta<-aggregate(delta_sat$rel_sat, by=list(delta_sat$O2, delta_sat$CO2, delta_sat$drug), FUN=sd, na.rm=T)
n_delta<-aggregate(delta_sat$rel_sat, by=list(delta_sat$O2, delta_sat$CO2, delta_sat$drug), FUN=length)
se_delta<-se_delta$x/sqrt(n_delta$x)
mean_delta<-data.frame(O2=mean_delta$Group.1, CO2= mean_delta$Group.2, drug=mean_delta$Group.3,
                       mean= mean_delta$x, SE=se_delta, n=n_delta$x)
mean_delta

rm(mean_delta,se_delta,n_delta)


#relevel drug
delta_sat$drug<-factor(delta_sat$drug, levels=c("DMSO", "ISO10", "Am1000", "CA1"),
                          labels=c("DMSO", "ISO", "ISO+Am", "ISO+CA"))

#absolute changes in delta Hb-O2 Sat.
# axis labels

xlab<-(expression(paste("PCO"[2], " (%)")))
ylab<-(expression(paste(Delta,"Hb-O"[2], " sat. (%)")))

F5A<-ggplot(data=subset(delta_sat, O2=="normoxia"), aes(x=CO2, y=rel_sat, colour=drug, fill=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  geom_smooth(formula=y~x, method="loess", size=0.5, show.legend = F, alpha=0.3)+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1, alpha=0.5, shape=21, show.legend = FALSE)+
  stat_summary(fun=mean, geom="point", shape=21, size=2, show.legend = T)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1, size=0.5, show.legend = T)+
 # stat_function(fun=normox_DMSO, colour="blue",xlim = c(0.29, 3))+
 # stat_function(fun=normox_ISO, colour="red", xlim = c(0.29, 3))+
 # stat_function(fun=normox_Am, colour="yellow", xlim = c(0.29, 3))+
 # stat_function(fun=normox_CA, colour="green", xlim = c(0.29, 3))+
  
  scale_x_continuous(name = xlab, labels=c("0", "", 0.5, 1, 1.5, 2, 3),
                     breaks = c(0, 0.3, 0.5, 1, 1.5, 2, 3), limits=c(0, 3.2))+
  scale_y_continuous(name = ylab, breaks=seq(-2, 12, 2), limits=c(-2, 13))+
  #   coord_cartesian(ylim=c(0.3, 1))+
  #  geom_hline(yintercept = 0, linetype="dashed", size=0.5, colour="black", show.legend = FALSE)+
  #geom_signif(comparisons=list(c("Ctrl","CA")), y_position=0.3, 
  #            annotations=sig_ctrl_ca, textsize = 4, show.legend = FALSE)+
  #  geom_text(label=sig_label, x=2.5 , y=90, colour="black", size=4,
  #            show.legend = FALSE, parse=TRUE, check_overlap = TRUE)+
  #  geom_text(data=sig_letters, aes(x=x, y=y, label=label), colour="black", size=4,
  #            show.legend = FALSE, check_overlap = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 12),
        axis.title.x = element_text(face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 10),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 10),
        strip.text.x = element_text(face="bold", colour="black", size = 10),
        legend.text = element_text(face="italic", colour="black", size = 8),
        legend.title = element_text(colour="black", size = 10),
        legend.position = c(0.9, 0.8),
        legend.margin = NULL,
        legend.key.height=unit(0.5,"line"),
        legend.key.width=unit(0.5,"line"))
F5A


F5B<-ggplot(data=subset(delta_sat, O2=="hypoxia"), aes(x=CO2, y=rel_sat, colour=drug, fill=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  geom_smooth(formula=y~x, method="loess", size=0.5, show.legend = F, alpha=0.3)+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1, alpha=0.5, shape=21, fill="white", show.legend = FALSE)+
  
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1, size=0.5, show.legend = T)+
  stat_summary(fun=mean, geom="point", shape=21, fill="white", size=2, show.legend = T)+
  # stat_function(fun=normox_DMSO, colour="blue",xlim = c(0.29, 3))+
  # stat_function(fun=normox_ISO, colour="red", xlim = c(0.29, 3))+
  # stat_function(fun=normox_Am, colour="yellow", xlim = c(0.29, 3))+
  # stat_function(fun=normox_CA, colour="green", xlim = c(0.29, 3))+
  
  scale_x_continuous(name = xlab, labels=c("0", "", 0.5, 1, 1.5, 2, 3),
                     breaks = c(0, 0.3, 0.5, 1, 1.5, 2, 3), limits=c(0, 3.2))+
  scale_y_continuous(name = ylab, breaks=seq(-2, 12, 2), limits=c(-4, 13))+
  #   coord_cartesian(ylim=c(0.3, 1))+
  #  geom_hline(yintercept = 0, linetype="dashed", size=0.5, colour="black", show.legend = FALSE)+
  #geom_signif(comparisons=list(c("Ctrl","CA")), y_position=0.3, 
  #            annotations=sig_ctrl_ca, textsize = 4, show.legend = FALSE)+
  #  geom_text(label=sig_label, x=2.5 , y=90, colour="black", size=4,
  #            show.legend = FALSE, parse=TRUE, check_overlap = TRUE)+
  #  geom_text(data=sig_letters, aes(x=x, y=y, label=label), colour="black", size=4,
  #            show.legend = FALSE, check_overlap = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 12),
        axis.title.x = element_text(face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 10),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 10),
        strip.text.x = element_text(face="bold", colour="black", size = 10),
        legend.text = element_text(face="italic", colour="black", size = 8),
        legend.title = element_text(colour="black", size = 10),
        legend.position = c(0.9, 0.8),
        legend.margin = NULL,
        legend.key.height=unit(0.5,"line"),
        legend.key.width=unit(0.5,"line"))
F5B


# Relative values to arterial Hb-O2 sat

delta_sat$rel_sat[delta_sat$O2=="normoxia"]<-delta_sat$rel_sat[delta_sat$O2=="normoxia"]/
                         (mean_DMSO$x[2])
delta_sat$rel_sat[delta_sat$O2=="hypoxia"]<-delta_sat$rel_sat[delta_sat$O2=="hypoxia"]/
  (mean_DMSO$x[1])

#for results text
mean_rel_sat<-aggregate(delta_sat$rel_sat, by=list(delta_sat$drug, delta_sat$O2, delta_sat$CO2), FUN=mean, na.rm=T)
se_rel_sat<-aggregate(delta_sat$sat, by=list(delta_sat$drug, delta_sat$O2, delta_sat$CO2), FUN=sd, na.rm=T)
n_rel_sat<-aggregate(delta_sat$sat, by=list(delta_sat$drug, delta_sat$O2, delta_sat$CO2), FUN=length)
se_rel_sat<-se_rel_sat$x/sqrt(n_rel_sat$x)
mean_rel_sat<-data.frame(CO2= mean_rel_sat$Group.1, O2=mean_rel_sat$Group.2, 
                         CO2=mean_rel_sat$Group.3, mean= mean_rel_sat$x,
                         SE=se_rel_sat, n=n_rel_sat$x)
mean_rel_sat

rm(mean_rel_sat,se_rel_sat,n_rel_sat)

mean_rel_sat<-mean(subset(delta_sat$rel_sat, delta_sat$drug=="ISO+Am"&delta_sat$O2=="normoxia"))
se_rel_sat<-sd(subset(delta_sat$rel_sat, delta_sat$drug=="ISO+Am"&delta_sat$O2=="normoxia"))
n_rel_sat<-length(subset(delta_sat$rel_sat, delta_sat$drug=="ISO+Am"&delta_sat$O2=="normoxia"))
se_rel_sat<-se_rel_sat/sqrt(n_rel_sat)
mean_rel_sat
se_rel_sat
rm(mean_rel_sat,se_rel_sat,n_rel_sat)


F5A<-ggplot(data=subset(delta_sat, O2=="normoxia"), aes(x=CO2, y=rel_sat, colour=drug, fill=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1.5, alpha=0.7, shape=21, show.legend = FALSE)+
  stat_summary(fun=mean, geom="point", shape=21, size=2, show.legend = T)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1, size=0.8, show.legend = T)+
  geom_smooth(formula=y~x, method="loess", size=0.5, show.legend = F, alpha=0.2)+
  # stat_function(fun=normox_DMSO, colour="blue",xlim = c(0.29, 3))+
  # stat_function(fun=normox_ISO, colour="red", xlim = c(0.29, 3))+
  # stat_function(fun=normox_Am, colour="yellow", xlim = c(0.29, 3))+
  # stat_function(fun=normox_CA, colour="green", xlim = c(0.29, 3))+
  
  scale_x_continuous(name = xlab, labels=c("0", "", 0.5, 1, 1.5, 2, 3),
                     breaks = c(0, 0.3, 0.5, 1, 1.5, 2, 3), limits=c(0, 3.2))+
  scale_y_continuous(name = ylab, breaks=seq(-2, 12, 2), limits=c(-2, 13))+
  geom_text(label="Normoxia", x=2.2 , y=12.5, colour="black", size=5,
            show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  
  #   coord_cartesian(ylim=c(0.3, 1))+
  #  geom_hline(yintercept = 0, linetype="dashed", size=0.5, colour="black", show.legend = FALSE)+
  #geom_signif(comparisons=list(c("Ctrl","CA")), y_position=0.3, 
  #            annotations=sig_ctrl_ca, textsize = 4, show.legend = FALSE)+
  #  geom_text(label=sig_label, x=2.5 , y=90, colour="black", size=4,
  #            show.legend = FALSE, parse=TRUE, check_overlap = TRUE)+
  #  geom_text(data=sig_letters, aes(x=x, y=y, label=label), colour="black", size=4,
  #            show.legend = FALSE, check_overlap = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 12),
        axis.title.x = element_text(face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 10),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 10),
        strip.text.x = element_text(face="bold", colour="black", size = 10),
        legend.text = element_text(face="italic", colour="black", size = 8),
        legend.title = element_text(colour="black", size = 10),
        legend.position=c(0.88, 0.8),
        legend.margin = NULL,
        legend.key.height=unit(0.5,"line"),
        legend.key.width=unit(0.5,"line"))
F5A

ggsave(filename="Figure 5A_WSB bNHE_delta Hb-Sat_normox_02162021.jpeg", unit="cm", width=9, height=7, plot=F5A)
ggsave(filename="Figure 5A_WSB bNHE_delta Hb-Sat_normox_02162021.pdf", unit="cm", width=9, height=7, plot=F5A, useDingbats=FALSE)

F5B<-ggplot(data=subset(delta_sat, O2=="hypoxia"), aes(x=CO2, y=rel_sat, colour=drug, fill=drug))+
  scale_fill_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5], col_pal[6]))+
  scale_colour_manual(name="Drug", values = c(col_pal[3], col_pal[4], col_pal[5],col_pal[6]))+
  
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.2),
              size=1.5, alpha=0.7, shape=21, fill="white", show.legend = FALSE)+
  
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1, size=0.8, show.legend = T)+
  stat_summary(fun=mean, geom="point", shape=21, fill="white", size=2, show.legend = T)+
  geom_smooth(formula=y~x, method="loess", size=0.5, show.legend = F, alpha=0.2)+
  # stat_function(fun=normox_DMSO, colour="blue",xlim = c(0.29, 3))+
  # stat_function(fun=normox_ISO, colour="red", xlim = c(0.29, 3))+
  # stat_function(fun=normox_Am, colour="yellow", xlim = c(0.29, 3))+
  # stat_function(fun=normox_CA, colour="green", xlim = c(0.29, 3))+
  
  scale_x_continuous(name = xlab, labels=c("0", "", 0.5, 1, 1.5, 2, 3),
                     breaks = c(0, 0.3, 0.5, 1, 1.5, 2, 3), limits=c(0, 3.2))+
  scale_y_continuous(name = ylab, breaks=seq(-5, 25, 5), limits=c(-7, 25))+
  
  geom_text(label="Hypoxia", x=2.2 , y=23.5, colour="black", size=5,
              show.legend = FALSE, parse=TRUE, check_overlap = TRUE, hjust=1)+
  
  #   coord_cartesian(ylim=c(0.3, 1))+
  #  geom_hline(yintercept = 0, linetype="dashed", size=0.5, colour="black", show.legend = FALSE)+
  #geom_signif(comparisons=list(c("Ctrl","CA")), y_position=0.3, 
  #            annotations=sig_ctrl_ca, textsize = 4, show.legend = FALSE)+
  #  geom_text(label=sig_label, x=2.5 , y=90, colour="black", size=4,
  #            show.legend = FALSE, parse=TRUE, check_overlap = TRUE)+
  #  geom_text(data=sig_letters, aes(x=x, y=y, label=label), colour="black", size=4,
  #            show.legend = FALSE, check_overlap = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(family="sans", face="bold", colour="black", size = 12),
        axis.title.x = element_text(face="bold", colour="black", size = 10),
        axis.text.x = element_text(family="sans", face="plain", colour="black", size = 10),
        axis.text.y = element_text(family="sans", face="plain", colour="black", size = 10),
        strip.text.x = element_text(face="bold", colour="black", size = 10),
        legend.text = element_text(face="italic", colour="black", size = 8),
        legend.title = element_text(colour="black", size = 10),
        legend.position =c(0.88, 0.8),
        legend.margin = NULL,
        legend.key.height=unit(0.5,"line"),
        legend.key.width=unit(0.5,"line"))
F5B

ggsave(filename="Figure 5B_WSB bNHE_delta Hb-Sat_hypox_02162021.tiff", unit="cm", width=9, height=7, plot=F5B)
ggsave(filename="Figure 5B_WSB bNHE_delta Hb-Sat_hypox_02162021.pdf", unit="cm", width=9, height=7, plot=F5B, useDingbats=FALSE)

#Combined graph
##------------------------------------------------------------------------------------------------
library(cowplot)

g_comb<-ggdraw() +
  draw_plot(F5A, x = 0, y = 1/2, width = 1, height = 1/2)+
  draw_plot(F5B, x = 0, y = 0, width = 1, height = 1/2)+
  
  draw_plot_label(label = c("A", "B"), size = 20,
                  x = c(0, 0), y = c(1, 1/2))
g_comb

ggsave(filename="WSB_bNHE_Figure 5_WSB bNHE_delta Hb-Sat_combined.jpeg", unit="cm", width = 9, height = 14, plot=g_comb)
ggsave(filename="WSB_bNHE_Figure 5_WSB bNHE_delta Hb-Sat_combined.pdf", unit="cm", width = 9, height = 14, plot=g_comb, useDingbats=FALSE)

rm(g_comb, F5A, F5B, delta_sat, xlab, ylab, mean_DMSO, col_pal)

