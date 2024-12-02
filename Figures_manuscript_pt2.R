##################################################
#Compare MODELS with and without diversity function
#On Spanish databse
##################################################

library(ggplot2)
library(dplyr)
library(maps)
library(scatterpie)
library(ggrepel)
library(readxl)
library(cowplot)
# Load tidyverse library
library(tidyverse)
library(ggpubr)

library(car)
library(FME)
library(coda)
library(ggmcmc)

######################################################################################################################################################

ROOT_DIR<- "/Users/ebruni/Desktop/HOLISOILS/DATI/MICROBES/"
ROOT_DIR_GRAPHS<- "/Users/ebruni/Desktop/HOLISOILS/MANUSCRIPTS/MICROB_DIVERSITY/FIGURES/"
data_mic_all<- read_xlsx(paste0(ROOT_DIR,"data_with_all_diversity.xlsx"))
source(paste0(ROOT_DIR,"calculate_R2_lme.R")) 

data_mic_all$Rh_C_corrected_day<-as.numeric(data_mic_all$Rh_C_corrected_day)

#a<- drop_na(data_mic_all$Rh_C_corrected_day)
######################################################################################################################################################
###################
# Explore the data
###################
#Rename variables for convenience
data_mic_all$Bacteria_Pielou_evenness<-as.numeric(data_mic_all$Eveness_bact)
data_mic_all$Bacteria_Shannon <-as.numeric(data_mic_all$Shannon_bact)
data_mic_all$Bacteria_invsimpson <-as.numeric(data_mic_all$Inv_Simpson_bact)

data_mic_all$Fungi_Pielou_evenness<-as.numeric(data_mic_all$Eveness_fungi)
data_mic_all$Fungi_Shannon <-as.numeric(data_mic_all$Shannon_fungi)
data_mic_all$Fungi_invsimpson <-as.numeric(data_mic_all$Inv_Simpson_fungi)

#Calculate bacteria diversity index normalized by biomass
#data_mic_all$Bacteria_Pielou_evenness_MICbiomass_ratio<-data_mic_all$Bacteria_Pielou_evenness #Eveness not normalized by biomass (even though the variable is called normalized for convenience of data manipulation)

data_mic_all$Bacteria_Pielou_evenness_MICbiomass_ratio<-data_mic_all$Bacteria_Pielou_evenness/data_mic_all$Microbial_Biomass_m
data_mic_all$Bacteria_Shannon_MICbiomass_ratio <-data_mic_all$Bacteria_Shannon/data_mic_all$Microbial_Biomass_m 
data_mic_all$Bacteria_invsimpson_MICbiomass_ratio <-data_mic_all$Bacteria_invsimpson/data_mic_all$Microbial_Biomass_m

#Calculate fungi diversity index normalized by biomass

#data_mic_all$Fungi_Pielou_evenness_MICbiomass_ratio<-data_mic_all$Fungi_Pielou_evenness #Eveness  not normalized by biomass (even though the variable is called normalized for convenience of data manipulation)
data_mic_all$Fungi_Pielou_evenness_MICbiomass_ratio<-data_mic_all$Fungi_Pielou_evenness/data_mic_all$Microbial_Biomass_m
data_mic_all$Fungi_Shannon_MICbiomass_ratio <-data_mic_all$Fungi_Shannon/data_mic_all$Microbial_Biomass_m 
data_mic_all$Fungi_invsimpson_MICbiomass_ratio <-data_mic_all$Fungi_invsimpson/data_mic_all$Microbial_Biomass_m

#Convert text variables to numeric
data_mic_all$Bulk_Density<-as.numeric(data_mic_all$Bulk_Density)
data_mic_all$Mineralizacion_mg_m<-as.numeric(data_mic_all$Mineralizacion_mg_m)
data_mic_all$Shannon_bact<-as.numeric(data_mic_all$Shannon_bact)

#Considering no stones
#C concentration (%)
data_mic_all$C_Org_conc<- data_mic_all$C_Org_kg_m/data_mic_all$Bulk_Density #% = kgC/m2 / g/cm3
#C content (g/kg)
data_mic_all$C_Org_cont<- data_mic_all$C_Org_conc *10

#Considering stones
#SOC (%) = SOC (tC/ha)/[BD(g/cm3)*depth(cm)*(1-vol.coarse (%/100))]
#sampl depth = 10.; conversion from SOC (kgC/m2) to (tC/ha) = 10
#=>SOC (%) = SOC (kgC/m2)/[BD(g/cm3)*(1-vol.coarse (%/100))]
#data_mic_all$C_Org_conc<- data_mic_all$C_Org_kg_m/(data_mic_all$Bulk_Density*(1-data_mic_all$Piedras_salud/100.)) #%

#Convert microbial biomass from gC/m2 (??) to kgC/m2
data_mic_all$Microbial_Biomass_kg_m2 <- data_mic_all$Microbial_Biomass_m/1000. #kgC/m2 ?
data_mic_all$Microbial_Biomass_m_conc<- data_mic_all$Microbial_Biomass_kg_m2/data_mic_all$Bulk_Density #%
#microbial C content (g/kg)
data_mic_all$Microbial_Biomass_m_cont<- data_mic_all$Microbial_Biomass_m_conc *10

#Mineralization
data_mic_all$Mineralizacion_kg_m2_day<-data_mic_all$Mineralizacion_mg_m/1e6 #from mgN/m2/day (?) to kgN/m2/day

#Suppose that Mineralizacion_mg_m is in kgN/m2/day
#data_mic_all$Mineralizacion_kg_m2_day<-data_mic_all$Mineralizacion_mg_m #from kgN/m2/day (?) to kgN/m2/day

data_mic_all$Mineralizacion_conc<-data_mic_all$Mineralizacion_kg_m2_day/data_mic_all$Bulk_Density#from kgN/m2/day to %
data_mic_all$Mineralizacion_cont<-data_mic_all$Mineralizacion_conc*10 #from % to gN/kg/day

#Respiration
#data_mic_all$Avg_Rh_kg_m2_day <-data_mic_all$Avg_Rh_m *4.4e-8 *0.27*24 # from umolCO2/m2/hr to kg C/m2/day
data_mic_all$Avg_Rh_kg_m2_day <-data_mic_all$Avg_Rh_m *4.4e-8 *0.27*24*3600 # from umolCO2/m2/sec to kg C/m2/day
data_mic_all$Avg_Rh_conc<- data_mic_all$Avg_Rh_kg_m2_day/data_mic_all$Bulk_Density ##from kg C/m2/day to %
data_mic_all$Avg_Rh_cont<- data_mic_all$Avg_Rh_conc*10  #from % to gC/kg/day

##############################################
# Aggregate by mean among replicates
##############################################
data_mic <- data_mic_all %>%
  group_by(Site, Type, .drop = FALSE,.add = TRUE) %>%
  summarise(across(everything(), ~ if(is.numeric(.)){mean(., na.rm = TRUE)} else {first(.)}), .groups = 'drop')

data_mic$C_Org_conc_intensity<-ifelse(data_mic$C_Org_conc<=median(data_mic$C_Org_conc,na.rm = TRUE),"low","high")
min_RhC<-min(data_mic$Rh_C_corrected_day,na.rm = TRUE)
max_RhC<-max(data_mic$Rh_C_corrected_day,na.rm = TRUE)

##############################################
respiration_variable<-data_mic$Avg_Rh_cont
mineralization_variable<-data_mic$Mineralizacion_cont
carbon_content_variable<-data_mic$C_Org_cont
micr_biom_content_variable<-data_mic$Microbial_Biomass_m_cont
##############################################


respiration_variable_expr <-expression("Average respiration (g C"~kg^{-1}~day^{-1}~")")
respiration_variable_expr_log<-expression("Average respiration (g C"~kg^{-1}~day^{-1}~"),"~log[10])

mineralization_variable_expr <-expression("Net mineralization (g N"~kg^{-1}~day^{-1}~")")

carbon_content_variable_expr <-expression("SOC content (g C"~kg^{-1}~")")
carbon_content_variable_expr_log <-expression("SOC content (g C"~kg^{-1}~"),"~log[10])

carbon_stock_variable_expr_log <-expression("SOC stock (kg C"~m^{-2}~"),"~log[10])

micr_biom_content_variable_expr <-expression("Microbial biomass (g C"~kg^{-1}~")")
micr_biom_content_variable_expr_log <-expression("Microbial biomass (g C"~kg^{-1}~"),"~log[10])

richness_bactfungi_expr<- expression("Microbial diversity ("~sqrt(S[fungi]~"*"~S[bacteria])~")")
richness_bactfungi_expr_log<- expression("Microbial diversity ("~sqrt(S[fungi]~"*"~S[bacteria])~"),"~log[10])


fdiv_func_expr <-expression(f[FD] == "max{"~e^(beta~"*"~f[0])~",1}")
div_variable_expr <-expression(f[0] == sqrt(S[fungi]~"*"~S[bacteria]))

model0_Rh_expr<-expression("(a) Model 0:"~R[h]==k[u]*C[s])
modelI_Rh_expr<-expression("(b) Model I:"~R[h](C[s])==k[u]*C[s]*C[b]^gamma)
modelII_Rh_expr<-expression("(c) Model II:"~R[h](C[s],diversity)==k[u]*C[s]*f(diversity)*C[b]^gamma)

model0_Min_expr<-expression("(d) Model 0:"~M==frac(R[h], (1-CUE))~"["~frac(1, (C:N)[s])-frac(CUE,(C:N)[b])~"]")
modelI_Min_expr<-expression("(e) Model I:"~M(C[s])==frac(R[h], (1-CUE))~"["~frac(1, (C:N)[s])-frac(CUE,(C:N)[b])~"]")
modelII_Min_expr<-expression("(f) Model II:"~M(C[s],diversity)==frac(R[h], (1-CUE))~"["~frac(1, (C:N)[s])-frac(CUE,(C:N)[b])~"]")


########################################
#Study difference between ecosystems in variables
#To use ANOVA:
#The data of each group must be  normally distributed

#Test normality of data
#Microbial biomass
shapiro.test(data_mic$Microbial_Biomass_m_cont)
shapiro.test(log(data_mic$Microbial_Biomass_m_cont)) #logtransf
plot(density(log(data_mic$Microbial_Biomass_m_cont)))
#p<0.05 => data NOT normal => use non-parametric tests

#Test variances of groups are equal
#(Levene’s test can be used to check this.)
leveneTest(Microbial_Biomass_m_cont~Manejo_Conteo,data = data_mic)
#p<0.05 => variances are NOT equal
ggqqplot(data_mic$Microbial_Biomass_m_cont)

#Test normality of data
#Bacterial richness
shapiro.test(data_mic$Richness_bact)
shapiro.test(log(data_mic$Richness_bact))#logtransf
plot(density(log(data_mic$Richness_bact)))
#p<0.05 => data NOT normal => use non-parametric tests

#Test variances of groups are equal
#(Levene’s test can be used to check this.)
leveneTest(Richness_bact~Manejo_Conteo,data = data_mic)
#p>0.05 => variances are equal
ggqqplot(data_mic$Richness_bact)


#Test normality of data
#Fungal richness
shapiro.test(data_mic$Richness_fungi)
#p>0.05 => data is normal => use anova

#Test variances of groups are equal
#(Levene’s test can be used to check this.)
leveneTest(Richness_fungi~Manejo_Conteo,data = data_mic)
#p>0.05 => variances are equal => use anova
ggqqplot(data_mic$Richness_fungi)

########################################
# TEST SIGNIFICANT DIFFERENCES
########################################
#Microbial biomass
# Visualize summary of groups
group_by(data_mic, Manejo_Conteo) %>%
  summarise(
    count = n(),
    mean = mean(Microbial_Biomass_m_cont, na.rm = TRUE),
    sd = sd(Microbial_Biomass_m_cont, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Microbial_Biomass_m_cont ~ Manejo_Conteo, data = data_mic)
# Summary of the analysis
summary(res.aov)
#p-value <0.05, => there are significant differences between the groups with “*"
#Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
#the difference between dehesa and the other two groups is significant (adjusted p-value <0.05)


#Non-parametric alternative to one-way ANOVA test for more than 2 samples
kruskal_result_bio<-kruskal.test(Microbial_Biomass_m_cont ~ Manejo_Conteo, data = data_mic)
#p-value <0.05, => there are significant differences between the groups with “*"

# Extract Test Statistic and p-value
krbstat<-as.numeric(kruskal_result_bio$statistic)  # Chi-squared statistic (test statistic)
krbpval<-kruskal_result_bio$p.value    # p-value

#Non-param Multiple pairwise-comparison between groups
pairwise.wilcox.test(data_mic$Microbial_Biomass_m_cont, data_mic$Manejo_Conteo,
                     p.adjust.method = "BH")
#the difference between dehesa and the other two groups is significant (adjusted p-value <0.05)


#Bacterial richness

# Visualize summary of groups
group_by(data_mic, Manejo_Conteo) %>%
  summarise(
    count = n(),
    mean = mean(Richness_bact, na.rm = TRUE),
    sd = sd(Richness_bact, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Richness_bact ~ Manejo_Conteo, data = data_mic)
# Summary of the analysis
summary(res.aov)
#p-value <0.05, => there are significant differences between the groups with “*"
#Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
#the difference between Open woodland and Dehesa is significant (adjusted p-value <0.05)

#Non-parametric alternative to one-way ANOVA test for more than 2 samples
kruskal_result_rb<-kruskal.test(Richness_bact ~ Manejo_Conteo, data = data_mic)
#p-value <0.05, => there are significant differences between the groups with “*"

# Extract Test Statistic and p-value
krRBstat<-kruskal_result_rb$statistic  # Chi-squared statistic (test statistic)
krRBpval<-kruskal_result_rb$p.value    # p-value

#Non-param Multiple pairwise-comparison between groups
pairwise.wilcox.test(data_mic$Richness_bact, data_mic$Manejo_Conteo,
                     p.adjust.method = "BH")
#the difference between dehesa and the other two groups is significant (adjusted p-value <0.05)



#Fungal richness

# Visualize summary of groups
group_by(data_mic, Manejo_Conteo) %>%
  summarise(
    count = n(),
    mean = mean(Richness_fungi, na.rm = TRUE),
    sd = sd(Richness_fungi, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Richness_fungi ~ Manejo_Conteo, data = data_mic)
# Summary of the analysis
aov_result_rf<-summary(res.aov)
#p-value >0.05, => there are NO significant differences between the groups

# Extract Test Statistic and p-value
aovRFstat<-aov_result_rf[[1]]$`F value`[1]  # F (test statistic)
aovRFpval<-aov_result_rf[[1]]$`Pr(>F)`[1]    # p-value

#Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
#the difference between Open woodland and Dehesa is significant (adjusted p-value <0.05)

#Non-parametric alternative to one-way ANOVA test for more than 2 samples
kruskal.test(Richness_fungi ~ Manejo_Conteo, data = data_mic)
#p-value >0.05, => there are NO significant differences between the groups

#Non-param Multiple pairwise-comparison between groups
pairwise.wilcox.test(data_mic$Richness_fungi, data_mic$Manejo_Conteo,
                     p.adjust.method = "BH")
#p-value >0.05, => there are NO significant differences between the groups


#Combined Fungal and bacterial richness

# Visualize summary of groups
group_by(data_mic, Manejo_Conteo) %>%
  summarise(
    count = n(),
    mean = mean(sqrt(Richness_fungi*Richness_bact), na.rm = TRUE),
    sd = sd(sqrt(Richness_fungi*Richness_bact), na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(sqrt(Richness_fungi*Richness_bact) ~ Manejo_Conteo, data = data_mic)
# Summary of the analysis
summary(res.aov)
#p-value >0.05, => there are NO significant differences between the groups
#Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
#there are no differences between groups (adjusted p-value >0.05)

#Non-parametric alternative to one-way ANOVA test for more than 2 samples
kruskal.test(sqrt(Richness_fungi*Richness_bact) ~ Manejo_Conteo, data = data_mic)
#p-value >0.05, => there are NO significant differences between the groups

#Non-param Multiple pairwise-comparison between groups
pairwise.wilcox.test(sqrt(data_mic$Richness_fungi*data_mic$Richness_bact), data_mic$Manejo_Conteo,
                     p.adjust.method = "BH")
#p-value >0.05, => there are NO significant differences between the groups

########################################

labelp1s<-bquote(chi^2== .(round(krbstat,2)))
labelp1p<-bquote(italic(p)== .(round(krbpval,2)))

labelp2s<-bquote(chi^2== .(round(krRBstat,2)))
labelp2p<-bquote(italic(p)== .(round(krRBpval,2)))

labelp3s<-bquote(F== .(round(aovRFstat,2)))
labelp3p<-bquote(italic(p)== .(round(aovRFpval,2)))
#################################
#boxplot version with combined fungal and bacterial richness
p1<-ggplot(data_mic,aes(x=Manejo_Conteo,y=micr_biom_content_variable)) + 
  geom_boxplot(aes(x=Manejo_Conteo,y=micr_biom_content_variable,fill=factor(Manejo_Conteo)))+
  geom_jitter(color="black",size = 2, width = 0.3, alpha = 0.5)+
  #scale_color_manual(values = c("Forest" = "#517664", "Dehesa" = "#F9B5AC", "Open woodland" = "#E69F00")) +  # Color points similarly to fill
  scale_fill_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                             "Open woodland"="#E69F00"))+
  labs(y=micr_biom_content_variable_expr)+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.title.x = element_blank(),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"))+
  ggtitle("(a)")+
  annotate("text", x = 1, y = max(micr_biom_content_variable) + 0.3, label = "a", size = 8) +
  annotate("text", x = 2, y = max(micr_biom_content_variable) + 0.3, label = "b", size = 8) +
  annotate("text", x = 3, y = max(micr_biom_content_variable) + 0.3, label = "b", size = 8)+
  annotate("text", x = 0.8, y = 5, label = labelp1s, size = 6)+
  annotate("text", x = 0.8, y = 4.2, label = labelp1p, size = 6)

p2<-ggplot(data_mic,aes(x=Manejo_Conteo,y=sqrt(Richness_bact*Richness_fungi))) + 
  geom_boxplot(aes(x=Manejo_Conteo,y=sqrt(Richness_bact*Richness_fungi),fill=factor(Manejo_Conteo)))+
  geom_jitter(color="black",size = 2, width = 0.3, alpha = 0.5)+
  #scale_color_manual(values = c("Forest" = "#517664", "Dehesa" = "#F9B5AC", "Open woodland" = "#E69F00")) +  # Color points similarly to fill
  scale_fill_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                             "Open woodland"="#E69F00"))+
  labs(x = "Land-use", y=richness_bactfungi_expr)+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"))+
  ggtitle("(b)")+
  annotate("text", x = 1, y = max(sqrt(data_mic$Richness_fungi*data_mic$Richness_bact)) + 0.3, label = "a", size = 8) +
  annotate("text", x = 2, y = max(sqrt(data_mic$Richness_fungi*data_mic$Richness_bact)) + 0.3, label = "a", size = 8) +
  annotate("text", x = 3, y = max(sqrt(data_mic$Richness_fungi*data_mic$Richness_bact)) + 0.3, label = "a", size = 8)

# Combine the plots into a 2x2 grid and add the shared legend
combined_plot2 <- plot_grid(p1,
                            p2,
                            ncol = 1, align = "hv")
print(combined_plot2)


#################################
#boxplot version with separated fungal and bacterial richness

p3<-ggplot(data_mic,aes(x=Manejo_Conteo,y=Richness_bact)) + 
  geom_boxplot(aes(x=Manejo_Conteo,y=Richness_bact,fill=factor(Manejo_Conteo)))+
  geom_jitter(color="black",size = 2, width = 0.3, alpha = 0.5)+
  #scale_color_manual(values = c("Forest" = "#517664", "Dehesa" = "#F9B5AC", "Open woodland" = "#E69F00")) +  # Color points similarly to fill
  scale_fill_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                             "Open woodland"="#E69F00"))+
  labs(x = " ", y="Bacterial richness")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.title.x = element_blank(),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"))+
  ggtitle("(b)")+
  annotate("text", x = 1, y = max(data_mic$Richness_bact) + 0.3, label = "a", size = 8) +
  annotate("text", x = 2, y = max(data_mic$Richness_bact) + 0.3, label = "b", size = 8) +
  annotate("text", x = 3, y = max(data_mic$Richness_bact) + 0.3, label = "b", size = 8) +
  annotate("text", x = 0.8, y = 950, label=labelp2s, size = 6)+
  annotate("text", x = 0.8, y = 870, label=labelp2p, size = 6)

p4<-ggplot(data_mic,aes(x=Manejo_Conteo,y=Richness_fungi)) + 
  geom_boxplot(aes(x=Manejo_Conteo,y=Richness_fungi,fill=factor(Manejo_Conteo)))+
  geom_jitter(color="black",size = 2, width = 0.3, alpha = 0.5)+
  #scale_color_manual(values = c("Forest" = "#517664", "Dehesa" = "#F9B5AC", "Open woodland" = "#E69F00")) +  # Color points similarly to fill
  scale_fill_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                             "Open woodland"="#E69F00"))+
  labs(x = "Land-use", y="Fungal richness")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"))+
  ggtitle("(c)")+
  annotate("text", x = 1, y = max(data_mic$Richness_fungi) + 0.3, label = "a", size = 8) +
  annotate("text", x = 2, y = max(data_mic$Richness_fungi) + 0.3, label = "a", size = 8) +
  annotate("text", x = 3, y = max(data_mic$Richness_fungi) + 0.3, label = "a", size = 8) +
  annotate("text", x = 0.8, y = 140, label=labelp3s, size = 6)+
  annotate("text", x = 0.8, y = 110, label=labelp3p, size = 6)

# Combine the plots into a 2x2 grid and add the shared legend
combined_plot3 <- plot_grid(p1,
                            p3, p4,
                            ncol = 1, align = "hv")
print(combined_plot3)


ggsave(plot=combined_plot3,paste0(ROOT_DIR_GRAPHS,"Figure_3.pdf"),
       width = 6, height = 13)



##################################################
#Compare MODELS WITH NITROGEN
#with and without microbial biomass 
# and with and without diversity function
##################################################

######################################################################################################################################################
#Diversity function
#f(diversity)=max(e^(beta*f0),1), f0=diversity_variable

#Respiration
#Rh=(1-CUE)*k0*f(diversity)*SOC*Cb^gamma

#Mineralization
#M=Rh/(1-CUE)* [1/CNs - CUE/CNb]

diversity_function<-function(beta_coeff,diversity_variable){
  #f(diversity)=e^(beta*f0), f0=diversity_variable
  func_unlimited<-exp(beta_coeff*diversity_variable)

  return(func_unlimited)
}

Rh_general_model<-function(k_pred,gamma_coeff,beta_coeff,diversity_variable,CUE_micr){
  #k_div = diversity_function
  #If beta = 0 => no effect of diversity
  #If gamma = 0 => no effect of mirob biomass
  
  #kdiv=max(k0*diversity_function,1)
  k_div_unlimited<-k_pred*diversity_function(beta_coeff,diversity_variable)
  # print("k_div_unlimited")
  # print(k_div_unlimited)
  k_div<-k_div_unlimited
  k_div[k_div > 1] <-1
  # print("k_div")
  # print(k_div)
  #Respiration
  #Rh=(1-CUE)*k0*f(diversity)*SOC*Cb^gamma
  Rh_pred<-carbon_content_variable*(1-CUE_micr)*micr_biom_content_variable^gamma_coeff*k_div
  
  return(Rh_pred)
}

Min_general_model<-function(Rhmodel,CUE_micr,CN_ratio_micr){
  #Mineralization
  #M=Rh/(1-CUE)* [1/CNs - CUE/CNb]
  Min_pred<- (Rhmodel/(1-CUE_micr))*((1/data_mic$C_N)-(CUE_micr/CN_ratio_micr))
  
  return(Min_pred)
}

#Model 0:
## gamma = 0
## beta =0
## optimized param: CUE, k0

#Model 1: biomass
## beta =0
## optimized param: CUE, k0, gamma

#Model 2: diversity + biomass
## optimized param: CUE, k0, gamma, beta

#Model 3: diversity
## optimized param: CUE, k0, beta

#NB: in manuscript model 2 = model 3 and viceversa

J_func_model0_MLE <- function(params){
  k_pred<-params[1]
  CUE_micr<-params[2]
  sd_in<-params[length(params)]
  #CN_ratio_micr=7
  
  gamma_coeff=0
  beta_coeff=0
  
  predicted_Rh_model_0<-Rh_general_model(k_pred,gamma_coeff,beta_coeff,diversity_variable_in,CUE_micr)
  predicted_Min_model_0<-Min_general_model(predicted_Rh_model_0,CUE_micr,CN_ratio_micr)
  
  
  obs_Rh_Min<-c(respiration_variable,mineralization_variable)
  
  # Negative log-likelihood 
  negloglike<--2*sum(dnorm(x = obs_Rh_Min, mean = c(predicted_Rh_model_0,predicted_Min_model_0),
                           sd = sd_in, log = TRUE))
  
  
  # print("Negologlike")
  # print(negloglike)
  # print("params")
  # print(params)
  
  
  return(negloglike)
}

J_func_modelI_MLE <- function(params){
  k_pred<-params[1]
  gamma_coeff<-params[2]
  CUE_micr<-params[3]
  sd_in<-params[length(params)]
  #CN_ratio_micr=7
  
  beta_coeff=0
  
  predicted_Rh_model_I<-Rh_general_model(k_pred,gamma_coeff,beta_coeff,diversity_variable_in,CUE_micr)
  predicted_Min_model_I<-Min_general_model(predicted_Rh_model_I,CUE_micr,CN_ratio_micr)
  
  
  obs_Rh_Min<-c(respiration_variable,mineralization_variable)
  
  # Negative log-likelihood 
  negloglike<--2*sum(dnorm(x = obs_Rh_Min, mean = c(predicted_Rh_model_I,predicted_Min_model_I),
                           sd = sd_in, log = TRUE))
  
  
  # print("Negologlike")
  # print(negloglike)
  # print("params")
  # print(params)
  
  
  return(negloglike)
}


J_func_modelII_MLE <- function(params){
  # params<-params_prior_modelII_MLE
  #params<-c(0.01, 0.50, 0.01, 0.10, 0.01)
  
  k_pred<-params[1]
  gamma_coeff<-params[2]
  beta_coeff<-params[3]
  CUE_micr<-params[4]
  sd_in<-params[length(params)]
  #CN_ratio_micr=7
  
  predicted_Rh_model_II<-Rh_general_model(k_pred,gamma_coeff,beta_coeff,diversity_variable_in,CUE_micr)
  # print("predicted_Rh_model_II")
  # print(predicted_Rh_model_II)
  predicted_Min_model_II<-Min_general_model(predicted_Rh_model_II,CUE_micr,CN_ratio_micr)
  # print("predicted_Min_model_II")
  # print(predicted_Min_model_II)
  
  obs_Rh_Min<-c(respiration_variable,mineralization_variable)
  
  # Negative log-likelihood 
  negloglike<--2*sum(dnorm(x = obs_Rh_Min, mean = c(predicted_Rh_model_II,predicted_Min_model_II),
                           sd = sd_in, log = TRUE))
  
  
  # print("Negologlike")
  # print(negloglike)
  # print("params")
  # print(params)
  
  
  return(negloglike)
}


J_func_modelIII_MLE <- function(params){
  #only diversity
  
  k_pred<-params[1]
  beta_coeff<-params[2]
  CUE_micr<-params[3]
  sd_in<-params[length(params)]
  
  gamma_coeff<-0
  #CN_ratio_micr=7
  
  predicted_Rh_model_III<-Rh_general_model(k_pred,gamma_coeff,beta_coeff,diversity_variable_in,CUE_micr)
  # print("predicted_Rh_model_III")
  # print(predicted_Rh_model_III)
  predicted_Min_model_III<-Min_general_model(predicted_Rh_model_III,CUE_micr,CN_ratio_micr)
  # print("predicted_Min_model_III")
  # print(predicted_Min_model_III)
  
  obs_Rh_Min<-c(respiration_variable,mineralization_variable)
  
  # Negative log-likelihood 
  negloglike<--2*sum(dnorm(x = obs_Rh_Min, mean = c(predicted_Rh_model_III,predicted_Min_model_III),
                           sd = sd_in, log = TRUE))
  
  
  # print("Negologlike")
  # print(negloglike)
  # print("params")
  # print(params)
  
  
  return(negloglike)
}

########################################################################
#Optimize with MLE to calculate BIC and Bayes
########################################################################

#Fixed parameters
CN_ratio_micr<-7. #CN ratio of microbes fixed to 7

##################################
#LIST OF DIVERSITY INDICES
#################################
diversity_variable_str_list<-c("sqrt(data_mic$Richness_fungi*data_mic$Richness_bact)",
                               "sqrt(data_mic$Shannon_fungi*data_mic$Shannon_bact)",
                               "sqrt(data_mic$Inv_Simpson_fungi*data_mic$Inv_Simpson_bact)",
                               "data_mic$Richness_fungi",
                               "data_mic$Shannon_fungi",
                               "data_mic$Inv_Simpson_fungi",
                               "data_mic$Richness_bact",
                               "data_mic$Shannon_bact",
                               "data_mic$Inv_Simpson_bact",
                               "data_mic$Observed_sapro",
                               "data_mic$Shannon_sapro",
                               "data_mic$InvSimpson_sapro",
                               "data_mic$Observed_ecto",
                               "data_mic$Shannon_ecto",
                               "data_mic$InvSimpson_ecto")

#Set max bound for priors 
#for all diversity indices
#since they have different numerical values
#they are upper bounded to different limits 
#(the function is exponential)
max_beta_prior_list<-rep(c(0.001,1,1),5)

#Start optimization 
#loop over all diversity indices
j<-1
for(diversity_variable_str in diversity_variable_str_list){
  
  #get the values of the diverisity index
  diversity_variable_in <- eval(parse(text=diversity_variable_str)) #Diversity variable used in model II
  #create a name of the diversity index for plots
  diversity_variable_name<-gsub("data_mic\\$", "", diversity_variable_str)
  
  print(paste0("Running optimization for variable: ",diversity_variable_name))
  
  #Set priors and bounds for optimization
  #Priors
  k_prior<-0.01
  gamma_prior<-0.5
  beta_prior<-0.0001
  CUE_prior<-0.1
  sd_prior<-0.01
  
  #Bounds
  min_k_prior<-1e-6
  min_gamma_prior<-0
  min_beta_prior<--0.1
  min_CUE_prior<-1e-6
  min_sd_prior<-1e-6
  
  max_k_prior<-1
  max_gamma_prior<-1
  #max_beta_prior<-0.001 #when have richness this needs to be low otherwise exp --> inf
  #max_beta_prior<-1 #when have shannon or inv simpson you can increase max beta
  max_beta_prior<-max_beta_prior_list[j]
  max_CUE_prior<-(1-1e-6)
  max_sd_prior<-1
  
  if(j==1){
    #Model 0
    params_prior_model0_MLE<-c(k_prior,CUE_prior,sd_prior)
    lower_lim_model0<-c(min_k_prior,min_CUE_prior,min_sd_prior)
    upper_lim_model0<-c(max_k_prior,max_CUE_prior,max_sd_prior)
    
    fit_model0 = optim(par = params_prior_model0_MLE, fn = J_func_model0_MLE, 
                       method = "L-BFGS-B",
                       lower = lower_lim_model0, 
                       upper = upper_lim_model0,
                       control=list(maxit=1000,parscale=abs(params_prior_model0_MLE)))
    print("--")
    print(paste0("Running model 0: ",fit_model0$value))
    
    #fit_model0
    fit_model0$par
    fit_model0$value
    if(fit_model0$convergence==0){
      print("convergence: yes")
    }else{
      print("convergence issues")
    }
    print("--")
    
    #Model I
    params_prior_modelI_MLE<-c(k_prior,gamma_prior,CUE_prior,sd_prior)
    lower_lim_modelI<-c(min_k_prior,min_gamma_prior,min_CUE_prior,min_sd_prior)
    upper_lim_modelI<-c(max_k_prior,max_gamma_prior,max_CUE_prior,max_sd_prior)
    
    fit_modelI = optim(par = params_prior_modelI_MLE, fn = J_func_modelI_MLE, 
                       method = "L-BFGS-B",
                       lower = lower_lim_modelI, 
                       upper = upper_lim_modelI,
                       control=list(maxit=1000,parscale=abs(params_prior_modelI_MLE)))
    print("--")
    print(paste0("Running model I: ",fit_modelI$value))
    #fit_modelI
    fit_modelI$par
    fit_modelI$value
    
    if(fit_modelI$convergence==0){
      print("convergence: yes")
    }else{
      print("convergence issues")
    }
  }
  
  
  print("--")
  #Model II
  params_prior_modelII_MLE<-c(k_prior, gamma_prior, beta_prior, CUE_prior, sd_prior)
  lower_lim_modelII<-c(min_k_prior, min_gamma_prior, min_beta_prior, min_CUE_prior,min_sd_prior)
  upper_lim_modelII<-c(max_k_prior, max_gamma_prior, max_beta_prior, max_CUE_prior,max_sd_prior)
  
  fit_modelII = optim(par = params_prior_modelII_MLE, fn = J_func_modelII_MLE,
                      method = "L-BFGS-B",
                      lower = lower_lim_modelII,
                      upper = upper_lim_modelII,
                      control=list(maxit=1000,parscale=abs(params_prior_modelII_MLE)))
  
  print("--")
  print(paste0("Running model II: ",fit_modelII$value))
  #fit_modelII
  fit_modelII$par
  fit_modelII$value
  if(fit_modelII$convergence==0){
    print("convergence: yes")
  }else{
    print("convergence issues")
  }
  
  print("--")
  
  #Model III
  params_prior_modelIII_MLE<-c(k_prior,beta_prior,CUE_prior,sd_prior)
  lower_lim_modelIII<-c(min_k_prior,min_beta_prior,min_CUE_prior,min_sd_prior)
  upper_lim_modelIII<-c(max_k_prior,max_beta_prior,max_CUE_prior,max_sd_prior)
  
  fit_modelIII = optim(par = params_prior_modelIII_MLE, fn = J_func_modelIII_MLE, 
                       method = "L-BFGS-B",
                       lower = lower_lim_modelIII, 
                       upper = upper_lim_modelIII,
                       control=list(maxit=1000,parscale=abs(params_prior_modelIII_MLE)))
  print("--")
  print(paste0("Running model III: Jnew=",fit_modelIII$value))
  #fit_modelIII
  fit_modelIII$par
  fit_modelIII$value
  
  if(fit_modelIII$convergence==0){
    print("convergence: yes")
  }else{
    print("convergence issues")
  }
  
  
  print("--")
  print("--------------------------")
  print("END deterministic optimization")
  print("--------------------------")
  
  #Set parameters for Bayesian optimization
  niter_Bayes=150000 #150000
  outputlength_n=1000
  burninlength_Bayes<-0.1*niter_Bayes
  updatecov_Bayes<-100
  ntrydr_Bayes<-2
  
  if(j==1){
    
    ########################################################################
    #Bayesian optim model 0
    ########################################################################
    
    ptm <- proc.time() #start clock to measure time
    
    print("#############################")  
    print("Start Bayes optim for model 0")
    print("#############################")  
    
    #Estimation of the best-fit parameters as a starting point
    iniParPORT_model0 <- fit_model0$par
    
    #Running 3 chains for convergence assessment
    iniParPORT_model0 <- data.frame(Chain1 = iniParPORT_model0,
                                    Chain2 = iniParPORT_model0,
                                    Chain3 = iniParPORT_model0,
                                    row.names = c("k_param","CUE","sd"))
    
    iniParPORT_model0 <- sweep(iniParPORT_model0, MARGIN = 2, STATS = c(1, 0.9, 1.1), FUN = "*")
    
    #Check limits
    for(row in 1:nrow(iniParPORT_model0)){
      iniParPORT_model0[row,][iniParPORT_model0[row,] < lower_lim_model0[row]] <- lower_lim_model0[row]
      iniParPORT_model0[row,][iniParPORT_model0[row,] > upper_lim_model0[row]] <- upper_lim_model0[row]
    }
    
    mcmcDRAM_model0 <- apply(iniParPORT_model0, MARGIN = 2, FUN = function(iIniParPORT_model0) {
      FME::modMCMC(f            = J_func_model0_MLE,
                   p            = iIniParPORT_model0,
                   lower        = lower_lim_model0, ## lower bounds for GR4J
                   upper        = upper_lim_model0, ## upper bounds for GR4J
                   niter        = niter_Bayes, #at least >5000
                   jump         = NULL, #0.01
                   outputlength = outputlength_n, 
                   burninlength = burninlength_Bayes, #needs to be >0 (typical 5-10% of iterations)
                   updatecov    = updatecov_Bayes, ## adaptative Metropolis (AM)
                   ntrydr       = ntrydr_Bayes)   ## delayed rejection (RD)
    })
    
    print("Elapsed time (minutes)")
    #print(proc.time() - ptm) #end and print time elapsed (seconds)
    print((proc.time() - ptm)/60) #minutes
    
    #summary(mcmcDRAM)
    
    #MCMC diagnostics and visualisation tools
    
    multDRAM_model0 <- coda::as.mcmc.list(lapply(mcmcDRAM_model0, FUN = function(x) {
      coda::as.mcmc(as.matrix(x$pars))
    }))
    
    gelRub_model0 <- coda::gelman.diag(multDRAM_model0, autoburnin = TRUE)$mpsrf
    
    print("Gelman Rubin")
    print(gelRub_model0) #if gelRub < 1.1 you can consider convergence
    
    if(gelRub_model0>1.1){
      warning("Careful: chains have not converged to the same posterior")
    }
    
    #Convergence
    parDRAM_model0 <- ggmcmc::ggs(multDRAM_model0) ## to convert object for using by all ggs_* graphical functions
    ggmcmc::ggs_traceplot(parDRAM_model0)
    
    #posterior density parameters
    burnParDRAM_model0 <- parDRAM_model0[parDRAM_model0$Iteration > 500, ] # to keep only the second half of the series
    
    burnParDRAM_model0$Parameter <- gsub("sd", "sigma", burnParDRAM_model0$Parameter) 
    burnParDRAM_model0$Parameter <- gsub("k_param", expression(k[0]), burnParDRAM_model0$Parameter) 
    burnParDRAM_model0$Parameter<-factor(burnParDRAM_model0$Parameter,levels = c("CUE", expression(k[0]),"sigma"))
    
    ggmcmc::ggs_density(burnParDRAM_model0, greek=TRUE)+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            legend.title = element_text(size=14),
            legend.text = element_text(size=14),
            strip.text = element_text(size = 14),
            panel.background = element_rect(fill = "white", color = "white"),
            panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
            panel.grid.minor = element_line(color = "white"))+
      ylab("Density")+xlab("Parameter value")
    
    #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_postPDF_model0_withN_",diversity_variable_name,".png"))
    
    #Parameter correlation
    ggmcmc::ggs_pairs(burnParDRAM_model0, lower = list(continuous = "density"))
    #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_corr_model0_withN_",diversity_variable_name,".png"))
    #----------
    
    #With mean parameter values across chains
    Optim_param_model0_out=t(aggregate(parDRAM_model0$value, list(parDRAM_model0$Parameter), FUN=mean))
    order_pars_out_model0<-Optim_param_model0_out["Group.1",]
    Optim_param_model0_out=as.numeric(Optim_param_model0_out["x",])
    
    Optim_param_std_model0_out=t(aggregate(parDRAM_model0$value, list(parDRAM_model0$Parameter), FUN=sd))
    order_pars_out_sd_model0<-Optim_param_std_model0_out["Group.1",]
    Optim_param_std_model0_out=as.numeric(Optim_param_std_model0_out["x",])
    
    print("Mean optimum parameters - Model 0")
    print(Optim_param_model0_out)
    print("Std optimum parameters - Model 0")
    print(Optim_param_std_model0_out)
    
    CI95_margin_model0<-(qt(0.975,df=outputlength_n-1)*Optim_param_std_model0_out/sqrt(outputlength_n))
    CI95_lower_model0<-Optim_param_model0_out-CI95_margin_model0
    CI95_upper_model0<-Optim_param_model0_out+CI95_margin_model0
    
    mean_bestfunp_model0<-mean(mcmcDRAM_model0$Chain1$bestfunp,
                               mcmcDRAM_model0$Chain2$bestfunp,
                               mcmcDRAM_model0$Chain3$bestfunp)
    
    mean_bestpar_model0<-c(mean(mcmcDRAM_model0$Chain1$bestpar["k_param"],
                                mcmcDRAM_model0$Chain2$bestpar["k_param"],
                                mcmcDRAM_model0$Chain3$bestpar["k_param"]),
                           mean(mcmcDRAM_model0$Chain1$bestpar["CUE"],
                                mcmcDRAM_model0$Chain2$bestpar["CUE"],
                                mcmcDRAM_model0$Chain3$bestpar["CUE"]),
                           mean(mcmcDRAM_model0$Chain1$bestpar["sd"],
                                mcmcDRAM_model0$Chain2$bestpar["sd"],
                                mcmcDRAM_model0$Chain3$bestpar["sd"]))
    
    #Create df with mean, sd, CI, and best
    df_bestpar_model0<-data.frame(rbind(Optim_param_model0_out,Optim_param_std_model0_out,
                                        CI95_margin_model0))
    colnames(df_bestpar_model0)<-order_pars_out_model0
    #Reorder parameters
    df_bestpar_model0<-df_bestpar_model0[c("k_param","CUE","sd")]
    #Add best (mean across chains)
    df_bestpar_model0["best",]<-mean_bestpar_model0
    
    #Add prior boundaries
    df_bestpar_model0["prior_lower",]<-lower_lim_model0
    df_bestpar_model0["prior_upper",]<-upper_lim_model0
    
    #Save df model0
    #write.csv(df_bestpar_model0,paste0(ROOT_DIR_GRAPHS,"model0_params_",diversity_variable_name,".csv"))
    
    
    ########################################################################
    #Bayesian optim model I
    ########################################################################
    
    ptm <- proc.time() #start clock to measure time
    
    print("#############################")  
    print("Start Bayes optim for model I")
    print("#############################")  
    
    #Estimation of the best-fit parameters as a starting point
    iniParPORT_modelI <- fit_modelI$par
    
    #Running 3 chains for convergence assessment
    iniParPORT_modelI <- data.frame(Chain1 = iniParPORT_modelI,
                                    Chain2 = iniParPORT_modelI,
                                    Chain3 = iniParPORT_modelI,
                                    row.names = c("k_param","gamma","CUE","sd"))
    
    iniParPORT_modelI <- sweep(iniParPORT_modelI, MARGIN = 2, STATS = c(1, 0.9, 1.1), FUN = "*")
    
    #Check limits
    for(row in 1:nrow(iniParPORT_modelI)){
      iniParPORT_modelI[row,][iniParPORT_modelI[row,] < lower_lim_modelI[row]] <- lower_lim_modelI[row]
      iniParPORT_modelI[row,][iniParPORT_modelI[row,] > upper_lim_modelI[row]] <- upper_lim_modelI[row]
    }
    
    mcmcDRAM_modelI <- apply(iniParPORT_modelI, MARGIN = 2, FUN = function(iIniParPORT_modelI) {
      FME::modMCMC(f            = J_func_modelI_MLE,
                   p            = iIniParPORT_modelI,
                   lower        = lower_lim_modelI, ## lower bounds for GR4J
                   upper        = upper_lim_modelI, ## upper bounds for GR4J
                   niter        = niter_Bayes, #at least >5000
                   jump         = NULL, #0.01
                   outputlength = outputlength_n, 
                   burninlength = burninlength_Bayes, #needs to be >0 (typical 5-10% of iterations)
                   updatecov    = updatecov_Bayes, ## adaptative Metropolis (AM)
                   ntrydr       = ntrydr_Bayes)   ## delayed rejection (RD)
    })
    
    print("Elapsed time (minutes)")
    #print(proc.time() - ptm) #end and print time elapsed (seconds)
    print((proc.time() - ptm)/60) #minutes
    
    #summary(mcmcDRAM)
    
    #MCMC diagnostics and visualisation tools
    
    multDRAM_modelI <- coda::as.mcmc.list(lapply(mcmcDRAM_modelI, FUN = function(x) {
      coda::as.mcmc(as.matrix(x$pars))
    }))
    
    gelRub_modelI <- coda::gelman.diag(multDRAM_modelI, autoburnin = TRUE)$mpsrf
    
    print("Gelman Rubin")
    print(gelRub_modelI) #if gelRub < 1.1 you can consider convergence
    
    if(gelRub_modelI>1.1){
      warning("Careful: chains have not converged to the same posterior")
    }
    
    #Convergence
    parDRAM_modelI <- ggmcmc::ggs(multDRAM_modelI) ## to convert object for using by all ggs_* graphical functions
    ggmcmc::ggs_traceplot(parDRAM_modelI)
    
    #posterior density parameters
    burnParDRAM_modelI <- parDRAM_modelI[parDRAM_modelI$Iteration > 500, ] # to keep only the second half of the series
    
    
    burnParDRAM_modelI$Parameter <- gsub("sd", "sigma", burnParDRAM_modelI$Parameter) 
    burnParDRAM_modelI$Parameter <- gsub("k_param", expression(k[0]), burnParDRAM_modelI$Parameter) 
    burnParDRAM_modelI$Parameter<-factor(burnParDRAM_modelI$Parameter,levels = c("CUE", expression(k[0]),"sigma","gamma"))
    
    
    ggmcmc::ggs_density(burnParDRAM_modelI, greek=TRUE)+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            legend.title = element_text(size=14),
            legend.text = element_text(size=14),
            strip.text = element_text(size = 14),
            panel.background = element_rect(fill = "white", color = "white"),
            panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
            panel.grid.minor = element_line(color = "white"))+
      ylab("Density")+xlab("Parameter value")
    
    #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_postPDF_modelI_withN_",diversity_variable_name,".png"))
    
    #Parameter correlation
    ggmcmc::ggs_pairs(burnParDRAM_modelI, lower = list(continuous = "density"))
    #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_corr_modelI_withN_",diversity_variable_name,".png"))
    #----------
    
    #With mean parameter values across chains
    Optim_param_modelI_out=t(aggregate(parDRAM_modelI$value, list(parDRAM_modelI$Parameter), FUN=mean))
    order_pars_out_modelI<-Optim_param_modelI_out["Group.1",]
    Optim_param_modelI_out=as.numeric(Optim_param_modelI_out["x",])
    
    Optim_param_std_modelI_out=t(aggregate(parDRAM_modelI$value, list(parDRAM_modelI$Parameter), FUN=sd))
    order_pars_out_sd_modelI<-Optim_param_std_modelI_out["Group.1",]
    Optim_param_std_modelI_out=as.numeric(Optim_param_std_modelI_out["x",])
    
    print("Mean optimum parameters - Model I")
    print(Optim_param_modelI_out)
    print("Std optimum parameters - Model I")
    print(Optim_param_std_modelI_out)
    
    CI95_margin_modelI<-(qt(0.975,df=outputlength_n-1)*Optim_param_std_modelI_out/sqrt(outputlength_n))
    CI95_lower_modelI<-Optim_param_modelI_out-CI95_margin_modelI
    CI95_upper_modelI<-Optim_param_modelI_out+CI95_margin_modelI
    
    mean_bestfunp_modelI<-mean(mcmcDRAM_modelI$Chain1$bestfunp,
                               mcmcDRAM_modelI$Chain2$bestfunp,
                               mcmcDRAM_modelI$Chain3$bestfunp)
    
    mean_bestpar_modelI<-c(mean(mcmcDRAM_modelI$Chain1$bestpar["k_param"],
                                mcmcDRAM_modelI$Chain2$bestpar["k_param"],
                                mcmcDRAM_modelI$Chain3$bestpar["k_param"]),
                           mean(mcmcDRAM_modelI$Chain1$bestpar["gamma"],
                                mcmcDRAM_modelI$Chain2$bestpar["gamma"],
                                mcmcDRAM_modelI$Chain3$bestpar["gamma"]),
                           mean(mcmcDRAM_modelI$Chain1$bestpar["CUE"],
                                mcmcDRAM_modelI$Chain2$bestpar["CUE"],
                                mcmcDRAM_modelI$Chain3$bestpar["CUE"]),
                           mean(mcmcDRAM_modelI$Chain1$bestpar["sd"],
                                mcmcDRAM_modelI$Chain2$bestpar["sd"],
                                mcmcDRAM_modelI$Chain3$bestpar["sd"]))
    
    #Create df with mean, sd, CI, and best
    df_bestpar_modelI<-data.frame(rbind(Optim_param_modelI_out,Optim_param_std_modelI_out,
                                        CI95_margin_modelI))
    colnames(df_bestpar_modelI)<-order_pars_out_modelI
    #Reorder parameters
    df_bestpar_modelI<-df_bestpar_modelI[c("k_param","gamma","CUE","sd")]
    #Add best (mean across chains)
    df_bestpar_modelI["best",]<-mean_bestpar_modelI
    
    #Add prior boundaries
    df_bestpar_modelI["prior_lower",]<-lower_lim_modelI
    df_bestpar_modelI["prior_upper",]<-upper_lim_modelI
    
    #Save df modelI
    #write.csv(df_bestpar_modelI,paste0(ROOT_DIR_GRAPHS,"modelI_params_",diversity_variable_name,".csv"))
  }
  
  ########################################################################
  #Bayesian optim model II
  ########################################################################
  
  ptm <- proc.time() #start clock to measure time
  print("#############################")  
  print("Start Bayes optim for model II")
  print(paste0("Diversity variable: ",diversity_variable_name))
  print("#############################")  
  
  #Estimation of the best-fit parameters as a starting point
  iniParPORT_modelII <- fit_modelII$par
  
  #Running 3 chains for convergence assessment
  iniParPORT_modelII <- data.frame(Chain1 = iniParPORT_modelII,
                                   Chain2 = iniParPORT_modelII,
                                   Chain3 = iniParPORT_modelII,
                                   row.names = c("k_param","gamma","beta","CUE","sd"))
  
  iniParPORT_modelII <- sweep(iniParPORT_modelII, MARGIN = 2, STATS = c(1, 0.9, 1.1), FUN = "*")
  
  #Check limits
  for(row in 1:nrow(iniParPORT_modelII)){
    iniParPORT_modelII[row,][iniParPORT_modelII[row,] < lower_lim_modelII[row]] <- lower_lim_modelII[row]
    iniParPORT_modelII[row,][iniParPORT_modelII[row,] > upper_lim_modelII[row]] <- upper_lim_modelII[row]
  }
  
  mcmcDRAM_modelII <- apply(iniParPORT_modelII, MARGIN = 2, FUN = function(iIniParPORT_modelII) {
    FME::modMCMC(f            = J_func_modelII_MLE,
                 p            = iIniParPORT_modelII,
                 lower        = lower_lim_modelII, ## lower bounds for GR4J
                 upper        = upper_lim_modelII, ## upper bounds for GR4J
                 niter        = niter_Bayes, #at least >5000
                 jump         = NULL, #0.01
                 outputlength = outputlength_n, 
                 burninlength = burninlength_Bayes, #needs to be >0 (typical 5-10% of iterations)
                 updatecov    = updatecov_Bayes, ## adaptative Metropolis (AM)
                 ntrydr       = ntrydr_Bayes)   ## delayed rejection (RD)
  })
  
  print("Elapsed time (minutes)")
  #print(proc.time() - ptm) #end and print time elapsed (seconds)
  print((proc.time() - ptm)/60) #minutes
  
  #summary(mcmcDRAM)
  
  #MCMC diagnostics and visualisation tools
  
  multDRAM_modelII <- coda::as.mcmc.list(lapply(mcmcDRAM_modelII, FUN = function(x) {
    coda::as.mcmc(as.matrix(x$pars))
  }))
  
  gelRub_modelII <- coda::gelman.diag(multDRAM_modelII, autoburnin = TRUE)$mpsrf
  
  print("Gelman Rubin")
  print(gelRub_modelII) #if gelRub < 1.1 you can consider convergence
  
  if(gelRub_modelII>1.1){
    warning("Careful: chains have not converged to the same posterior")
    print(paste0("Max prior used for beta: ",max_beta_prior))
    #stop()
  }
  
  #Convergence
  parDRAM_modelII <- ggmcmc::ggs(multDRAM_modelII) ## to convert object for using by all ggs_* graphical functions
  ggmcmc::ggs_traceplot(parDRAM_modelII)
  
  #posterior density parameters
  burnParDRAM_modelII <- parDRAM_modelII[parDRAM_modelII$Iteration > 500, ] # to keep only the second half of the series
  
  
  burnParDRAM_modelII$Parameter <- gsub("sd", "sigma", burnParDRAM_modelII$Parameter) 
  burnParDRAM_modelII$Parameter <- gsub("k_param", expression(k[0]), burnParDRAM_modelII$Parameter) 
  burnParDRAM_modelII$Parameter<-factor(burnParDRAM_modelII$Parameter,levels = c("CUE", expression(k[0]),"sigma","gamma","beta"))
  
  
  ggmcmc::ggs_density(burnParDRAM_modelII, greek=TRUE)+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          strip.text = element_text(size = 14),
          panel.background = element_rect(fill = "white", color = "white"),
          panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
          panel.grid.minor = element_line(color = "white"))+
    ylab("Density")+xlab("Parameter value")
  
  #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_postPDF_modelII_withN_",diversity_variable_name,".png"))
  
  #Parameter correlation
  ggmcmc::ggs_pairs(burnParDRAM_modelII, lower = list(continuous = "density"))
  #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_corr_modelII_withN_",diversity_variable_name,".png"))
  #----------
  
  #With mean parameter values across chains
  Optim_param_modelII_out=t(aggregate(parDRAM_modelII$value, list(parDRAM_modelII$Parameter), FUN=mean))
  order_pars_out_modelII<-Optim_param_modelII_out["Group.1",]
  Optim_param_modelII_out=as.numeric(Optim_param_modelII_out["x",])
  
  Optim_param_std_modelII_out=t(aggregate(parDRAM_modelII$value, list(parDRAM_modelII$Parameter), FUN=sd))
  order_pars_out_sd_modelII<-Optim_param_std_modelII_out["Group.1",]
  Optim_param_std_modelII_out=as.numeric(Optim_param_std_modelII_out["x",])
  
  print("Mean optimum parameters - Model II")
  print(Optim_param_modelII_out)
  print("Std optimum parameters - Model II")
  print(Optim_param_std_modelII_out)
  
  CI95_margin_modelII<-(qt(0.975,df=outputlength_n-1)*Optim_param_std_modelII_out/sqrt(outputlength_n))
  CI95_lower_modelII<-Optim_param_modelII_out-CI95_margin_modelII
  CI95_upper_modelII<-Optim_param_modelII_out+CI95_margin_modelII
  
  mean_bestfunp_modelII<-mean(mcmcDRAM_modelII$Chain1$bestfunp,
                              mcmcDRAM_modelII$Chain2$bestfunp,
                              mcmcDRAM_modelII$Chain3$bestfunp)
  
  mean_bestpar_modelII<-c(mean(mcmcDRAM_modelII$Chain1$bestpar["k_param"],
                               mcmcDRAM_modelII$Chain2$bestpar["k_param"],
                               mcmcDRAM_modelII$Chain3$bestpar["k_param"]),
                          mean(mcmcDRAM_modelII$Chain1$bestpar["gamma"],
                               mcmcDRAM_modelII$Chain2$bestpar["gamma"],
                               mcmcDRAM_modelII$Chain3$bestpar["gamma"]),
                          mean(mcmcDRAM_modelII$Chain1$bestpar["beta"],
                               mcmcDRAM_modelII$Chain2$bestpar["beta"],
                               mcmcDRAM_modelII$Chain3$bestpar["beta"]),
                          mean(mcmcDRAM_modelII$Chain1$bestpar["CUE"],
                               mcmcDRAM_modelII$Chain2$bestpar["CUE"],
                               mcmcDRAM_modelII$Chain3$bestpar["CUE"]),
                          mean(mcmcDRAM_modelII$Chain1$bestpar["sd"],
                               mcmcDRAM_modelII$Chain2$bestpar["sd"],
                               mcmcDRAM_modelII$Chain3$bestpar["sd"]))
  
  
  #Create df with mean, sd, CI, and best
  df_bestpar_modelII<-data.frame(rbind(Optim_param_modelII_out,Optim_param_std_modelII_out,
                                       CI95_margin_modelII))
  colnames(df_bestpar_modelII)<-order_pars_out_modelII
  #Reorder parameters
  df_bestpar_modelII<-df_bestpar_modelII[c("k_param","gamma","beta","CUE","sd")]
  #Add best (mean across chains)
  df_bestpar_modelII["best",]<-mean_bestpar_modelII
  
  #Add prior boundaries
  df_bestpar_modelII["prior_lower",]<-lower_lim_modelII
  df_bestpar_modelII["prior_upper",]<-upper_lim_modelII
  
  #Save df modelI
  #write.csv(df_bestpar_modelII,paste0(ROOT_DIR_GRAPHS,"modelII_params_",diversity_variable_name,".csv"))
  
  
  
  ########################################################################
  #Bayesian optim model III
  ########################################################################
  
  ptm <- proc.time() #start clock to measure time
  
  print("#############################")  
  print("Start Bayes optim for model III")
  print("#############################")  
  
  #Estimation of the best-fit parameters as a starting point
  iniParPORT_modelIII <- fit_modelIII$par
  
  #Running 3 chains for convergence assessment
  iniParPORT_modelIII <- data.frame(Chain1 = iniParPORT_modelIII,
                                    Chain2 = iniParPORT_modelIII,
                                    Chain3 = iniParPORT_modelIII,
                                    row.names = c("k_param","beta","CUE","sd"))
  
  iniParPORT_modelIII <- sweep(iniParPORT_modelIII, MARGIN = 2, STATS = c(1, 0.9, 1.1), FUN = "*")
  
  #Check limits
  for(row in 1:nrow(iniParPORT_modelIII)){
    iniParPORT_modelIII[row,][iniParPORT_modelIII[row,] < lower_lim_modelIII[row]] <- lower_lim_modelIII[row]
    iniParPORT_modelIII[row,][iniParPORT_modelIII[row,] > upper_lim_modelIII[row]] <- upper_lim_modelIII[row]
  }
  
  mcmcDRAM_modelIII <- apply(iniParPORT_modelIII, MARGIN = 2, FUN = function(iIniParPORT_modelIII) {
    FME::modMCMC(f            = J_func_modelIII_MLE,
                 p            = iIniParPORT_modelIII,
                 lower        = lower_lim_modelIII, ## lower bounds for GR4J
                 upper        = upper_lim_modelIII, ## upper bounds for GR4J
                 niter        = niter_Bayes, #at least >5000
                 jump         = NULL, #0.01
                 outputlength = outputlength_n, 
                 burninlength = burninlength_Bayes, #needs to be >0 (typical 5-10% of iterations)
                 updatecov    = updatecov_Bayes, ## adaptative Metropolis (AM)
                 ntrydr       = ntrydr_Bayes)   ## delayed rejection (RD)
  })
  
  print("Elapsed time (minutes)")
  #print(proc.time() - ptm) #end and print time elapsed (seconds)
  print((proc.time() - ptm)/60) #minutes
  
  #summary(mcmcDRAM)
  
  #MCMC diagnostics and visualisation tools
  
  multDRAM_modelIII <- coda::as.mcmc.list(lapply(mcmcDRAM_modelIII, FUN = function(x) {
    coda::as.mcmc(as.matrix(x$pars))
  }))
  
  gelRub_modelIII <- coda::gelman.diag(multDRAM_modelIII, autoburnin = TRUE)$mpsrf
  
  print("Gelman Rubin")
  print(gelRub_modelIII) #if gelRub < 1.1 you can consider convergence
  
  if(gelRub_modelIII>1.1){
    # print("Careful: chains have not converge to the same posterior")
    # stop()
    warning("Careful: chains have not converged to the same posterior")
  }
  
  #Convergence
  parDRAM_modelIII <- ggmcmc::ggs(multDRAM_modelIII) ## to convert object for using by all ggs_* graphical functions
  ggmcmc::ggs_traceplot(parDRAM_modelIII)
  
  #posterior density parameters
  burnParDRAM_modelIII <- parDRAM_modelIII[parDRAM_modelIII$Iteration > 500, ] # to keep only the second half of the series
  
  
  burnParDRAM_modelIII$Parameter <- gsub("sd", "sigma", burnParDRAM_modelIII$Parameter) 
  burnParDRAM_modelIII$Parameter <- gsub("k_param", expression(k[0]), burnParDRAM_modelIII$Parameter) 
  burnParDRAM_modelIII$Parameter<-factor(burnParDRAM_modelIII$Parameter,levels = c("CUE", expression(k[0]),"sigma","beta"))
  
  ggmcmc::ggs_density(burnParDRAM_modelIII, greek = TRUE)+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          strip.text = element_text(size = 14),
          panel.background = element_rect(fill = "white", color = "white"),
          panel.grid.major = element_line(color = "lightgray"),  # Major grid lines
          panel.grid.minor = element_line(color = "white"))+
    ylab("Density")+xlab("Parameter value")
  
  #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_postPDF_modelIII_withN_",diversity_variable_name,".png"))
  
  
  
  #Parameter correlation
  ggmcmc::ggs_pairs(burnParDRAM_modelIII, lower = list(continuous = "density"))
  #ggsave(paste0(ROOT_DIR_GRAPHS,"Bayes_corr_modelIII_withN_",diversity_variable_name,".png"))
  #----------
  
  #With mean parameter values across chains
  Optim_param_modelIII_out=t(aggregate(parDRAM_modelIII$value, list(parDRAM_modelIII$Parameter), FUN=mean))
  order_pars_out_modelIII<-Optim_param_modelIII_out["Group.1",]
  Optim_param_modelIII_out=as.numeric(Optim_param_modelIII_out["x",])
  
  Optim_param_std_modelIII_out=t(aggregate(parDRAM_modelIII$value, list(parDRAM_modelIII$Parameter), FUN=sd))
  order_pars_out_sd_modelIII<-Optim_param_std_modelIII_out["Group.1",]
  Optim_param_std_modelIII_out=as.numeric(Optim_param_std_modelIII_out["x",])
  
  print("Mean optimum parameters - Model III")
  print(Optim_param_modelIII_out)
  print("Std optimum parameters - Model III")
  print(Optim_param_std_modelIII_out)
  
  CI95_margin_modelIII<-(qt(0.975,df=outputlength_n-1)*Optim_param_std_modelIII_out/sqrt(outputlength_n))
  CI95_lower_modelIII<-Optim_param_modelIII_out-CI95_margin_modelIII
  CI95_upper_modelIII<-Optim_param_modelIII_out+CI95_margin_modelIII
  
  mean_bestfunp_modelIII<-mean(mcmcDRAM_modelIII$Chain1$bestfunp,
                               mcmcDRAM_modelIII$Chain2$bestfunp,
                               mcmcDRAM_modelIII$Chain3$bestfunp)
  
  mean_bestpar_modelIII<-c(mean(mcmcDRAM_modelIII$Chain1$bestpar["k_param"],
                                mcmcDRAM_modelIII$Chain2$bestpar["k_param"],
                                mcmcDRAM_modelIII$Chain3$bestpar["k_param"]),
                           mean(mcmcDRAM_modelIII$Chain1$bestpar["beta"],
                                mcmcDRAM_modelIII$Chain2$bestpar["beta"],
                                mcmcDRAM_modelIII$Chain3$bestpar["beta"]),
                           mean(mcmcDRAM_modelIII$Chain1$bestpar["CUE"],
                                mcmcDRAM_modelIII$Chain2$bestpar["CUE"],
                                mcmcDRAM_modelIII$Chain3$bestpar["CUE"]),
                           mean(mcmcDRAM_modelIII$Chain1$bestpar["sd"],
                                mcmcDRAM_modelIII$Chain2$bestpar["sd"],
                                mcmcDRAM_modelIII$Chain3$bestpar["sd"]))
  
  #Create df with mean, sd, CI, and best
  df_bestpar_modelIII<-data.frame(rbind(Optim_param_modelIII_out,Optim_param_std_modelIII_out,
                                        CI95_margin_modelIII))
  colnames(df_bestpar_modelIII)<-order_pars_out_modelIII
  #Reorder parameters
  df_bestpar_modelIII<-df_bestpar_modelIII[c("k_param","beta","CUE","sd")]
  #Add best (mean across chains)
  df_bestpar_modelIII["best",]<-mean_bestpar_modelIII
  
  #Add prior boundaries
  df_bestpar_modelIII["prior_lower",]<-lower_lim_modelIII
  df_bestpar_modelIII["prior_upper",]<-upper_lim_modelIII
  
  #Save df modelIII
  #write.csv(df_bestpar_modelIII,paste0(ROOT_DIR_GRAPHS,"modelIII_params_",diversity_variable_name,".csv"))
  
  
  
  
  #calculate CI beta
  beta_allchains<-subset(parDRAM_modelII,parDRAM_modelII$Parameter=="beta")$value
  dens_beta_allchains<-density(beta_allchains)
  plot(dens_beta_allchains)
  beta_quart_1st<-summary(beta_allchains)["1st Qu."]
  beta_quart_3rd<-summary(beta_allchains)["3rd Qu."]
  
  Optim_param_model0=mean_bestpar_model0
  Optim_param_modelI=mean_bestpar_modelI
  Optim_param_modelII=mean_bestpar_modelII
  Optim_param_modelIII=mean_bestpar_modelIII
  
  ##############################################
  #Simulations with optimal parameters
  ##############################################
  #Rh_general_model(k_pred,gamma_coeff,beta_coeff,diversity_variable,CUE_micr)
  #Min_general_model(Rhmodel,CUE_micr,CN_ratio_micr)
  
  if(j==1){
    ##########################
    #Model 0 (k, CUE)
    ##########################
    Rh_model_0_optim<-Rh_general_model(k_pred=Optim_param_model0[1],gamma_coeff=0,beta_coeff=0,
                                       diversity_variable=diversity_variable_in,CUE_micr=Optim_param_model0[2])
    Min_model_0_optim<-Min_general_model(Rhmodel=Rh_model_0_optim,CUE_micr=Optim_param_model0[2],CN_ratio_micr=CN_ratio_micr)
    
    #Calculate R2
    r_sq_model0 <-cor(Rh_model_0_optim,respiration_variable)^2
    r_sq_model0_expr<-bquote(R^2 == .(round(r_sq_model0, 2)))
    print(r_sq_model0_expr)
    r_sq_model0_Min <-cor(Min_model_0_optim,mineralization_variable)^2
    r_sq_model0_expr_Min<-bquote(R^2 == .(round(r_sq_model0_Min, 3)))
    print(r_sq_model0_expr_Min)
    ##########################
    
    ##########################
    #Model I (k, gamma, CUE)
    ##########################
    Rh_model_I_optim<-Rh_general_model(k_pred=Optim_param_modelI[1],gamma_coeff=Optim_param_modelI[2],beta_coeff=0,
                                       diversity_variable=diversity_variable_in,CUE_micr=Optim_param_modelI[3])
    Min_model_I_optim<-Min_general_model(Rhmodel=Rh_model_I_optim,CUE_micr=Optim_param_modelI[3],CN_ratio_micr=CN_ratio_micr)
    
    #Calculate R2
    r_sq_modelI <-cor(Rh_model_I_optim,respiration_variable)^2
    r_sq_modelI_expr<-bquote(R^2 == .(round(r_sq_modelI, 2)))
    print(r_sq_modelI_expr)
    r_sq_modelI_Min <-cor(Min_model_I_optim,mineralization_variable)^2
    r_sq_modelI_expr_Min<-bquote(R^2 == .(round(r_sq_modelI_Min, 3)))
    print(r_sq_modelI_expr_Min)  
  }
  
  
  ##########################
  #Model II (k, gamma, beta, CUE)
  ##########################
  Rh_model_II_optim<-Rh_general_model(k_pred=Optim_param_modelII[1],gamma_coeff=Optim_param_modelII[2],beta_coeff=Optim_param_modelII[3],
                                      diversity_variable=diversity_variable_in,CUE_micr=Optim_param_modelII[4])
  Min_model_II_optim<-Min_general_model(Rhmodel=Rh_model_II_optim,CUE_micr=Optim_param_modelII[4],CN_ratio_micr=CN_ratio_micr)
  
  
  #Calculate R2
  r_sq_modelII <-cor(Rh_model_II_optim,respiration_variable)^2
  r_sq_modelII_expr<-bquote(R^2 == .(round(r_sq_modelII, 2)))
  print(r_sq_modelII_expr)
  r_sq_modelII_Min <-cor(Min_model_II_optim,mineralization_variable)^2
  r_sq_modelII_expr_Min<-bquote(R^2 == .(round(r_sq_modelII_Min, 3)))
  print(r_sq_modelII_expr_Min)
  
  
  ##########################
  #Model III (k, beta, CUE)
  ##########################
  Rh_model_III_optim<-Rh_general_model(k_pred=Optim_param_modelIII[1],gamma_coeff=0,beta_coeff=Optim_param_modelIII[2],
                                       diversity_variable=diversity_variable_in,CUE_micr=Optim_param_modelIII[3])
  Min_model_III_optim<-Min_general_model(Rhmodel=Rh_model_III_optim,CUE_micr=Optim_param_modelIII[3],CN_ratio_micr=CN_ratio_micr)
  
  #Calculate R2
  r_sq_modelIII <-cor(Rh_model_III_optim,respiration_variable)^2
  r_sq_modelIII_expr<-bquote(R^2 == .(round(r_sq_modelIII, 2)))
  print(r_sq_modelIII_expr)
  r_sq_modelIII_Min <-cor(Min_model_III_optim,mineralization_variable)^2
  r_sq_modelIII_expr_Min<-bquote(R^2 == .(round(r_sq_modelIII_Min, 3)))
  print(r_sq_modelIII_expr_Min) 
  
  #Now you can calculate the BIC, AIC etc to compare models
  #BIC=−2×log(L)+k×log(n)
  
  if(j==1){
    
    BIC_model0<-2*mean_bestfunp_model0+length(lower_lim_model0)*log(length(c(respiration_variable,mineralization_variable)))
    print("BIC model 0")
    print(BIC_model0)
    
    BIC_modelI<-2*mean_bestfunp_modelI+length(lower_lim_modelI)*log(length(c(respiration_variable,mineralization_variable)))
    print("BIC model I")
    print(BIC_modelI)
  }
  
  
  
  BIC_modelII<-2*mean_bestfunp_modelII+length(lower_lim_modelII)*log(length(c(respiration_variable,mineralization_variable)))
  print("BIC model II")
  print(BIC_modelII)
  
  BIC_modelIII<-2*mean_bestfunp_modelIII+length(lower_lim_modelIII)*log(length(c(respiration_variable,mineralization_variable)))
  print("BIC model III")
  print(BIC_modelIII)
  
  #RMSE
  if(j==1){
    RMSE_model0_Rh = sqrt(sum((Rh_model_0_optim-respiration_variable)^2)/length(Rh_model_0_optim))
    RMSE_model0_Rh
    RMSE_model0_Min = sqrt(sum((Min_model_0_optim-mineralization_variable)^2)/length(Min_model_0_optim))
    RMSE_model0_Min
    
    RMSE_modelI_Rh = sqrt(sum((Rh_model_I_optim-respiration_variable)^2)/length(Rh_model_I_optim))
    RMSE_modelI_Rh
    RMSE_modelI_Min = sqrt(sum((Min_model_I_optim-mineralization_variable)^2)/length(Min_model_I_optim))
    RMSE_modelI_Min
  }
  
  RMSE_modelII_Rh = sqrt(sum((Rh_model_II_optim-respiration_variable)^2)/length(Rh_model_II_optim))
  RMSE_modelII_Rh
  RMSE_modelII_Min = sqrt(sum((Min_model_II_optim-mineralization_variable)^2)/length(Min_model_II_optim))
  RMSE_modelII_Min
  
  RMSE_modelIII_Rh = sqrt(sum((Rh_model_III_optim-respiration_variable)^2)/length(Rh_model_III_optim))
  print(paste0("RMSE Rh model III: ",RMSE_modelIII_Rh))
  RMSE_modelIII_Min = sqrt(sum((Min_model_III_optim-mineralization_variable)^2)/length(Min_model_III_optim))
  print(paste0("RMSE Min model III: ",RMSE_modelIII_Min))
  
  
  
  #Plot Model0, Model II, Model II, Model III against obs
  #only for combined richness
  if(j==1){
    pdf(file=paste0(ROOT_DIR_GRAPHS,"Figure_5.pdf"),
        width=20, height=30)
    
    #par(mfrow=c(4,2), mar=c(14,20,12,6),mgp = c(10, 3, 0))
    par(mfrow=c(4,2), mar=c(14,15,1,3),mgp = c(9, 3, 0), oma=c(10,10,5,0))
    cex_text=5.
    
    plot(respiration_variable,Rh_model_0_optim,
         ylim=c(0,0.2),xlim=c(0,0.2),
         #main="(a)",
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkgreen",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    #title(main = "(a)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(a) Model 0", adj=0, cex.main=cex_text)
    ##title(main = "Model 0",line=5,adj=1.5,cex.main=cex_text,font.main=1)
    lines(c(0,0.2),c(0,0.2),lty=2, lwd=2)
    abline(lm(Rh_model_0_optim~respiration_variable), col="darkgreen")
    
    text(0.15,0.09,r_sq_model0_expr,cex=cex_text)
    text(0.15,0.05,paste0("RMSE = ",round(RMSE_model0_Rh,2)),cex=cex_text)
    text(0.15,0.01,paste0("BIC = ",round(BIC_model0,2)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(a)", pos = 2, cex = cex_text, font = 1)
    
    plot(mineralization_variable,Min_model_0_optim,
         ylim=c(-1e-3,1e-2),xlim=c(-1e-3,1e-2),
         #main=model0_Min_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkgreen",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(b)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(b) Model 0", adj=0, cex.main=cex_text)
    lines(c(-1e-3,1e-2),c(-1e-3,1e-2),lty=2, lwd=2)
    abline(lm(Min_model_0_optim~mineralization_variable), col="darkgreen")
    
    text(6.8e-3,0.003,r_sq_model0_expr_Min,cex=cex_text)
    text(6.8e-3,0.001,paste0("RMSE = ",round(RMSE_model0_Min,3)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(b)", pos = 2, cex = cex_text, font = 1)
    
    plot(respiration_variable,Rh_model_I_optim,
         ylim=c(0,0.2),xlim=c(0,0.2),
         #main=modelI_Rh_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkblue",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(c)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(c) Model I", adj=0, cex.main=cex_text)
    lines(c(0,0.2),c(0,0.2),lty=2, lwd=2)
    abline(lm(Rh_model_I_optim~respiration_variable), col="darkblue")
    text(0.15,0.09,r_sq_modelI_expr,cex=cex_text)
    text(0.15,0.05,paste0("RMSE = ",round(RMSE_modelI_Rh,2)),cex=cex_text)
    text(0.15,0.01,paste0("BIC = ",round(BIC_modelI,2)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(c)", pos = 2, cex = cex_text, font = 1)
    
    
    plot(mineralization_variable,Min_model_I_optim,
         ylim=c(-1e-3,1e-2),xlim=c(-1e-3,1e-2),
         #main=modelI_Min_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkblue",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(d)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(d) Model I", adj=0, cex.main=cex_text)
    lines(c(-1e-3,1e-2),c(-1e-3,1e-2),lty=2, lwd=2)
    abline(lm(Min_model_I_optim~mineralization_variable), col="darkblue")
    text(6.8e-3,0.003,r_sq_modelI_expr_Min,cex=cex_text)
    text(6.8e-3,0.001,paste0("RMSE = ",round(RMSE_modelI_Min,3)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(d)", pos = 2, cex = cex_text, font = 1)
    
    
    #model III = model II in the manuscript
    plot(respiration_variable,Rh_model_III_optim,
         ylim=c(0,0.2),xlim=c(0,0.2),
         #main=modelIII_Rh_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("black",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(e)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(e) Model II", adj=0, cex.main=cex_text)
    lines(c(0,0.2),c(0,0.2),lty=2, lwd=2)
    abline(lm(Rh_model_III_optim~respiration_variable), col="black")
    text(0.15,0.09,r_sq_modelIII_expr,cex=cex_text)
    text(0.15,0.05,paste0("RMSE = ",round(RMSE_modelIII_Rh,2)),cex=cex_text)
    text(0.15,0.01,paste0("BIC = ",round(BIC_modelIII,2)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(e)", pos = 2, cex = cex_text, font = 1)
    
    
    plot(mineralization_variable,Min_model_III_optim,
         ylim=c(-1e-3,1e-2),xlim=c(-1e-3,1e-2),
         #main=modelIII_Min_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("black",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(f)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(f) Model II", adj=0, cex.main=cex_text)
    lines(c(-1e-3,1e-2),c(-1e-3,1e-2),lty=2, lwd=2)
    abline(lm(Min_model_III_optim~mineralization_variable), col="black")
    text(6.8e-3,0.003,r_sq_modelIII_expr_Min,cex=cex_text)
    text(6.8e-3,0.001,paste0("RMSE = ",round(RMSE_modelIII_Min,3)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(f)", pos = 2, cex = cex_text, font = 1)
    
    
    #model II = model III in the ms
    plot(respiration_variable,Rh_model_II_optim,
         ylim=c(0,0.2),xlim=c(0,0.2),
         #main=modelII_Rh_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkred",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(g)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(g) Model III", adj=0, cex.main=cex_text)
    lines(c(0,0.2),c(0,0.2),lty=2, lwd=2)
    abline(lm(Rh_model_II_optim~respiration_variable), col="darkred")
    text(0.15,0.09,r_sq_modelII_expr,cex=cex_text)
    text(0.15,0.05,paste0("RMSE = ",round(RMSE_modelII_Rh,2)),cex=cex_text)
    text(0.15,0.01,paste0("BIC = ",round(BIC_modelII,2)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(g)", pos = 2, cex = cex_text, font = 1)
    
    plot(mineralization_variable,Min_model_II_optim,
         ylim=c(-1e-3,1e-2),xlim=c(-1e-3,1e-2),
         #main=modelII_Min_expr,
         ylab=expression("Predicted"),
         xlab=expression("Observed"),
         pch=19,col=grDevices::adjustcolor("darkred",alpha=0.5),
         cex=5,cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text)
    ##title(main = "(h)", adj=0, line=1.5,cex.main=cex_text,font.main=1)
    #title(main="(h) Model III", adj=0, cex.main=cex_text)
    lines(c(-1e-3,1e-2),c(-1e-3,1e-2),lty=2, lwd=2)
    abline(lm(Min_model_II_optim~mineralization_variable), col="darkred")
    text(6.8e-3,0.003,r_sq_modelII_expr_Min,cex=cex_text)
    text(6.8e-3,0.001,paste0("RMSE = ",round(RMSE_modelII_Min,3)),cex=cex_text)
    text(x = par("usr")[2] * 0.1, y = par("usr")[4] * 0.9, 
         "(h)", pos = 2, cex = cex_text, font = 1)
    
    
    mtext("Model 0", outer = TRUE, side = 2, line = 2, cex = 4,font = 2, adj=0.95)
    mtext("Model I", outer = TRUE, side = 2, line = 2, cex = 4,font = 2, adj=0.67)
    mtext("Model II", outer = TRUE, side = 2, line = 2, cex = 4,font = 2, adj=0.39)
    mtext("Model III", outer = TRUE, side = 2, line = 2, cex = 4,font = 2, adj=0.12)
    
    
    mtext(bquote(bold(" R (g C"~kg^{-1}~d^{-1}~")")), side = 1, outer = TRUE, line = 4, cex = 4, adj=0.2, font=2)
    mtext(bquote(bold(" M (g N"~kg^{-1}~d^{-1}~")")), side = 1, outer = TRUE, line = 4, cex = 4, adj=0.9, font=2)
    dev.off()
    
    
    #Plot diversity function
    
    pdf(file=paste0(ROOT_DIR_GRAPHS,"Figure_S3.pdf"),
        width=5, height=5)
    
    #par(mfrow=c(1,1), mar=c(10,12,2,2),mgp = c(8, 3, 0))
    par(mfrow=c(1,1), mar=c(7,8,3,2),mgp = c(5, 2, 0))
    cex_text = 2.
    
    sorted_diversity_variable_in<-sort(diversity_variable_in)
    
    
    div_func_opt_unlimited<-Optim_param_modelII[1]*diversity_function(Optim_param_modelII[3],sorted_diversity_variable_in)
    div_func_opt<-div_func_opt_unlimited
    div_func_opt[div_func_opt > 1] <-1
    
    div_func_1stq_unlimited<-Optim_param_modelII[1]*diversity_function(as.numeric(beta_quart_1st),sorted_diversity_variable_in)
    div_func_1stq<-div_func_1stq_unlimited
    div_func_1stq[div_func_1stq > 1] <-1
    
    div_func_3rdq_unlimited<-Optim_param_modelII[1]*diversity_function(as.numeric(beta_quart_3rd),sorted_diversity_variable_in)
    div_func_3rdq<-div_func_3rdq_unlimited
    div_func_3rdq[div_func_3rdq > 1] <-1
    
    if(diversity_variable_name=="sqrt(Richness_fungi*Richness_bact)"){
      expr_div_function<-expression(f[0]==sqrt(S[fungi]*S[bacteria]))
    }else{
      expr_div_function<-diversity_variable_name
    }
    
    plot(sorted_diversity_variable_in,div_func_opt,
         ylim=c(0,0.001),
         ylab=expression("Optimized"~k[U]),
         xlab=expr_div_function,
         pch=19,col="black",
         cex.axis=cex_text,cex.lab=cex_text,cex.main=cex_text,type="b")
    points(sorted_diversity_variable_in,div_func_1stq,pch=17,col="grey")
    lines(sorted_diversity_variable_in,div_func_1stq,col="grey")
    points(sorted_diversity_variable_in,div_func_3rdq,pch=15,col="grey")
    lines(sorted_diversity_variable_in,div_func_3rdq,col="grey")
    
    # text(500,0.35,fdiv_func_expr,cex=cex_text)
    # text(500,0.2,div_variable_expr,cex=cex_text)
    #text(450,0.4,"Diversity index: Richness",cex=cex_text)
    
    # Add legend
    legend("topright", legend = c("Mean", "1stQ", "3rdQ"),
           col = c("black", "grey", "grey"), lty = 1, pch = c(19,17,15), box.lwd = 0.1,
           title = expression(beta), cex = 1.5)
    
    dev.off()
  }

  
  
  if(j==1){
    #Make plots for all models
    cc0<-ggs_crosscorrelation(burnParDRAM_model0, absolute_scale=FALSE, greek=TRUE)+
      ggtitle("(a) Model 0")
    cc1<-ggs_crosscorrelation(burnParDRAM_modelI, absolute_scale=FALSE, greek=TRUE)+
      ggtitle("(b) Model I")
    cc2<-ggs_crosscorrelation(burnParDRAM_modelIII, absolute_scale=FALSE, greek=TRUE)+ #mind the change in model names (II = III in ms)
      ggtitle("(c) Model II") 
    cc3<-ggs_crosscorrelation(burnParDRAM_modelII, absolute_scale=FALSE, greek=TRUE)+
      ggtitle("(d) Model III")
    
    # Combine the plots into a 2x2 grid 
    combined_plot_cc <- plot_grid(cc0,cc1,
                                  cc2,cc3,
                                  ncol = 2, align = "hv")
    print(combined_plot_cc)
    ggsave(plot=combined_plot_cc,paste0(ROOT_DIR_GRAPHS,"Figure_S5.pdf"),
           width = 5, height = 4)

    cp0<-ggs_caterpillar(burnParDRAM_model0, greek=TRUE)+
      ggtitle("(a) Model 0")
    cp1<-ggs_caterpillar(burnParDRAM_modelI, greek=TRUE)+
      ggtitle("(b) Model I")
    cp2<-ggs_caterpillar(burnParDRAM_modelIII, greek=TRUE)+ #mind the change in model names (II = III in ms)
      ggtitle("(c) Model II") 
    cp3<-ggs_caterpillar(burnParDRAM_modelII, greek=TRUE)+
      ggtitle("(d) Model III")
    
    # Combine the plots into a 2x2 grid 
    combined_plot_cp <- plot_grid(cp0,cp1,
                                  cp2,cp3,
                                  ncol = 2, align = "hv")
    print(combined_plot_cp)
    # ggsave(plot=combined_plot_cp,paste0(ROOT_DIR,"NEW_GRAPHS/NEW/","caterpillar_models_params.png"),
    #        width = 5, height = 4)
    
    
    #Calculate correlation matrixes
    
    # Convert parDRAM_model0 to wide format
    mcmc_wide0 <- burnParDRAM_model0 %>%
      pivot_wider(names_from = Parameter, values_from = value)
    
    # Remove columns that are not parameters (like Iteration or Chain, if present)
    mcmc_wide0 <- mcmc_wide0 %>% select(-Iteration, -Chain)
    
    # Calculate the cross-correlation matrix
    cross_cor_matrix0 <- cor(mcmc_wide0, use = "pairwise.complete.obs")
    
    
    # Convert parDRAM_modelI to wide format
    mcmc_wideI <- burnParDRAM_modelI %>%
      pivot_wider(names_from = Parameter, values_from = value)
    
    # Remove columns that are not parameters (like Iteration or Chain, if present)
    mcmc_wideI <- mcmc_wideI %>% select(-Iteration, -Chain)
    
    # Calculate the cross-correlation matrix
    cross_cor_matrixI <- cor(mcmc_wideI, use = "pairwise.complete.obs")    
    
    
    # Convert parDRAM_modelII to wide format
    mcmc_wideII <- burnParDRAM_modelIII %>%
      pivot_wider(names_from = Parameter, values_from = value)
    
    # Remove columns that are not parameters (like Iteration or Chain, if present)
    mcmc_wideII <- mcmc_wideII %>% select(-Iteration, -Chain)
    
    # Calculate the cross-correlation matrix
    cross_cor_matrixII <- cor(mcmc_wideII, use = "pairwise.complete.obs")
    
    
    # Convert parDRAM_modelIII to wide format
    mcmc_wideIII <- burnParDRAM_modelII %>%
      pivot_wider(names_from = Parameter, values_from = value)
    
    # Remove columns that are not parameters (like Iteration or Chain, if present)
    mcmc_wideIII <- mcmc_wideIII %>% select(-Iteration, -Chain)
    
    # Calculate the cross-correlation matrix
    cross_cor_matrixIII <- cor(mcmc_wideIII, use = "pairwise.complete.obs")
    
    cross_cor_matrix0
    cross_cor_matrixI
    cross_cor_matrixII
    cross_cor_matrixIII
    
    
    abs(min(cross_cor_matrix0))
    max(cross_cor_matrix0[row(cross_cor_matrix0)!=col(cross_cor_matrix0)])
    
    
    abs(min(cross_cor_matrixI))
    max(cross_cor_matrixI[row(cross_cor_matrixI)!=col(cross_cor_matrixI)])
    
    
    abs(min(cross_cor_matrixII))
    max(cross_cor_matrixII[row(cross_cor_matrixII)!=col(cross_cor_matrixII)])
    
    
    abs(min(cross_cor_matrixIII))
    max(cross_cor_matrixIII[row(cross_cor_matrixIII)!=col(cross_cor_matrixIII)])

    #Posterior distributions of the parameters - all models together
    cairo_pdf(file=paste0(ROOT_DIR_GRAPHS,"Figure_6.pdf"),
        width=30, height=30)

    #Plot the distribution density of the parameters
    par(mfrow=c(4,4), mar=c(1.4, 1.4, 0.7, 0.2), mgp = c(1.7, 0.4, 0),oma = c(2, 3, 0, 0),cex=5)
    
    #par(mfrow=c(4,4), mar=c(2.8, 2.8, 1.1, 0.2), mgp = c(1.7, 0.6, 0),cex=5)
    #Model 0
    plot(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==1)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    #title(main="(a) Model 0", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(a)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==2)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==3)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="k[0]" & burnParDRAM_model0$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_model0[1], col = "black", lty = 2, lwd = 2)
    #par(mar=c(3, 2, 1, 1))
    
    plot(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==1)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    #title(main="(b) Model 0", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(b)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_model0,burnParDRAM_model0$Parameter=="CUE" & burnParDRAM_model0$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_model0[2], col = "black", lty = 2, lwd = 2)
    
    plot.new()
    plot.new()
    
    
    #par(mar=c(3, 4, 1, 1))
    #Model I
    plot(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==1)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    #title(main="(c) Model I", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(c)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==2)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==3)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="k[0]" & burnParDRAM_modelI$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelI[1], col = "black", lty = 2, lwd = 2)
    
    #par(mar=c(3, 2, 1, 1))
    
    plot(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    #title(main="(d) Model I", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(d)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==3)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==1)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="CUE" & burnParDRAM_modelI$Chain==1)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelI[3], col = "black", lty = 2, lwd = 2)
    
    plot(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==1)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0.6,1))
    #title(main="(e) Model I", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.7, y = par("usr")[4] * 0.9, 
         "(e)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0.6,1))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0.6,1))
    polygon(density(subset(burnParDRAM_modelI,burnParDRAM_modelI$Parameter=="gamma" & burnParDRAM_modelI$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelI[2], col = "black", lty = 2, lwd = 2)
    
    plot.new()
    
    
    
    #par(mar=c(3, 4, 1, 1))
    #Model II
    plot(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==3)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    #title(main="(f) Model II", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(f)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==3)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==2)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==1)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="k[0]" & burnParDRAM_modelIII$Chain==1)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelIII[1], col = "black", lty = 2, lwd = 2)
    #par(mar=c(3, 2, 1, 1))
    
    plot(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    #title(main="(g) Model II", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(g)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==3)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==1)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="CUE" & burnParDRAM_modelIII$Chain==1)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelIII[3], col = "black", lty = 2, lwd = 2)
    
    plot.new()
    
    plot(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(-0.005,0.001))
    #title(main="(h) Model II", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(h)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==3)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(-0.005,0.001))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==1)$value),main=" ", xlab=" ", ylab=" ", xlim=c(-0.005,0.001))
    polygon(density(subset(burnParDRAM_modelIII,burnParDRAM_modelIII$Parameter=="beta" & burnParDRAM_modelIII$Chain==1)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelIII[2], col = "black", lty = 2, lwd = 2)
    
    
    #par(mar=c(3, 4, 1, 1))
    #Model III
    plot(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==1)$value),
         main=" ",  xlim=c(0,0.005),xlab =bquote(bold(k[0])))
    #title(main="(i) Model III", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(i)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==2)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==3)$value),main=" ", xlab=" ", xlim=c(0,0.005))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="k[0]" & burnParDRAM_modelII$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelII[1], col = "black", lty = 2, lwd = 2)
    #par(mar=c(3, 2, 1, 1))
    
    plot(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==1)$value),ylab=" ",
         main=" ", xlim=c(0,0.8),xlab = "CUE", font.lab = 2)
    #title(main="(j) Model III", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(j)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0,0.8))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="CUE" & burnParDRAM_modelII$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelII[4], col = "black", lty = 2, lwd = 2)
    
    plot(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==1)$value),ylab=" ",
         main=" ", xlim=c(0.6,1))
    #title(main="(k) Model III", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.7, y = par("usr")[4] * 0.9, 
         "(k)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0.6,1))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(0.6,1))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="gamma" & burnParDRAM_modelII$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelII[2], col = "black", lty = 2, lwd = 2)
    
    plot(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==1)$value),ylab=" ",
         main=" ", xlim=c(-0.005,0.001))
    #title(main="(l) Model III", font.main=1, cex.main=1)
    text(x = par("usr")[2] * 0.95, y = par("usr")[4] * 0.9, 
         "(l)", pos = 2, cex = 1, font = 1)
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==1)$value),border="#506fb8",col=grDevices::adjustcolor("#506fb8",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==2)$value),main=" ", xlab=" ", ylab=" ", xlim=c(-0.005,0.001))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==2)$value),border="#844ea6",col=grDevices::adjustcolor("#844ea6",alpha=0.5))
    
    lines(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==3)$value),main=" ", xlab=" ", ylab=" ", xlim=c(-0.005,0.001))
    polygon(density(subset(burnParDRAM_modelII,burnParDRAM_modelII$Parameter=="beta" & burnParDRAM_modelII$Chain==3)$value),border="lightpink",col=grDevices::adjustcolor("lightpink",alpha=0.5))
    
    # Add a vertical line corresponding to mean best param (i.e., mean across chain of param that achieves highest probability density)
    abline(v = mean_bestpar_modelII[3], col = "black", lty = 2, lwd = 2)
    
    
    # Add x and y labels in the outer margin
    mtext(bquote(bold(k[0])), side = 1, outer = TRUE, line = 0.4, cex = cex_text*2.5, adj=0.12)
    mtext("CUE", side = 1, outer = TRUE, line = 0.4, cex = cex_text*2.5, adj=0.39, font=2)
    mtext(bquote(bold("\u03b3")), side = 1, outer = TRUE, line = 0.4, cex = cex_text*2.5, adj=0.65)
    mtext(bquote(bold("\u03b2")), side = 1, outer = TRUE, line = 0.4, cex = cex_text*2.5, adj=0.9)

    mtext("Density", side = 2, outer = TRUE, line = 2, cex = cex_text*2.5)
    
    mtext("Model 0", font=2, side = 2, outer = TRUE, line = 0.5, cex = cex_text*2.5, adj=0.92)
    mtext("Model I", font=2, side = 2, outer = TRUE, line = 0.5, cex = cex_text*2.5, adj=0.65)
    mtext("Model II", font=2, side = 2, outer = TRUE, line = 0.5, cex = cex_text*2.5, adj=0.37)
    mtext("Model III", font=2, side = 2, outer = TRUE, line = 0.5, cex = cex_text*2.5, adj=0.1)
    
    dev.off()
    
    
    #Save model predictions to csv file
    
    db_for_model_incremental_expl <- cbind("SOC"=data_mic$C_Org_cont,
                                           micr_biom_content_variable,
                                           diversity_variable_in,
                                           respiration_variable,
                                           Rh_model_0_optim,Rh_model_I_optim,Rh_model_III_optim,Rh_model_II_optim,
                                           mineralization_variable,
                                           Min_model_0_optim,Min_model_I_optim,Min_model_III_optim,Min_model_II_optim)
    
    #write.csv(db_for_model_incremental_expl,paste0(ROOT_DIR_GRAPHS,"predicted_observed_Rh_Min.csv"))
    
  }
  
  
  
  
  
  #Create dataframe with  model performances across all diversity indexes
  
  if(j==1){
    dataframe_modeloutputs<-data.frame(cbind(c(BIC_model0, RMSE_model0_Rh,RMSE_model0_Min,r_sq_model0,r_sq_model0_Min),
                                             c(BIC_modelI, RMSE_modelI_Rh,RMSE_modelI_Min, r_sq_modelI,r_sq_modelI_Min),
                                             c(BIC_modelII, RMSE_modelII_Rh,RMSE_modelII_Min, r_sq_modelII,r_sq_modelII_Min),
                                             c(BIC_modelIII, RMSE_modelIII_Rh,RMSE_modelIII_Min, r_sq_modelIII,r_sq_modelIII_Min)))
    colnames(dataframe_modeloutputs)<-c("Model0","ModelI",paste0("ModelII_",diversity_variable_name),paste0("ModelIII_",diversity_variable_name))
    rownames(dataframe_modeloutputs)<-c("BIC","RMSE_Rh","RMSE_Min","R2_Rh","R2_Min")
    
  }else{
    
    dataframe_modeloutputs$new_diversity<-c(BIC_modelII, RMSE_modelII_Rh,RMSE_modelII_Min, r_sq_modelII,r_sq_modelII_Min)
    colnames(dataframe_modeloutputs)[colnames(dataframe_modeloutputs)=="new_diversity"]<-paste0("ModelII_",diversity_variable_name)
    dataframe_modeloutputs$new_diversity_modelIII<-c(BIC_modelIII, RMSE_modelIII_Rh,RMSE_modelIII_Min, r_sq_modelIII,r_sq_modelIII_Min)
    colnames(dataframe_modeloutputs)[colnames(dataframe_modeloutputs)=="new_diversity_only"]<-paste0("ModelIII_",diversity_variable_name)
  }
  
  
  #stop()
  
  j<-j+1
  
}

write.csv(round(dataframe_modeloutputs,3),paste0(ROOT_DIR_GRAPHS,"Table_3_and_S1.csv"))
