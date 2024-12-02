#####################
#Figures in manuscript
#CN model
#####################

#---------------
#Upload libraries
#---------------
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

library(ggspatial)
library(sf)
#---------------
#Set working directories
#---------------
ROOT_DIR<- "/Users/ebruni/Desktop/HOLISOILS/DATI/MICROBES/"
ROOT_DIR_GRAPHS<- "/Users/ebruni/Desktop/HOLISOILS/MANUSCRIPTS/MICROB_DIVERSITY/FIGURES/"

#Upload data
data_mic_all<- read_csv(paste0(ROOT_DIR,"base_datos_metro_186&diversity.csv"))

#Upload function to calculate R2 in lme model
source(paste0(ROOT_DIR,"calculate_R2_lme.R")) 

#Convert Rh_C data to numeric
data_mic_all$Rh_C_corrected_day<-as.numeric(data_mic_all$Rh_C_corrected_day)

###################
# Explore the data
###################
#Rename variables for convenience
data_mic_all$Bacteria_Pielou_evenness<-data_mic_all$Eveness_bact
data_mic_all$Bacteria_Shannon <-data_mic_all$Shannon_bact
data_mic_all$Bacteria_invsimpson <-data_mic_all$Inv_Simpson_bact

data_mic_all$Fungi_Pielou_evenness<-data_mic_all$Eveness_fungi
data_mic_all$Fungi_Shannon <-data_mic_all$Shannon_fungi
data_mic_all$Fungi_invsimpson <-data_mic_all$Inv_Simpson_fungi

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

data_mic_all$C_Org_conc<- data_mic_all$C_Org_kg_m/data_mic_all$Bulk_Density #%
#C content (g/kg)
data_mic_all$C_Org_cont<- data_mic_all$C_Org_conc *10

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


respiration_variable_expr <-expression("Heterotrophic respiration (g C"~kg^{-1}~day^{-1}~")")
respiration_variable_expr_log<-expression("Heterotrophic respiration (g C"~kg^{-1}~day^{-1}~"),"~log[10])

mineralization_variable_expr <-expression("Net N mineralization (g N"~kg^{-1}~day^{-1}~")")

carbon_content_variable_expr <-expression("SOC content (g C"~kg^{-1}~")")
carbon_content_variable_expr_log <-expression("SOC content (g C"~kg^{-1}~"),"~log[10])

carbon_stock_variable_expr_log <-expression("SOC stock (kg C"~m^{-2}~"),"~log[10])

micr_biom_content_variable_expr <-expression("Microbial biomass (g C"~kg^{-1}~")")
micr_biom_content_variable_expr_log <-expression("Microbial biomass (g C"~kg^{-1}~"),"~log[10])

richness_bactfungi_expr<- expression("Microbial diversity ("~sqrt(S[fungi]~"*"~S[bacteria])~")")
richness_bactfungi_expr_log<- expression("Microbial diversity ("~sqrt(S[fungi]~"*"~S[bacteria])~"),"~log[10])

##############################################


###########
#FIGURE 1
##########

# Get world map data
world_map <- map_data("world")

# Filter data for Spain
spainmap <- subset(world_map, region %in% "Spain")

# Assign CRS if missing (e.g., WGS84)
spainmap_sf <- sf::st_as_sf(spainmap, coords = c("long", "lat"), crs = 4326)

# Convert to sf object with projection WGS 84
spainmap_sf <- sf::st_as_sf(spainmap, coords = c("long", "lat"), crs = 4326)

# Separate mainland and islands based on the 'subregion' column
mainland_spain <- spainmap_sf %>% 
  filter(is.na(subregion))  # Mainland (NA in subregion)

# # Create a list of islands based on subregion names
# island_names <- c("Ibiza", "Majorca", "Minorca", "Formentera")  # Add more island names as needed

# islands_spain <- spainmap_sf %>% 
#   filter(subregion %in% island_names)  # Filter for islands

# Convert mainland into polygons
mainland_polygon <- mainland_spain %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

#Convert names to uniform with paper Garcia Angulo and add accents
data_mic$Site_formap<-ifelse(data_mic$Site=="Alicante","Alcoy",
                             ifelse(data_mic$Site=="Madrid","Tres Cantos",
                                    ifelse(data_mic$Site=="Chapineria","Chapinería",
                                           ifelse(data_mic$Site=="Talavera","Talavera de la Reina",  
                                                  ifelse(data_mic$Site=="Navarra","Pamplona",
                                                         ifelse(data_mic$Site=="Teruel","Formiche Bajo",
                                                                ifelse(data_mic$Site=="Almeria","Almería",
                                                                       ifelse(data_mic$Site=="Huesca","Arascués",
                                                                              ifelse(data_mic$Site=="Sev-Usera","Sevilla-Usera",
                                                                                     ifelse(data_mic$Site=="Leon","León",data_mic$Site
                                                                                            
                                                                                     ))))))))))


pdf(file=paste0(ROOT_DIR_GRAPHS,"Figure_1.pdf"),
    width=10, height=6)

par(mar=c(3,4,3,2),oma=c(5,0,1,4))
ggplot() +
  #geom_sf(data = spainmap_sf, fill = "lightgray", col = "black", size = 0.5)+
  # Plot mainland and islands separately
  geom_sf(data = mainland_polygon, color = "black", fill = "lightgray") +  # Mainland
  #geom_sf(data = island_polygons, color = "black", fill = "lightgray") +  # Islands
  expand_limits(x = spainmap$long, y = spainmap$lat)+
  geom_point(aes(x=Longitude, y=Latitude,color=Manejo_Conteo, shape=Manejo_Conteo), 
             data = data_mic,size = 5)+
  geom_text_repel(aes(x=Longitude, y=Latitude,label=Site_formap),
                  data=subset(data_mic,data_mic$Type=="Healthy"),
                  max.overlaps=Inf, colour="black", segment.colour="black",force=2,
                  box.padding = 1.2,point.padding =1.2, show.legend = FALSE,size=6)+
  guides(color=guide_legend(title="Land-use"),shape=guide_legend(title="Land-use"))+
  labs(x = "Longitude", y=" Latitude",size=5)+
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+ 
  scale_shape_manual(values = c("Forest" = 16,  # Solid circle
                                "Dehesa" = 17,  # Triangle
                                "Open woodland" = 15)) +  # Solid square
  theme(legend.text=element_text(size=22),legend.title=element_text(size=22),
        axis.text = element_text(size = 22),axis.title = element_text(size = 22),
        plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "whitesmoke"),)+
  annotation_north_arrow(location = "br", which_north = "true") +  # Compass rose
  coord_sf(datum = st_crs(spainmap_sf))  # Add this to handle geographic coordinates
dev.off()


######################
#Graphs for manuscript
######################
p1<-ggplot(data_mic,aes(x=micr_biom_content_variable,y=respiration_variable)) + 
  #scale_x_continuous(trans = 'log10') +
  geom_point(aes(x=micr_biom_content_variable,y=respiration_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  geom_smooth(method='lm',formula ='y~x' ,color="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8)+
  #facet_wrap(~ Manejo_Conteo, scales = "free")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        plot.title= element_text(size=22))+
  labs(color = "Land-use", shape = "Land-use")+
  xlab(micr_biom_content_variable_expr) + ylab(respiration_variable_expr)+ 
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(a)")

p2<-ggplot(data_mic,aes(x=micr_biom_content_variable,y=mineralization_variable)) + 
  #scale_x_continuous(trans = 'log10') +
  geom_point(aes(x=micr_biom_content_variable,y=mineralization_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  geom_smooth(method='lm',formula ='y~x' ,color="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8)+
  #facet_wrap(~ Manejo_Conteo, scales = "free")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  #labs(color = "Land-use")+
  xlab(micr_biom_content_variable_expr) + ylab(mineralization_variable_expr)+ 
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(b)")


# # Respiration vs Richness
p3<-ggplot(data_mic,aes(x=sqrt(Richness_bact*Richness_fungi),y=respiration_variable)) +
  geom_point(aes(x=sqrt(Richness_bact*Richness_fungi),y=respiration_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  # scale_x_continuous(trans = 'log10') +
  # scale_y_continuous(trans = 'log10')+
  geom_smooth(method='lm',formula ='y~x' ,color="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8)+
  xlab(richness_bactfungi_expr) + ylab(respiration_variable_expr)+
  #guides(color=guide_legend(title="Land-use"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(c)")

# # Mineralization vs Richness
p4<-ggplot(data_mic,aes(x=sqrt(Richness_bact*Richness_fungi),y=mineralization_variable)) +
  geom_point(aes(x=sqrt(Richness_bact*Richness_fungi),y=mineralization_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  # scale_x_continuous(trans = 'log10') +
  # scale_y_continuous(trans = 'log10')+
  geom_smooth(method='lm',formula ='y~x' ,color="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8)+
  xlab(richness_bactfungi_expr) + ylab(mineralization_variable_expr)+
  #guides(color=guide_legend(title="Land-use"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(d)")


#summary(lm(respiration_variable~sqrt(data_mic$Richness_bact*data_mic$Richness_fungi)))

# # Combine the plots in a 2x2 layout
# grid.arrange(p1, p2, p3, p4, ncol = 2)

# Extract legend from one of the plots
legend_plots <- get_legend(p1)

# Combine the plots into a 2x2 grid and add the shared legend
combined_plot <- plot_grid(p1 + theme(legend.position = "none"), 
                           p2, p3, p4, 
                           ncol = 2, align = "hv")

# Add the legend below the combined plot
final_plot <- plot_grid(combined_plot, legend_plots, ncol = 2, rel_widths = c(5, 1))
print(final_plot)
ggsave(plot=final_plot,paste0(ROOT_DIR_GRAPHS,"Figure_4.pdf"),
       width = 14, height = 12)



###


p1<-ggplot(data_mic,aes(x=micr_biom_content_variable,y=respiration_variable)) + 
  #scale_x_continuous(trans = 'log10') +
  geom_point(aes(x=micr_biom_content_variable,y=respiration_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  #geom_smooth(method='lm',formula ='y~x' ,color="black")+
  geom_smooth(aes(color=Manejo_Conteo), method='lm',formula ='y~x')+
  stat_cor(aes(color=Manejo_Conteo,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8, show.legend = FALSE)+
  #facet_wrap(~ Manejo_Conteo, scales = "free")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        plot.title= element_text(size=22))+
  labs(color = "Land-use", shape = "Land-use")+
  xlab(micr_biom_content_variable_expr) + ylab(respiration_variable_expr)+ 
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(a)")

p2<-ggplot(data_mic,aes(x=micr_biom_content_variable,y=mineralization_variable)) + 
  #scale_x_continuous(trans = 'log10') +
  geom_point(aes(x=micr_biom_content_variable,y=mineralization_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  #geom_smooth(method='lm',formula ='y~x' ,color="black")+
  geom_smooth(aes(color=Manejo_Conteo), method='lm',formula ='y~x')+
  stat_cor(aes(color=Manejo_Conteo,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8, show.legend = FALSE)+
  #facet_wrap(~ Manejo_Conteo, scales = "free")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  #labs(color = "Land-use")+
  xlab(micr_biom_content_variable_expr) + ylab(mineralization_variable_expr)+ 
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(b)")


# # Respiration vs Richness
p3<-ggplot(data_mic,aes(x=sqrt(Richness_bact*Richness_fungi),y=respiration_variable)) +
  geom_point(aes(x=sqrt(Richness_bact*Richness_fungi),y=respiration_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  # scale_x_continuous(trans = 'log10') +
  # scale_y_continuous(trans = 'log10')+
  #geom_smooth(method='lm',formula ='y~x' ,color="black")+
  geom_smooth(aes(color=Manejo_Conteo), method='lm',formula ='y~x')+
  stat_cor(aes(color=Manejo_Conteo,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8, show.legend = FALSE)+
  xlab(richness_bactfungi_expr) + ylab(respiration_variable_expr)+
  #guides(color=guide_legend(title="Land-use"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(c)")

# # Mineralization vs Richness
p4<-ggplot(data_mic,aes(x=sqrt(Richness_bact*Richness_fungi),y=mineralization_variable)) +
  geom_point(aes(x=sqrt(Richness_bact*Richness_fungi),y=mineralization_variable,color=Manejo_Conteo,shape=Manejo_Conteo),size=5, stroke = 2)+
  # scale_x_continuous(trans = 'log10') +
  # scale_y_continuous(trans = 'log10')+
  #geom_smooth(method='lm',formula ='y~x' ,color="black")+
  geom_smooth(aes(color=Manejo_Conteo), method='lm',formula ='y~x')+
  stat_cor(aes(color=Manejo_Conteo,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size=8, show.legend = FALSE)+
  xlab(richness_bactfungi_expr) + ylab(mineralization_variable_expr)+
  #guides(color=guide_legend(title="Land-use"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),  # Major grid lines
        panel.grid.minor = element_line(color = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",   # Remove legend from individual plot
        plot.title= element_text(size=22))+
  scale_color_manual(values=c("Forest"= "#517664","Dehesa"="#F9B5AC",
                              "Open woodland"="#E69F00"))+
  scale_shape_manual(values = c("Forest" = 16, "Dehesa" = 17, "Open woodland" = 3))+  # Circle, triangle, plus
  ggtitle("(d)")


#summary(lm(respiration_variable~sqrt(data_mic$Richness_bact*data_mic$Richness_fungi)))

# # Combine the plots in a 2x2 layout
# grid.arrange(p1, p2, p3, p4, ncol = 2)

# Extract legend from one of the plots
legend_plots <- get_legend(p1)

# Combine the plots into a 2x2 grid and add the shared legend
combined_plot <- plot_grid(p1 + theme(legend.position = "none"), 
                           p2, p3, p4, 
                           ncol = 2, align = "hv")

# Add the legend below the combined plot
final_plot <- plot_grid(combined_plot, legend_plots, ncol = 2, rel_widths = c(5, 1))
print(final_plot)
ggsave(plot=final_plot,paste0(ROOT_DIR_GRAPHS,"Figure_S2.pdf"),
       width = 14, height = 12)

#################
# FIGURE 4, 5 6 and supp in Spanish_database_Jorge_model_comparison_WITH_N_alldiv_withModel3.R
################
