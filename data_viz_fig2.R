#Loading packages
library(tidyverse)
library(stringr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(ggspatial)
library(rfishbase)
library(rredlist)
library(ggrepel)
library(fmsb)
library(scales)

#Setting directory
setwd("/Users/YourPathFile")

#Setting metadata
df<- read.csv("primers_specific_unique-cleaned.csv")

#Filtering only rapid tools primers
df<- df %>%
  filter(purpose == "rapid_tools") ##228 primers

df_sr<- df %>%
  filter(!Shark_Ray_Chimaera == "na") %>% 
  separate_rows(Shark_Ray_Chimaera, sep = "&")%>%
  group_by(Shark_Ray_Chimaera)%>%
  summarise(count=n()) #193 sharks and 35 rays

df_spp<- df %>%
  separate_rows(Species, sep = "; ")%>%
  filter(!Species == "na")%>%
  group_by(Species, Shark_Ray_Chimaera)%>%
  summarise(count=n())%>%
  arrange(desc(count))#58 species - 47 sharks and 11 rays 

df_genus<- df %>%
  filter(!Genus == "na")%>%
  group_by(Genus)%>%
  summarise(count=n()) #23 genus

df_fam<- df %>%
  filter(!Family == "na")%>%
  group_by(Family, Shark_Ray_Chimaera)%>%
  summarise(count=n()) #15 families

df_ord<- df %>%
  group_by(Order, Shark_Ray_Chimaera)%>%
  summarise(count=n()) #9 orders

df_tax <- df %>%
  filter(!is.na(Family)) %>% 
  filter(!Family == "na")%>%
  group_by(molecular_marker, Family) %>% 
  summarise(count = n())
unique(df_tax$molecular_marker) #8 markers

df_rp<- df %>%
  separate_rows(Rapid.DNA.Tool, sep = "; ")
unique(df_rp$Rapid.DNA.Tool) #8 different rapid tools

df_rp_c<- df_rp %>%
  group_by(Rapid.DNA.Tool, molecular_marker)%>%
  summarise(count=n()) 

####Figure 2A - Spider plot marker vs application####

#Filtering data - all loci
mk_spider<- df%>%
  filter(molecular_marker %in% c("12S", "COI", "CytB", "D-loop", "ITS2", "ND2", "ND4", "ND5"))%>%
  drop_na()

#couting markers
tool_summary<- df %>%
  group_by(molecular_marker) %>%
  summarize(count_by_tool = n())

write.csv(tool_summary, "marker_count.csv")

#Counting papers by application
tool_spider_summary<- mk_spider %>%
  filter(!Application == "na") %>%
  group_by(molecular_marker, Application) %>%
  summarize(count_by_tool = n())

write.csv(tool_spider_summary, "marker_application_count.csv")

#Reshaping data - long to wide
tool_spider_rs<- tool_spider_summary %>%
  spread(key=molecular_marker, value=count_by_tool) %>%
  replace(is.na(.), 0) 

#Adjusting before plotting
tool_spider_rs <- tool_spider_rs %>%
  remove_rownames %>% 
  column_to_rownames(var="Application")
# Adding maximum and minimum to the data 
tool_spi <- rbind(rep(100, 8), rep(0, 8), tool_spider_rs)

#Spider plot by filter size
colors_bord<- c("deeppink4", "lightpink")
colors_in<- c("deeppink4", "lightpink")

pdf(file = "spideplot_app_1.pdf", width = 10, height = 7)

radarchart(tool_spi, axistype = 4, 
           maxmin = TRUE,
           pcol = colors_bord, 
           pfcol = scales::alpha(colors_in, 0.1), 
           plwd = 3, 
           plty = 1,
           # Custom grid with adjusted labels and manual lines
           cglcol = "black", 
           cglty = 3, 
           axislabcol = "black", 
           cglwd = 1.2,
           caxislabels = c("0", "25", "50", "75", "100"),  # Directly set the labels as text
           vlcex = 1.2)

legend(x = "bottom", legend = rownames(tool_spi[-c(1, 2),]), 
       horiz = TRUE, bty = "n", pch = 20, 
       col = colors_in, 
       text.col = "black", 
       cex = 1.3, pt.cex = 1.8, 
       x.intersp = 0.3, y.intersp = 0.5) 

dev.off()


####Figure 2B - Donut plot - sharks vs rays####

#percentages and cumulative 
df_sr$fraction <- df_sr$count / sum(df_sr$count)
df_sr$ymax <- cumsum(df_sr$fraction)
df_sr$ymin <- c(0, head(df_sr$ymax, n=-1))

#create donut plot
tax_don<-ggplot(df_sr, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Shark_Ray_Chimaera)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2, 5)) + 
  theme_void() +
  scale_fill_manual(values = c("Shark" = "gray34", "Ray" = "lightcyan4"),
                     labels = c("Shark", "Ray")) +
  theme(legend.position = "none")

tax_don

ggsave('donut_sharkray.pdf', tax_don,
       width = 9, height = 6, units = c('cm'),
       dpi = 600)

####Figure 2C - Heat map - family vs marker####

df_tax <- df_tax %>%
  mutate(Family = factor(Family, levels = sort(unique(Family))))

#reshaping data format
df_wide <- dcast(df_tax, Family ~ molecular_marker, value.var = "count", fill = 0)

#creating matrix
mat <- as.matrix(df_wide[,-1])
rownames(mat) <- df_wide$Family

df_melt <- melt(mat)
colnames(df_melt) <- c("Family", "Molecular_Marker", "Count")

df_melt$Family <- trimws(df_melt$Family)

df_melt$Family <- factor(df_melt$Family, levels = sort(unique(df_melt$Family)))

fac_order <- c("Trygonorrhinidae", "Rhinobatidae", "Mobulidae", "Dasyatidae", 
               "Rajidae", "Pristidae", "Rhincodontidae", "Squatinidae", "Heterodontidae",
               "Cetorhinidae", "Triakidae", "Lamnidae", "Sphyrnidae", "Carcharhinidae", "Alopiidae")

#converting Rapid_DNAtools_c to a factor with the specified levels
df_melt <- df_melt %>%
  mutate(Family2 = factor(Family, levels = fac_order))

#plottinh heatmap
tax_hm<-ggplot(df_melt, aes(x = Molecular_Marker, y = Family2, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(
    x = "Molecular Marker",
    y = "Family",
    fill = "Primer count" 
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axix.text.y = element_text(size = 10))

plot(tax_hm)

ggsave('heatmap_taxon.pdf', tax_hm,
       width = 15, height = 12, units = c('cm'),
       dpi = 600)
