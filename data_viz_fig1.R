#Loading packages
library(tidyverse)
library(stringr)
library(sf)
library(ggplot2)
library(packcircles)
library(maps)
library(treemap)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(patchwork)

#Setting directory
setwd("/Users/ingridbunholi/Desktop/Bioinformatics/DNArapidtools/metadata")

#Setting metadata
df<- read.csv("allsets_rapid_tools_v2.csv")

####Figure 1A - Map + donut plots####

#Getting shape file for both maps - maps package
world_map<- map_data("world") #including Antarctica

#Separating papers with more than 1 country (one each row)
df.country.first<- df%>%
  mutate(country = strsplit(as.character(Country_first_author_c), "; ")) %>%
  unnest(country)
unique(df.country.first$country) 

df.country.first<- df.country.first%>%
  mutate(country = recode(country,
                          `United States` = "USA"))

#Counting papers
map.country.first<- df.country.first %>%
  group_by(country) %>%
  summarise(count = n())

#Example data preparation for donut plots
donut_data <- df.country.first %>%
  filter(country %in% c("USA", "Brazil", "Australia", "United Kingdom")) %>%
  count(country, Application_primer_2) %>%
  group_by(country) %>%
  mutate(percentage = n / sum(n) * 100)

donut_data2 <- df.country.first %>%
  #filter(country %in% c("USA", "Brazil", "Australia", "China", "United Kingdom")) %>%
  count(country, Application_primer_2) %>%
  group_by(country) %>%
  mutate(percentage = n / sum(n) * 100)

#Create a list to store grobs (donut plots)
donut_grobs <- list()

category_colors <- c("lightpink", "deeppink4")

#Create a grob for each country
for (country in unique(donut_data$country)) {
  p <- ggplot(donut_data %>% filter(country == !!country), aes(x = "", y = percentage, fill = Application_primer_2)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = category_colors) +
    theme_void() +
    theme(legend.position = "none")  # Hide legend if you want to customize it separately
  
  # Convert the plot to a grob
  grob <- ggplotGrob(p)
  
  # Store the grob with its corresponding country
  donut_grobs[[country]] <- grob
}

#Create the map plot
map_plot <- ggplot(df.country.first) +
  geom_map(
    data = world_map, map = world_map, aes(map_id = region),
    fill = "lightgray", color = "white", size = 0.25) +
  geom_map(data = map.country.first, map = world_map, aes(map_id = country, fill = count)) +
  scale_fill_gradient(name = "Studies per country", low = "lightblue", high = "darkblue") +
  expand_limits(x = world_map$long, y = world_map$lat) +
  theme_void() +  # Remove the default background and grid
  theme(legend.title = element_text(size=10, face= "bold"),
        legend.text = element_text(size=8),
        legend.position=c(.14,0.4),
        legend.key.size=unit(0.6, 'cm'),
        plot.margin = unit(c(0.5, 0, -0.3, -0.2), 'cm'))

#Define positions for the donut plots
positions <- data.frame(
  country = c("USA", "Brazil", "Australia", "United Kingdom"),
  x = c(-100, -55, 135, 0), #this can be changed at illustrator 
  y = c(40, -10, -25, 55))  

#Add donut plots to the map plot
for (i in 1:nrow(positions)) {
  map_plot <- map_plot +
    annotation_custom(donut_grobs[[positions$country[i]]],
                      xmin = positions$x[i] - 5, xmax = positions$x[i] + 5,
                      ymin = positions$y[i] - 5, ymax = positions$y[i] + 5)
}

print(map_plot)

ggsave('map_author_new.pdf', map_plot,
       width = 21, height = 13, units = c('cm'),
       dpi = 600)

####Figure 1B - Time series - eDNA eRNA studies over time####

data<- df

##Setting time series graphical elements
set<- theme(axis.text.x=element_text(angle=45,vjust = 0.9, hjust=1, size = 11),
            axis.text.y=element_text(size = 11),
            axis.ticks = element_line(linewidth = 0.5),
            axis.title.y = element_text(size = 13, vjust = 2),
            legend.title = element_blank(),
            legend.text = element_blank())

#Checking years
unique(data$year)
data$year<- as.character(data$year)

#Counting DNA_RNA studies per year
data_time<- data%>%
  group_by(year, Application_primer_2)%>%
  summarise(count = n())

# Convert year to numeric
data_time <- data_time %>%
  mutate(year = as.numeric(year))

#Plot time series eDNA/eRNA studies over time
time<- ggplot(data_time, aes(x=year, y=count, color=Application_primer_2, group=Application_primer_2))+
  geom_line(linewidth=1.5)+
  geom_point()+
  labs(x="", y="Number of studies")+
  theme_classic()+
  scale_x_continuous(breaks = c(2001, 2005, 2010, 2015, 2020, 2024)) +
  #scale_x_discrete(breaks=data_time$year)+
  scale_color_manual(values = c("deeppink4", "lightpink"))+
  #scale_y_continuous(breaks =c(0,10,20,30,40,50,60,70,80))+
  set
time
#Saving plot
ggsave('timeseries.pdf', time,
       width = 9, height = 6, units = c('cm'),
       dpi = 600)


####Figure 1C - lollipop plot DNA rapid tools and application ####

#Count occurrences of each "Rapid_DNAtools" for each "Application_primer"
df_unique<- df%>%
  separate_rows(Rapid_DNAtools_c, sep = "; ")

data_counts_rp <- df_unique %>%
  group_by(Rapid_DNAtools_c, Application_primer_2) %>%
  summarise(count = n()) %>%
  ungroup()

#Reshape the data for plotting
data_wide_rp <- data_counts_rp %>%
  pivot_wider(names_from = Application_primer_2, values_from = count, values_fill = list(count = 0))

data_wide_rp <- data_wide_rp %>%
  rename(V1 = eDNA, V2 = `Tissue DNA`)

desired_order<- c("on-site identification","RT-HRM-PCR", "rep-PCR", "Taxon-Specific PCR ", 
                  "PCR RFLP", "Droplet Digital PCR (ddPCR)", "Real Time PCR (qPCR)", "Multiplex qPCR", "Multiplex PCR")

#Convert Rapid_DNAtools_c to a factor with the specified levels
data_wide_rp2 <- data_wide_rp %>%
  mutate(Rapid_DNAtools_c = factor(Rapid_DNAtools_c, levels = desired_order))

#lollipop plot
rapid_loli_rp2 <- ggplot(data_wide_rp2, aes(x = Rapid_DNAtools_c)) +
  geom_segment(aes(xend = Rapid_DNAtools_c, y = V2, yend = V1), color = "grey", size = 1) +
  geom_point(aes(y = V1, color = "eDNA"), size = 5) +
  geom_point(aes(y = V2, color = "Tissue"), size = 5) +
  geom_text(aes(y = V1, label = V1), vjust = -1, size = 4) +
  geom_text(aes(y = V2, label = V2), vjust = -1, size = 4) +
  scale_color_manual(values = c("eDNA" = "deeppink4", "Tissue" = "lightpink"),
                     labels = c("eDNA", "Tissue DNA")) +
  coord_flip() +
  theme_minimal() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  ylim(0, max(data_wide_rp$V1, data_wide_rp$V2) * 1.1) +
  ylab("Paper Count") +
  xlab("")

ggsave('rapid_application_lollipop.pdf', rapid_loli_rp2,
       width = 20, height = 12, units = c('cm'),
       dpi = 600)
