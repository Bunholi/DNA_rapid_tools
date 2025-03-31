##############################################################
#### Figure 1 - world map/time series/lolipop rapid tools ####
##############################################################

#Loading packages
library(tidyverse)
library(stringr)
library(ggplot2)
library(sf)
library(maps)
library(treemap)
library(patchwork)
library(RColorBrewer)
library(cowplot)

#Setting directory - path to your directory
setwd("/Users/YourDirectory")

#Setting metadata
df<- read.csv("allsets_rapid_tools_v2.csv")

####Figure 1A - Map + donut plots####

#Getting shape file for both maps - maps package
world_map<- map_data("world") #including Antarctica

#Separating papers with more than 1 country (one each row)
df.country.first<- df%>%
  mutate(country = strsplit(as.character(country_first_author_c), "; ")) %>%
  unnest(country)
unique(df.country.first$country) 

df.country.first<- df.country.first%>%
  mutate(country = recode(country,
                          `United States` = "USA"))

#Counting papers
map.country.first<- df.country.first %>%
  group_by(country) %>%
  summarise(count = n())

map.country.first.2<- df.country.first %>%
  group_by(country, application_primer_2) %>%
  summarise(count = n())

write.csv(map.country.first.2, "map.country.first.2.csv")

#Example data preparation for donut plots
donut_data <- df.country.first %>%
  filter(country %in% c("USA", "Brazil", "Australia", "Mexico", "United Kingdom")) %>%
  count(country, application_primer_2) %>%
  group_by(country) %>%
  mutate(percentage = n / sum(n) * 100)

donut_data2 <- df.country.first %>%
  #filter(country %in% c("USA", "Brazil", "Australia", "China", "United Kingdom")) %>%
  count(country, application_primer_2) %>%
  group_by(country) %>%
  mutate(percentage = n / sum(n) * 100)

#Create a list to store grobs (donut plots)
donut_grobs <- list()

category_colors <- c("lightpink", "deeppink4")

#Create a grob for each country
for (country in unique(donut_data$country)) {
  p <- ggplot(donut_data %>% filter(country == !!country), aes(x = "", y = percentage, fill = application_primer_2)) +
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
  country = c("USA", "Brazil", "Australia", "Mexico", "United Kingdom"),
  x = c(-100, -55, 135, -102, 0), #this can be changed in Illustrator
  y = c(40, -10, -25, 23, 55))  


#Add donut plots to the map plot
for (i in 1:nrow(positions)) {
  map_plot <- map_plot +
    annotation_custom(donut_grobs[[positions$country[i]]],
                      xmin = positions$x[i] - 5, xmax = positions$x[i] + 5,
                      ymin = positions$y[i] - 5, ymax = positions$y[i] + 5)
}

print(map_plot)

ggsave('map_author_033125.pdf', map_plot,
       width = 21, height = 13, units = c('cm'),
       dpi = 600)

####Figure 1B - Time series - eDNA eRNA studies over time####

data<- df

#Counting DNA_RNA studies per year
data_time<- data%>%
  group_by(year, application_primer_2)%>%
  summarise(count = n())

#Convert year to numeric
data_time$year <- as.numeric(as.character(data_time$year))

#UNGROUP in case it's grouped from previous operations
data_time <- data_time %>% ungroup()

#Define full range of years
year_range <- seq(min(data_time$year), max(data_time$year), by = 1)

#Make sure both columns are present before completing
data_time <- data_time %>%
  complete(
    year = year_range,
    application_primer_2 = unique(data_time$application_primer_2),
    fill = list(count = 0)
  )

#Find the first non-zero year per group
first_years <- data_time %>%
  filter(count > 0) %>%
  group_by(application_primer_2) %>%
  summarise(min_year = min(year), .groups = "drop")

#Replace early zeros with NA to prevent flatlines
data_time_filtered <- data_time %>%
  left_join(first_years, by = "application_primer_2") %>%
  mutate(count = ifelse(year < min_year, NA, count))

#Setting time series graphical elements
set<- theme(axis.text.x=element_text(angle=45,vjust = 0.9, hjust=1, size = 11),
            axis.text.y=element_text(size = 12),
            axis.title.y = element_text(size = 13, vjust = 2),
            legend.title = element_blank(),
            legend.text = element_blank())

#Plotting time series eDNA and Tissue DNA studies
time<- ggplot(data_time_filtered, aes(x = year, y = count, color = application_primer_2)) +
  geom_line(size = 1.2, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = c("deeppink4", "lightpink")) +
  scale_x_continuous(
    limits = c(2001, 2025),
    breaks = c(2001, seq(2005, 2025, by = 5))
  ) +
  labs(x = NULL, y = "Number of studies", color = NULL) +
  theme_minimal(base_size = 11) +
  set

time
#Saving plot
ggsave('timeseries.pdf', time,
       width = 12, height = 9, units = c('cm'),
       dpi = 600)


####Figure 1C - Lollipop plot - DNA rapid tools and application ####

#Count occurrences of each "Rapid_DNAtools" for each "Application_primer"
df_unique<- df%>%
  separate_rows(rapid_DNAtools_c, sep = "; ")

data_counts_rp <- df_unique %>%
  group_by(rapid_DNAtools_c, application_primer_2) %>%
  summarise(count = n()) %>%
  ungroup()

#Reshape the data for plotting
data_wide_rp <- data_counts_rp %>%
  pivot_wider(names_from = application_primer_2, values_from = count, values_fill = list(count = 0))

data_wide_rp <- data_wide_rp %>%
  rename(V1 = eDNA, V2 = `Tissue DNA`)

desired_order<- c("on-site identification","RT-HRM-PCR", "rep-PCR", "Taxon-Specific PCR ", 
                  "PCR RFLP", "Droplet Digital PCR (ddPCR)", "Real Time PCR (qPCR)", "Multiplex qPCR", "Multiplex PCR")

#Convert Rapid_DNAtools_c to a factor with the specified levels
data_wide_rp2 <- data_wide_rp %>%
  mutate(rapid_DNAtools_c = factor(rapid_DNAtools_c, levels = desired_order))

#lollipop plot
rapid_loli_rp2 <- ggplot(data_wide_rp2, aes(x = rapid_DNAtools_c)) +
  geom_segment(aes(xend = rapid_DNAtools_c, y = V2, yend = V1), color = "grey", size = 1) +
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
    legend.text = element_blank()
  ) +
  ylim(0, max(data_wide_rp$V1, data_wide_rp$V2) * 1.1) +
  ylab("Paper Count") +
  xlab("")

rapid_loli_rp2

ggsave('rapid_application_lollipop.pdf', rapid_loli_rp2,
       width = 16, height = 10, units = c('cm'),
       dpi = 600)
