## Figure 2 - Peru sequence plots per week (fill and stacked) 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ape)

Peru_cov <- read.csv("Final_Peru_metadata.csv") ## Read metadata
variants<-read.csv("Variants.csv") # Read variants file with colours

variant_tip<-data.frame(Tip=Peru_cov$strain,Variant="Other",Color="lightgrey") # new dataframe

for (i in 1:nrow(variants)){
  variant_tip$Variant[grep(variants$Lineages[i],Peru_cov$pangolin_lineage)]<-variants$Variant[i]
  variant_tip$Color[grep(variants$Lineages[i],Peru_cov$pangolin_lineage)]<-variants$Colour[i]
}

# Change variant labels
variant_tip$Variant[which(variant_tip$Variant=="BA.1")]<-"Omicron_BA.1"
variant_tip$Variant[which(variant_tip$Variant=="BA.2")]<-"Omicron_BA.2"
variant_tip$Variant[which(variant_tip$Variant=="BA.4")]<-"Omicron_BA.4"
variant_tip$Variant[which(variant_tip$Variant=="BA.5")]<-"Omicron_BA.5"
variant_tip$Variant[which(variant_tip$Variant=="Delta")]<-"Other Delta"
variant_tip$Variant<-as.factor(variant_tip$Variant)
colors<-unique(variant_tip$Color)
names(colors)<-unique(variant_tip$Variant)

## Add date to variant_tip dataframe
variant_tip<-variant_tip[order(match(variant_tip$Tip,Peru_cov$strain)),]
variant_tip$date<-Peru_cov$date
variant_tip$date<-as.Date(variant_tip$date)

weekly_summary <- variant_tip %>%
  group_by(week = floor_date(date, "week"), Variant) %>%
  summarise(count = n(), .groups = 'drop')
p1<-ggplot(weekly_summary, aes(x = week, y = count, fill = Variant)) +
  geom_area(alpha=0.75,position = 'identity') + scale_fill_manual(values=colors) +
  labs(
    x = "Week",
    y = "Count",
    color = "Variant") +
  theme_minimal()

## Stacked proportion
weekly_summary_prop <- variant_tip %>%
  group_by(week = floor_date(date, "week"), Variant) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(week) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p2<-ggplot(weekly_summary_prop, aes(x = week, y = proportion, fill = Variant)) +
  geom_area(position = "fill", alpha = 0.8) +  scale_fill_manual(values=colors) +
  labs(
    x = "Week",
    y = "Proportion",
    fill = "Variant") +
  theme_minimal()

tiff("Weekly_variants.tiff",width = 20,height=16,units = "cm",res = 300)
ggarrange(p1,p2,ncol=1,nrow=2,common.legend = T,legend = 'right')
dev.off()

