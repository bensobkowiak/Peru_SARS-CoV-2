### Ancestral reconstruction of Peru regions with ace and alluvial plots and colored phylogeny

## Install packages
packages <- c("ape", "ggplot2", "treeio","tidytree","lubridate",
              "tidyverse","dplyr", "RColorBrewer","ggtree","reshape2","ggalluvial")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])}
invisible(lapply(packages, library, character.only = TRUE))
rm(installed_packages, packages)

## Read tree and clean names (drop outgroup)
perutree <- read.nexus("timetree.nexus")
perutree$tip.label<-gsub("'",replacement = "",perutree$tip.label)
perutree<-drop.tip(perutree,tip=c("MN908947.3"))

# Read in metadata and clean region names
fullmeta <- read.csv("Final_Peru_metadata.csv")
fullmeta<-fullmeta[order(match(fullmeta$strain,perutree$tip.label)),]
fullmeta$division[which(fullmeta$division=="Amazona")]<-"Amazonas"
fullmeta$division[which(fullmeta$division=="Cuzco")]<-"Cusco"
fullmeta$division[which(fullmeta$division=="Hu\xe1nuco")]<-"Huanuco"
fullmeta$division[which(fullmeta$division=="Madre de Dios")]<-"Madre De Dios"
fullmeta$division[which(fullmeta$division=="HuancavelIca")]<-"Huancavelica"
fullmeta$division[which(fullmeta$division=="La libertad")]<-"La Libertad"

# Set colours
prpal<-read.csv("Region_color_scheme.csv")
color_vec<-prpal$Color
names(color_vec)<-prpal$Region

# Ancestral state reconstruction with ace
region_vec<-fullmeta$division
names(region_vec)<-fullmeta$strain

perutree<-multi2di(perutree) # Make bifuricating and add non-zero branches
zerobranchesr<-which(perutree$edge.length==0)
if (length(zerobranchesr)>0){
  perutree$edge.length[zerobranchesr]<-0.0000000001
}

acl_region<-ace(region_vec,perutree, type="d") # Anc rec 
write_rds(acl_region,"Acl_region_new2.RData") # Save result


#Get the most likely country for each node
regionorder<-colnames(acl_region$lik.anc)
rmax <- regionorder[apply(acl_region$lik.anc, 1, which.max)]
names(rmax)<-row.names(acl_region$lik.anc)

# check treedata object with inferred node location
r_assignment_df<-as_tibble(perutree)
r_assignment_df$reg_tip<-region_vec[match(r_assignment_df$label,names(region_vec))]
r_assignment_df$reg_node<-rmax[match(r_assignment_df$node,names(rmax))]
r_assignment_df<-r_assignment_df %>% 
  mutate(reg_all = coalesce(reg_tip,reg_node))
r_assignment_df$reg_tip<-NULL
r_assignment_df$reg_node<-NULL
node_times <- max(node.depth.edgelength(perutree)) - node.depth.edgelength(perutree)
r_assignment_df$Date<-decimal_date(as.Date(max(fullmeta$date)))-node_times
rlda_data<-as.treedata(r_assignment_df)


### Plot phylogeny with tips coloured by variant
peru_reg<-ggtree(rlda_data,right=TRUE, mrsd=max(fullmeta$date),aes(color=reg_all))+
  scale_color_manual(values=color_vec,
                     labels=names(color_vec))+
  theme_tree2(legend.title=element_blank(),text = element_text(size = 20)) + 
  geom_point(size=0) +
  guides(colour = guide_legend(override.aes = list(size=5,linetype=0)))

pdf("Peru_timed_tree.pdf",height = 14,width = 12)
plot(peru_reg)
dev.off()

## Make data frame counting any movements between regions from nodes
rasf<-as_tibble(r_assignment_df)
regcrossings<-r_assignment_df[0,]
colnames(regcrossings)<-c("Child","ChildCountry","ChildDate","Parent","ParentCountry","ParentDate")
for (i in 1:nrow(r_assignment_df)){
  #get tip.label, associated node, and country
  tip<-as.character(r_assignment_df[i,4]) #label column is 4
  tipctry<-as.character(r_assignment_df[i,5]) #node column is 2
  tipdate<-as.numeric(r_assignment_df[i,6])
  #get parent node and country
  pnode<-as.character(r_assignment_df[i,1]) #parent column is 1
  pctry<-as.character(r_assignment_df[pnode,5])
  pdate<-as.numeric(r_assignment_df[pnode,6])
  vec<-c(tip,tipdate,tipctry,pnode,pctry,pdate)
  #if they are not the same, store all of that into a dataframe
  if(tipctry!=pctry){
    regcrossings=rbind(regcrossings,vec)
  }
}
colnames(regcrossings)<-c("Child","ChildRegion","ChildDate","Parent","ParentRegion","ParentDate")

write.table(regcrossings[,c(2,3,5,6)],"Fig3_AllPeru_RegionDateMovements.txt",quote = F,sep = "\t",row.names = F)

### Count all movements from each department
rcr<-regcrossings[,c(3,5)]
rcr_freq <- as.data.frame(table(rcr))
colnames(rcr_freq)<-c("Destination","Origin","Freq")
write.csv(rcr_freq,"Departmental_movements.csv",row.names = F)

## Calculate proportion of movements from each department
depts_movements<-data.frame(Dept=as.character(unique(rcr_freq$Origin)),Origin=0,Dest=0)
for (i in 1:nrow(depts_movements)){
  depts_movements$Origin[i]<-round(sum(rcr_freq$Freq[which(rcr_freq$Origin==depts_movements$Dept[i])])/sum(rcr_freq$Freq),3)*100
  depts_movements$Dest[i]<-round(sum(rcr_freq$Freq[which(rcr_freq$Destination==depts_movements$Dept[i])])/sum(rcr_freq$Freq),3)*100
}
write.csv(depts_movements,"Departmental_proportion_movements.csv",row.names = F)

## Movements from lima
rcr_freq_lima<-rcr_freq[which(rcr_freq$Origin=="Lima"),]
depts_movements_lima<-data.frame(Dept=as.character(unique(rcr_freq_lima$Destination)),Dest=0,seqsProp=0)
for (i in 1:nrow(depts_movements_lima)){
  depts_movements_lima$Dest[i]<-round(sum(rcr_freq_lima$Freq[which(rcr_freq_lima$Destination==depts_movements_lima$Dept[i])])/sum(rcr_freq_lima$Freq),3)*100
  depts_movements_lima$seqsProp[i]<-round(sum(rcr_freq_lima$Freq[which(rcr_freq_lima$Destination==depts_movements_lima$Dept[i])])/length(which(fullmeta$division==depts_movements_lima$Dept[i])),3)*100
}
write.csv(depts_movements_lima,"Departmental_proportion_movements_lima.csv",row.names = F)


## Convert movements data frame into alluvial object
rcr_freq$Origin<-factor(rcr_freq$Origin,levels=names(color_vec))
rcr_freq$Destination<-factor(rcr_freq$Destination,levels=names(color_vec)) #625
rcr_freq_clean<-rcr_freq[rcr_freq$Freq>0,] #477
rcr_freq_clean$subject<-row.names(rcr_freq_clean)
rcrm<-reshape2::melt(rcr_freq_clean,id.vars=c("Freq","subject"))
rcrm$value<-factor(rcrm$value,levels=names(color_vec)) 
is_alluvia_form(rcrm, axes = 2, silent = TRUE)

### Make alluvial plot
all_reg2<-ggplot(rcrm,
                 aes(x=variable, stratum = value, y = Freq, 
                     fill=value, alluvium=subject, label=value)) +
  geom_stratum(width = 1/6, alpha=0.8) +
  geom_flow(width = 1/6)+
  geom_text(stat = "stratum", size=4, min.y = 100) +
  scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.01, .01)) +
  scale_fill_manual(values = color_vec, labels=names(color_vec))+
  ylab("Count of regional border crossings")+
  xlab("")+
  labs(fill="Region")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(vjust=10,size=15),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 15))

tiff("Peru_regionalmovements_alluvial.tiff",width = 25,height = 30,units = "cm",res = 300)
all_reg2
dev.off()


