#install the below packages using install.packages(c("package1", "package2"))

library(tidyverse)
library(readr)
library(openxlsx)
library(plotrix)
library(rstatix)
library(Hmisc)
library(ggpubr)
library(ggrepel)

read_csv("PTS3886_AHA_Surf4_protein_filtered.csv") ->PTS3886_Annot
read_csv("PTS3886_Subcell_annot.csv") -> PTS3886

#Merge with annotations and write
inner_join(PTS3886_Annot, PTS3886, by = "Protein_names") -> PTS3886full_annot

#PTS3886 plot
read_csv("PTS3886final.csv")->data

#FDR-adjust p-value
data$adj_p_val_fdr <- p.adjust(data$t_test, "fdr")

#add treshold for original p_values
data$significant1 = as.factor(data$t_test < 0.05)

#labels to display
data$absLog2FC <- abs(data$Log2FC)

#to avoid displaying cytosolic protein names do not use abs
SILAC_labels <- filter(data, Log2FC > 0.5 , t_test<0.05)
SILAC_labels_fdr <- filter(data, Log2FC > 0.5 , adj_p_val_fdr<0.05)

#add colour thresholds
PTS3886colour_cyt <- data %>% filter(str_detect(Location, "Cytoplasm"))
PTS3886colour_nuc <- data %>% filter(str_detect(Location, "Nucleus"))
PTS3886colour_secreted <- data %>% filter(str_detect(Location, "Secreted"))
PTS3886colour_ER <- data %>% filter(str_detect(Location, "reticulum"))
PTS3886colour_Golgi <- data %>% filter(str_detect(Location, "Golgi"))
PTS3886colour_lysosome <- data %>% filter(str_detect(Location, "Lysosome"))
PTS3886colour_PM <- data %>% filter(str_detect(Location, "Cell membrane"))

#Plot for non-adjusted p_values
ggplot(data, 
       aes(Log2FC , -log10(t_test)))+
  geom_vline(xintercept = 0.5, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype="dashed", color = "grey50")+
  geom_hline(yintercept = 1.3, linetype="dashed", color = "grey50") + 
  geom_text_repel(data = SILAC_labels,
                  aes(x=Log2FC,
                      y=-log10(t_test),
                      label = `Gene_names`))+
  geom_point(data = PTS3886colour_PM,            #have to use variable names here!
             aes(x = Log2FC, y = -log10(t_test)),
             colour="deepskyblue3") +
  geom_point(data = PTS3886colour_cyt,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="pink") +
  geom_point(data = PTS3886colour_nuc,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="pink") +
  geom_point(data = PTS3886colour_ER,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="bisque2") +
  geom_point(data = PTS3886colour_lysosome,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="darksalmon") +
  geom_point(data = PTS3886colour_Golgi,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="darkseagreen3") +
  geom_point(data = PTS3886colour_secreted,            
             aes(x = Log2FC, y = -log10(t_test)),
             colour="turquoise")

#save the last plot 
ggsave(filename = "PTS3886_bylocation.jpeg", device = "jpeg", path = NULL, width = 200,
       height = 150, units = c("mm"), dpi = 600)

#save the new csv file
write.csv(PTS3886full_annot, "PTS3886final.csv")
