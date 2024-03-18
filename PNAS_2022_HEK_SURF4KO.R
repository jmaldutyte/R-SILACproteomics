#install the below packages using install.packages(c("package1", "package2"))

library(plotrix)
library(tidyverse)
library(rstatix)
library(Hmisc)
library(ggpubr)
library(ggrepel)
library(Rcpp)
library(sjlabelled)

PTS5107 <- read_csv("Perseus-one-sample-t-testR.csv")

#adjust p-value
PTS5107$adj_p_val_fdr <- p.adjust(PTS5107$t_test, "fdr")

#add threshold for original p_values
PTS5107$significant1 = as.factor(PTS5107$t_test < 0.05)
PTS5107$significant_adj = as.factor(PTS5107$adj_p_val_fdr < 0.05)

#labels to display
PTS5107$absLog2_FC <- abs(PTS5107$Log2_FC)
PTS5107labels <- filter(PTS5107, Log2_FC>1.25 , t_test<0.05)

#test plot, original p-values
ggplot(PTS5107, 
       aes(Log2_FC , -log10(t_test),
           colour=significant1)) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = 0.5, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype="dashed", color = "grey50")+
  geom_hline(yintercept = 1.3, linetype="dashed", color = "grey50") +
  aes(x=Log2_FC,
      y=-log10(t_test))

#test plot, FDR-adjusted p-values
ggplot(PTS5107, 
       aes(Log2_FC , -log10(adj_p_val_fdr),
           colour=significant_adj)) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = 0.5, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype="dashed", color = "grey50")+
  geom_hline(yintercept = 1.3, linetype="dashed", color = "grey50") +
  aes(x=Log2_FC,
      y=-log10(adj_p_val_fdr))

#add colour thresholds
PTS5107colour_cyt <- PTS5107 %>% filter(str_detect(Location, "Cytoplasm"))
PTS5107colour_nuc <- PTS5107 %>% filter(str_detect(Location, "Nucleus"))
PTS5107colour_secreted <- PTS5107 %>% filter(str_detect(Location, "Secreted"))
PTS5107colour_ER <- PTS5107 %>% filter(str_detect(Location, "reticulum"))
PTS5107colour_Golgi <- PTS5107 %>% filter(str_detect(Location, "Golgi"))
PTS5107colour_lysosome <- PTS5107 %>% filter(str_detect(Location, "Lysosome"))
PTS5107colour_PM <- PTS5107 %>% filter(str_detect(Location, "Cell membrane"))
PTS5107colour_unknown <- PTS5107 %>% filter(str_detect(Location, "Unknown"))
PTS5107colour_stressgranule <- PTS5107 %>% filter(str_detect(Location, "Stress granule"))
PTS5107colour_mito <- PTS5107 %>% filter(str_detect(Location, "Mitochondrion"))

#plot with original FDR-adjusted p-values, dashed lines indicates statistical and fold-change significance boundaries
ggplot(PTS5107, 
       aes(Log2_FC , -log10(adj_p_val_fdr())))+
  geom_vline(xintercept = 1.25, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = -1.25, linetype="dashed", color = "grey50")+
  geom_hline(yintercept = 1.3, linetype="dashed", color = "grey50") + 
  geom_text_repel(data = PTS5107labels,
                  size=3,
                  aes(x=Log2_FC,
                      y=-log10(adj_p_val_fdr),
                      label = `Gene names`)) +
  geom_point(data = PTS5107colour_PM,            #have to use variable names here!
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="deepskyblue3",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_cyt,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="pink",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_nuc,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="pink",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_ER,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="bisque2",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_lysosome,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="darksalmon",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_Golgi,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="darkseagreen3",
             alpha = 0.5) +
  geom_point(data = PTS5107colour_secreted,            
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="turquoise",
             alpha = 0.5) +
  geom_point(data=PTS5107colour_unknown,
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="grey",
             alpha = 0.5) +
  geom_point(data=PTS5107colour_stressgranule,
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="grey",
             alpha = 0.5) +
  geom_point(data=PTS5107colour_mito,
             aes(x = Log2_FC, y = -log10(adj_p_val_fdr)),
             colour="sandybrown",
             alpha = 0.5)

#filter out labels to exclude cytosolic proteins, RPL19 is cut-off
PTS5107labels_adj_RPL19 <- filter(PTS5107, Log2_FC>1.25 , adj_p_val_fdr<0.027)

#unadjusted p-values but labels displayed only for FDR-adjusted significant ones
#that are also Log2_FC>1.25. Secreted, Secretory pathway, Other, white background
ggplot(PTS5107, 
       aes(Log2_FC , -log10(t_test())))+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_text_repel(data = PTS5107labels_adj_RPL19,
                  size=3,
                  nudge_x = 0.3,
                  aes(x=Log2_FC,
                      y=-log10(t_test),
                      label = `Gene names`)) +
  geom_point(data = PTS5107colour_PM,           
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="darkseagreen3",
             alpha = 1) +
  geom_point(data = PTS5107colour_cyt,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="pink",
             alpha = 1) +
  geom_point(data = PTS5107colour_nuc,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="pink",
             alpha = 1) +
  geom_point(data = PTS5107colour_ER,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="darkseagreen3",
             alpha = 1) +
  geom_point(data = PTS5107colour_lysosome,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="darkseagreen3",
             alpha = 1) +
  geom_point(data = PTS5107colour_Golgi,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="darkseagreen3",
             alpha = 1) +
  geom_point(data=PTS5107colour_unknown,
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="pink",
             alpha = 1) +
  geom_point(data=PTS5107colour_stressgranule,
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="pink",
             alpha = 1) +
  geom_point(data=PTS5107colour_mito,
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="pink",
             alpha = 1) +
  geom_point(data = PTS5107colour_secreted,            
             aes(x = Log2_FC, y = -log10(t_test)),
             colour="turquoise",
             alpha = 1) +
  theme_pubr()

#saves the last plot and a new, updated csv file
ggsave(filename = "Final_HEK_SURF4KO.jpeg", device = "jpeg", path = NULL, width = 200, height = 150, units = c("mm"), dpi=700)
write.csv(PTS5107,"C:/Users/jmald/silac/PTS_5107_September_2020_3-day_split/PTS5107_export.csv", row.names = FALSE)

