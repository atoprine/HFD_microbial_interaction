#-----------------------------------------------------------------------#
# Figure 2               
#-----------------------------------------------------------------------#
library(tidyverse)
library(cowplot)

uvmr_list <- read_rds("./2.2.1_panCancer_genus_MR_result.rds")

uvmr <- lapply(names(uvmr_list), function(x){
  tmp <- uvmr_list[[x]] %>% 
    filter(pval < 0.05) %>% 
    arrange(exposure) %>% 
    mutate(type = x,
           trait = paste(cancer_type, exposure, sep = "|")) %>% 
    select(type, exposure, trait, pval, everything())
}) %>% 
  do.call(rbind, .)

ind <- as.data.frame(table(uvmr$trait)) %>% 
  rename(trait = Var1) %>%
  arrange(desc(Freq)) %>% 
  filter(Freq > 3) %>% 
  mutate(trait = as.character(trait)) %>% 
  pull(trait)


df <- do.call(rbind, uvmr_list) %>% 
  mutate(trait = paste(cancer_type, exposure, sep = "|")) %>% 
  select(trait, or, pval, method) %>% 
  filter(trait %in% ind) %>% 
  mutate(lgOR  = log2(or)) %>% 
  mutate(sig = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 & pval >= 0.001 ~ "**",
    pval < 0.05 & pval >= 0.01 ~ "*",
    T ~ ""
  ))

rownames(df) <- NULL
ind <- sort(ind)

p.list <- lapply(ind, function(x){
  df %>% 
    filter(trait == x) %>% 
    ggplot(aes(x = method, y = or)) +
    geom_col(aes(fill = method), width = 0.8) +
    scale_fill_manual(values = c("#0078E9", "#FFA500", "purple", "#C40000", "#700045")) +
    geom_hline(yintercept = 1, color = "white", linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_text(aes(y = or + 0.02, label = sig), size = 8) +
    #coord_fixed(ratio = 3.5) +
    ylab(paste(gsub("(.*)\\.+(.*?)", "\\2", x), "(OR)", sep = " ")) +
    ggtitle(label = gsub("(.*?)\\|(.*)", "\\1", x)) +
    xlab("")
})

p.list[[1]]

plot_grid(plotlist = p.list, ncol = 6)
ggsave(filename = "figure_2_outside_barplot.pdf" 
       ,width = 15, height = 10)

#-------------------------------------------------------------#
# figure 2 inner circos plot for 27 representative UVMR results
# no less than two cancer types p < 0.05
#-------------------------------------------------------------#
# circos input data preparing
# circos plot 

# cd circos/
# perl "circos-0.69/bin/circos" -config  "./fig2_circos/circos.config"






