########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()

############################ Load R package and data ###########################
library(readxl)    #读取excel表
library(dplyr)
library(writexl)    #保存excel表
library(lme4)    #混合线性模型
library(lmerTest)    #为lmer模型拟合对象补充了固定效应的P值
library(broom.mixed)    #整理混合线性模型结果
library(tidyverse)    #unnest函数
library(ggplot2)
library(scales)
library(tidytext)   # 提供 reorder_within 和 scale_x_reordered
library(ggh4x)    #facet_wrap2,修改分面标签背景颜色
library(tidyr)
library(circlize)        # 提供colorRamp2函数
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplotify)    #把热图转成ggplot对象
library(cowplot)
library(Hmisc)    #用于rcorr 函数（计算相关系数和 p 值）
library(corrplot)     #用于可视化相关性矩阵


tlens <- read_xlsx("tlens_chr_filtered.xlsx")
clinical_info <- read_xlsx("31family_info.xlsx", range = "A1:G125")    #clinical phenotype data: family_ID, sample_ID, clinical_group, sex, age


################## Combine clinical information with ChArmTLs ##################
######## calculate the mTL of the individual and the average TL of the chromosome arm
#### calculate the average of two alleles by chromosome arm
tlens_dip_arm <- tlens %>%
  group_by(sample, chr) %>%
  summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "Diploid") %>%
  select(sample, level, chr, TL_p75)

#### mTL of the sample
tlens_dip_mTL <- tlens %>%
  group_by(sample) %>%
  summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
  mutate(chr = "mTL", level = "Diploid") %>%
  select(sample, level, chr, TL_p75)


######## Calculate the mTL of the offspring haplotype and extract the ChArmTLs
#### extract ChArmTLs
tlens_hap_arm <- tlens %>%
  filter(clinical_group %in% c("ASD", "Sibling")) %>%
  select(sample, parent_of_origin, chr, TL_p75) %>%
  rename(level = parent_of_origin)

#### mTL of paternal and maternal haplotypes
tlens_hap_mTL <- tlens %>%
  filter(clinical_group %in% c("ASD", "Sibling")) %>%
  group_by(sample, parent_of_origin) %>%
  summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
  rename(level = parent_of_origin) %>%
  mutate(chr = "mTL") %>%
  select(sample, level, chr, TL_p75)


######## merge data
tlens_mean <- bind_rows(tlens_dip_arm, tlens_dip_mTL, tlens_hap_arm, tlens_hap_mTL) %>%
  mutate(level = factor(level, levels = c("Diploid","Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(sample, level, chr)

tlens_mean <- merge(tlens_mean, clinical_info, by = "sample")     # combined with clinical information by the sample
tlens_mean <- tlens_mean[,c(5,1,6:10,2:4)]     # modify the order of the columns

write_xlsx(tlens_mean, "tlens_chr_filtered_mean.xlsx")


################################ Sex differences ################################
tlens_mean <- tlens_mean %>%
  mutate(diagnosis = if_else(clinical_group == "ASD", "ASD", "normal")) %>%    # add diagnostic status
  relocate(diagnosis, .after = clinical_group)

######## define a linear mixed-effects model, adjusting for age, disease status, and family ID
fit_sex_model <- function(data) {
  data <- data %>% filter(!is.na(TL_p75))    # Filter out missing values
  
  data <- data %>%    
    mutate(sex = factor(sex, levels = c("Male", "Female")),    # Ensure sex is a factor, with Male as the reference
           diagnosis = factor(diagnosis, levels = c("normal","ASD")))    # Ensure diagnosis is a factor, with normal as the reference
  
  model <- tryCatch({
    lmer(TL_p75 ~ sex + age + diagnosis + (1 | family), data = data)
  }, warning = function(w) {
    suppressWarnings(lmer(TL_p75 ~ sex + age + diagnosis + (1 | family), data = data))    # Catch warnings but continue execution
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(NULL)
  
  model_output <- tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%    # Extract the coefficient table for all fixed effects (sex, age, diagnosis)
    mutate(n_obs = nrow(data), singular_flag = isSingular(model)) %>%
    filter(term != "(Intercept)")     # Remove the intercept term
  
  return(model_output)
}

sex_diff <- tlens_mean %>%
  group_by(level, chr) %>%
  nest() %>%
  mutate(model_summary = map(data, fit_sex_model)) %>%
  select(-data) %>%
  unnest(cols = c(model_summary))

sex_diff <- sex_diff %>%
  group_by(level, term) %>%
  mutate(p_adj_fdr = case_when(
    chr == "mTL" ~ p.value,
    TRUE ~ p.adjust(p.value, method = "fdr")    # Adjust p-values using FDR method
  )) %>%
  ungroup() %>%
  mutate(level = factor(level, levels = c("Diploid","Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(level, chr) %>%
  select(level,chr,n_obs,term,estimate,std.error,conf.low,conf.high,statistic,p.value,p_adj_fdr,singular_flag,everything())     # Reorder columns for better readability

write_xlsx(sex_diff, "sex_diff.xlsx")


######################### Sex difference visualization #########################
sex_diff_dip_plot <- sex_diff %>%
  filter(level == "Diploid", term == "sexFemale") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         estimate_sign = ifelse(estimate > 0, "Female > Male", "Female < Male"),
         significance = ifelse(p.value < 0.05, "Significant", "Not Significant"),
         logP = ifelse(p.value < 0.05, -log10(p.value), 1))      # Effect direction

sex_diff_dip_plot$estimate_sign <- factor(sex_diff_dip_plot$estimate_sign, levels = c("Female > Male", "Female < Male"))
sex_diff_dip_plot$significance <- factor(sex_diff_dip_plot$significance, levels = c("Significant", "Not Significant"))

p1 <- ggplot(sex_diff_dip_plot, aes(x = chr_number, y = abs(estimate), color = estimate_sign, shape = significance, size = logP))+
  geom_point(stroke = 0.35)+    # Adjust the stroke thickness of the points
  scale_color_manual(values = c("#EB6368","#82ADD0"))+
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1))+
  theme_bw()+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  scale_y_continuous(limits = c(0,900), breaks = c(0,300,600,900))+
  scale_size_continuous(range = c(0.5,3), breaks = c(1, -log10(0.05), 2, 3, 4), labels = c("ns","0.05","0.01","0.001", "0.0001"))+    # Control the minimum and maximum sizes of the bubbles
  labs(x = NULL, y = "Mean difference (bp)", size = "P", color = NULL)+
  guides(color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1), 
         size = guide_legend(override.aes = list(shape = 1, stroke = 0.35), position = "top", order = 2),
         shape = "none")+    # Modify the scatter point size in the legend and the legend position
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),    # Distance between the legend box and the plot area
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


##################### Sex difference direction consistency #####################
# Perform sign test grouped by level
level <- unique(sex_diff$level)

sex_diff_direction <- list()
for (g in level) {
  sub <- sex_diff %>% filter(level == g, chr != "mTL", term == "sexFemale")
  pos <- sum(sub$estimate > 0)
  neg <- sum(sub$estimate < 0)
  total <- pos + neg
  # If estimate == 0 exists, it can be ignored or allocated as needed; assuming no 0 values here
  test <- binom.test(pos, total, p = 0.5, alternative = "greater")
  sex_diff_direction[[g]] <- data.frame(
    group = g,
    positive = pos,
    negative = neg,
    total = total,
    prop_positive = pos / total,
    p_value = test$p.value,
    lower_ci = test$conf.int[1],
    upper_ci = test$conf.int[2]
  )
}

rm(level, g, sub, pos, neg, total, test)    # Free up memory

# Combine results
sex_diff_direction <- do.call(rbind, sex_diff_direction)
write_xlsx(sex_diff_direction, "sex_diff_direction.xlsx")


################################ Age correlation ################################
age_cor <- sex_diff %>% filter(term == "age")    # Extract directly from the sex analysis results

write_xlsx(age_cor, "age_correlation.xlsx")


######################### Age correlation visualization #########################
age_cor_dip_plot <- age_cor %>%
  filter(level == "Diploid") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         logP = ifelse(p.value > 1e-10, -log10(p.value), -log10(p.value)-10))    # Note the p-value range, the p-value for mTL is vastly different from the chromosome arms

p2 <- ggplot(age_cor_dip_plot, aes(x = chr_number, y = estimate, size = logP, color = logP))+
  geom_point()+
  scale_color_gradientn(colors = c("lightgrey", "#EB6368", "darkred"), values = rescale(c(0, 1.23, 8.8)), name = "-log10(P)")+
  theme_bw()+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  scale_y_reverse(expand =c(0,0), limits = c(-20,-60), breaks = c(-30,-40,-50), labels = c(30,40,50))+
  scale_size_continuous(range = c(0.5,3))+    # Control the minimum and maximum sizes of the bubbles
  labs(x = NULL, y = "Annual rate of TL loss (bp/year)", size = NULL)+
  guides(size = "none")+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),    # Distance between the legend box and the plot area
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


####################### Conserved rank order of ChArmTLs #######################
tlens$chr <- factor(tlens$chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                          "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                          "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                          "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q"))


########sample
tlens$clinical_group <- factor(tlens$clinical_group, levels = c("Sibling","ASD","Father","Mother"))     # Convert to factor

sample_TL_stats <- tlens %>%    # Calculate the mean and sd for each sample
  group_by(sample) %>%
  summarise(sample_TL_mean = mean(TL_p75, na.rm = TRUE),
            sample_TL_sd = sd(TL_p75, na.rm = TRUE),
            .groups = "drop")

tlens_sample_rank <- tlens %>%    # Merge statistics back into original data, and calculate relative deviation and z-score
  left_join(sample_TL_stats, by = "sample") %>%
  mutate(relative_TL = TL_p75 - sample_TL_mean,
         z_score = relative_TL / sample_TL_sd) %>%
  select(-sample_TL_mean, -sample_TL_sd)

sample_arm_rank <- tlens_sample_rank %>%    # Group by clinical_group and chr, and calculate the mean of relative deviation and z-score
  group_by(clinical_group, chr) %>%
  summarise(mean_relative_TL = mean(relative_TL, na.rm = TRUE),
            mean_z_score = mean(z_score, na.rm = TRUE),
            sd_z_score = sd(z_score, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(clinical_group) %>%    # Calculate rank values in ascending order of the mean within each clinical_group
  mutate(rank_relative_TL = rank(mean_relative_TL, ties.method = "min"),
         rank_z_score = rank(mean_z_score, ties.method = "min"))

tlens_sample_rank <- tlens_sample_rank %>%    # Merge original data and statistics for plotting
  left_join(sample_arm_rank, by = c("clinical_group", "chr"))

tlens_sample_rank$chr_number <- substr(tlens_sample_rank$chr,4,nchar(as.character(tlens_sample_rank$chr)))    # Remove "chr" from chromosome arms
sample_arm_rank$chr_number <- substr(sample_arm_rank$chr,4,nchar(as.character(sample_arm_rank$chr)))

p3 <- ggplot()+
  geom_jitter(data = tlens_sample_rank, aes(x = reorder_within(chr_number, rank_z_score, clinical_group), y = z_score, color = clinical_group),
              width = 0.25, size = 0.75, show.legend = FALSE)+    # Jitter points, color by clinical_group (distinguishable within facets)
  geom_errorbar(data = sample_arm_rank, aes(x = reorder_within(chr_number, rank_z_score, clinical_group), ymin = mean_z_score - sd_z_score, ymax = mean_z_score + sd_z_score),
                width = 0.25, linewidth = 0.35)+    # Error bars: mean ± sd
  stat_summary(data = tlens_sample_rank, aes(x = reorder_within(chr_number, rank_z_score, clinical_group), y = z_score), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+    # Mean
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  facet_wrap2(~ clinical_group, nrow = 4, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#82ADD080","#EB636880","#FCD28c80","#A1B5A380"))))+    # Facet, free x-axis for each facet (but the order is already fixed by reorder_within)
  scale_color_manual(values = c("#82ADD0","#EB6368","#FCD28c","#A1B5A3"))+
  theme_bw()+
  scale_x_reordered()+    # Use scale_x_reordered to correctly display the original chr labels
  labs(x = NULL, y = "ChArmTL Z-score")+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 4.pdf", plot = p1, width = 17, height = 15, units = "cm", dpi = 600)


########offspring haplotypes
tlens_offspring <- tlens %>% filter(clinical_group %in% c("ASD", "Sibling"))    # Extract offspring groups

tlens_offspring <- tlens_offspring %>%    # Add phase_group and convert to factor
  mutate(phase_group = case_when(
    clinical_group == "ASD" & parent_of_origin == "Pat" ~ "ASD-Pat",
    clinical_group == "ASD" & parent_of_origin == "Mat" ~ "ASD-Mat",
    clinical_group == "Sibling" & parent_of_origin == "Pat" ~ "Sib-Pat",
    clinical_group == "Sibling" & parent_of_origin == "Mat" ~ "Sib-Mat",
    TRUE ~ NA_character_
  )) %>%
  mutate(phase_group = factor(phase_group, levels = c("Sib-Pat","Sib-Mat","ASD-Pat","ASD-Mat")))

hap_TL_stats <- tlens_offspring %>%    # Calculate the mean and sd for paternal (Pat) and maternal (Mat) within each offspring
  group_by(sample, parent_of_origin) %>%
  summarise(hap_TL_mean = mean(TL_p75, na.rm = TRUE),
            hap_TL_sd = sd(TL_p75, na.rm = TRUE),
            .groups = "drop")

tlens_hap_rank <- tlens_offspring %>%    # Merge statistics back into original data, and calculate relative deviation and z-score
  left_join(hap_TL_stats, by = c("sample", "parent_of_origin")) %>%
  mutate(relative_TL = TL_p75 - hap_TL_mean,
         z_score = relative_TL / hap_TL_sd) %>%
  select(-hap_TL_mean, -hap_TL_sd)

hap_arm_rank <- tlens_hap_rank %>%    # Group by phase_group and chr, and calculate the mean of relative deviation and z-score
  group_by(phase_group, chr) %>%
  summarise(mean_relative_TL = mean(relative_TL, na.rm = TRUE),
            mean_z_score = mean(z_score, na.rm = TRUE),
            sd_z_score = sd(z_score, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(phase_group) %>%    # Calculate rank values in ascending order of the mean within each phase_group
  mutate(rank_relative_TL = rank(mean_relative_TL, ties.method = "min"),
         rank_z_score = rank(mean_z_score, ties.method = "min"))

tlens_hap_rank <- tlens_hap_rank %>%    # Merge original data and statistics for plotting
  left_join(hap_arm_rank, by = c("phase_group", "chr"))

tlens_hap_rank$chr_number <- substr(tlens_hap_rank$chr,4,nchar(as.character(tlens_hap_rank$chr)))    # Remove "chr" from chromosome arms
hap_arm_rank$chr_number <- substr(hap_arm_rank$chr,4,nchar(as.character(hap_arm_rank$chr)))

p4 <- ggplot()+
  geom_jitter(data = tlens_hap_rank, aes(x = reorder_within(chr_number, rank_z_score, phase_group), y = z_score, color = phase_group),
              width = 0.25, size = 0.75, show.legend = FALSE)+    # Jitter points, color by phase_group (distinguishable within facets)
  geom_errorbar(data = hap_arm_rank, aes(x = reorder_within(chr_number, rank_z_score, phase_group), ymin = mean_z_score - sd_z_score, ymax = mean_z_score + sd_z_score),
                width = 0.25, linewidth = 0.35)+    # Error bars: mean ± sd
  stat_summary(data = tlens_hap_rank, aes(x = reorder_within(chr_number, rank_z_score, phase_group), y = z_score), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+    # Mean
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  facet_wrap2(~ phase_group, nrow = 4, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#82ADD0","#82ADD080","#EB6368","#EB636880"))))+    # Facet, free x-axis for each facet (but the order is already fixed by reorder_within)
  scale_color_manual(values = c("#82ADD0","#82ADD080","#EB6368","#EB636880"))+
  theme_bw()+
  scale_x_reordered()+    # Use scale_x_reordered to correctly display the original chr labels
  labs(x = NULL, y = "ChArmTL Z-score")+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 5.pdf", plot = p2, width = 17, height = 15, units = "cm", dpi = 600)


##################### Heatmap of ChArmTL rank conservation ######################
######## Combine sample and haplotype data
hap_arm_conservatism <- hap_arm_rank %>%
  select(chr, phase_group, rank_z_score) %>%
  pivot_wider(names_from = phase_group, values_from = rank_z_score)

sample_arm_conservatism <- sample_arm_rank %>%
  select(chr, clinical_group, rank_z_score) %>%
  pivot_wider(names_from = clinical_group, values_from = rank_z_score)

arm_conservatism <- sample_arm_conservatism %>%
  full_join(hap_arm_conservatism, by = "chr")


######## Add two datasets from Nature
arm_conservatism$`Cord blood`<- c(26,NA,32,6,36,13,20,34,31,21,12,29,9,NA,33,10,14,5,28,22,18,24,3,35,30,
                                  25,19,8,1,15,NA,27,11,23,17,2,7,4)
arm_conservatism$`147 adults`<- c(31,23,35,8,38,6,25,36,32,12,19,29,10,33,30,18,15,7,22,28,21,20,3,37,
                                  34,24,16,13,1,4,NA,26,9,27,11,2,14,5)

arm_conservatism$rank <- apply(arm_conservatism[,c(2:11)], 1, mean, na.rm = T)    # Calculate by row, used for sorting (mean rank value)

rank_order <- as.vector(arm_conservatism[order(arm_conservatism$rank),]$chr)   # Sort from shortest to longest
rank_order <- sub("chr", "", rank_order)    # Chromosome arm order


######## Convert data to matrix
arm_conservatism_matrix <- as.matrix(arm_conservatism[,c(2:11)])    # Convert to matrix format
rownames(arm_conservatism_matrix) <- gsub("^chr", "", arm_conservatism$chr)     # Use chr as row names
arm_conservatism_matrix_t <- t(arm_conservatism_matrix)    # Transpose
arm_conservatism_matrix_t <- arm_conservatism_matrix_t[,rank_order]    # Rearrange chromosome arm order from shortest to longest


######## Define heatmap colors
col_fun <- colorRamp2(seq(min(arm_conservatism_matrix_t, na.rm = TRUE), 
                          max(arm_conservatism_matrix_t, na.rm = TRUE), length.out = 7),
                      rev(brewer.pal(7, "RdYlBu")))


######## Create row split factor (groups of 4 rows, last group with 2 rows)
row_split <- factor(rep(c("Group1", "Group2", "Group3"), times = c(4, 4, 2)))

p5 <- Heatmap(arm_conservatism_matrix_t,
              column_title = "Rank order of ChArmTL",    # Column title
              column_title_side = "top",    # Column title position
              column_title_gp = gpar(fontsize = 8, fontface = "bold"),    # Column title font size, bold
              cluster_rows = FALSE, cluster_columns = FALSE,    # Do not cluster rows, do not cluster columns
              col = col_fun, na_col = "grey50",    # Color mapping, NA value color
              show_row_names = TRUE,            # Show row names (sample types)
              show_column_names = TRUE,         # Show column names (chromosome arms)
              column_names_side = "bottom",     # Place column names at the bottom
              row_names_side = "left",          # Place row names on the left
              column_names_gp = gpar(fontsize = 6), column_names_rot = 70,  # Adjust column name font and angle
              row_names_gp = gpar(fontsize = 6),     # Adjust row name font
              row_split = row_split,            # Split by row
              row_title = NULL,    # Hide row split title
              heatmap_legend_param = list(
                title = "Rank",
                title_position = "topleft",
                title_gp = gpar(fontsize = 7),       # Legend title font size
                labels_gp = gpar(fontsize = 6),       # Legend label font size
                grid_width = unit(2, "mm"))      # Adjust the thickness of the color bar
)

p5_1 <- as.ggplot(p5) + theme(plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1, unit = "cm"))


######## Combine plots
p1_2 <- plot_grid(p1, p2, nrow = 2, labels = c('a', 'b'), label_size = 8, rel_heights = c(1,0.9), align = 'v', axis = 'lr')
p5_1_NULL <- plot_grid(p5_1, NULL, nrow = 1, labels = c('c', ''), label_size = 8, rel_widths = c(1,0.05), align = 'h', axis = 'tb')

p1_2_3 <- plot_grid(p1_2, p5_1_NULL, ncol = 1, rel_heights = c(1,0.6), align = 'v')

ggsave(filename = "Figure 1.pdf", plot = p1_2_3, width = 17, height = 17, units = "cm", dpi = 600)


############### Correlation of ChArmTL rank order across groups ################
arm_conservatism_cor <- rcorr(arm_conservatism_matrix, type = "spearman")    # Correlation matrix

write.csv(arm_conservatism_cor$r, "arm_conservatism_cor.csv")
write.csv(arm_conservatism_cor$P, "arm_conservatism_cor_P.csv")

pdf("Supplementary Figure 6.pdf", width = 5, height = 5)
corrplot(arm_conservatism_cor$r,
         method = "color",
         col = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
         type = "upper", 
         addCoef.col ="black",
         number.cex =0.6,
         tl.cex = 0.7,
         tl.col = "black",
         tl.srt = 45)
dev.off()