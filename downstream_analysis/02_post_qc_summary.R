########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
install.packages(c("readxl", "dplyr", "ggplot2", "ggh4x", 
                   "writexl", "cowplot", "rstatix", "car", 
                   "ggsignif", "patchwork"))

library(readxl)
library(dplyr)
library(ggplot2)
library(ggh4x)    # facet_wrap2, modify facet strip background colors
library(writexl)
library(cowplot)
library(rstatix)    # shapiro_test() function
library(car)    # Levene test
library(ggsignif)    #Add significance
library(patchwork)


tlens <- read_xlsx("tlens_unphased_filtered.xlsx")    # The result file of ChArmTLs with unphased data filtered out.  8263 rows


############### Statistics and visualization of supporting reads ################
########Mean values and quantiles
mean(tlens$read_number)    #6.030134

quantile(tlens$read_number, .25)    #4
median(tlens$read_number)    #6
quantile(tlens$read_number, .75)    #8

tlens %>%
  group_by(clinical_group) %>%
  summarise(
    n = n(),
    mean = mean(read_number),
    Q1 = quantile(read_number, .25),
    median = median(read_number),
    Q3 = quantile(read_number, .75)
  )

tlens$clinical_group <- factor(tlens$clinical_group, levels = c("Sibling","ASD","Father","Mother"))

p1 <- ggplot(tlens, aes(x = chr, y = read_number, color = clinical_group))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, outlier.alpha = 0, show.legend = FALSE)+
  facet_wrap2(clinical_group ~ ., nrow = 4, strip.position = "right",
             strip = strip_themed(background_y = elem_list_rect(fill = c("#82ADD080","#EB636880","#FCD28c80","#A1B5A380"))))+
  scale_color_manual(values = c("#82ADD0","#EB6368","#FCD28c","#A1B5A3"))+
  theme_bw()+
  scale_x_discrete(limits = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                              "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                              "chr11p","chr11q","chr12p","chr12q","chr13p","chr13q","chr14p","chr14q","chr15p","chr15q",
                              "chr16p","chr16q","chr17p","chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q",
                              "chr21p","chr21q","chr22p","chr22q"),
                   labels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q", "8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q","19p","19q",
                              "20p","20q","21p","21q","22p","22q"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,20), breaks = c(5,10,15,20))+
  labs(x = NULL, y = "Number of supporting reads")+
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

ggsave(filename = "Supplementary Figure 1.pdf", plot = p1, width = 17, height = 15, units = "cm", dpi = 600)


######################## Haplotype coverage of ChArmTLs ########################
arm_coverage <- as.data.frame(table(tlens$chr, tlens$clinical_group))
colnames(arm_coverage) <- c("chr", "clinical_group", "counts")    #Rename the default column names
arm_coverage$coverage <- arm_coverage$counts/248     #124 samples * 2 haplotypes
arm_coverage$chr_number <- substr(arm_coverage$chr, 4, nchar(as.character(arm_coverage$chr)))    #Sort by coverage descending when plotting


p2_top <- ggplot(arm_coverage, aes(x = reorder(chr_number, -coverage), y = coverage, fill = clinical_group))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#82ADD0","#EB6368","#FCD28c","#A1B5A3"))+
  theme_bw()+
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0.25,0.5,0.75,1),
                     sec.axis = sec_axis(transform = ~ .*248, name = "Number of alleles", breaks = c(50,100,150,200)))+
  labs(x = NULL, y = "Coverage (%)", fill = NULL)+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.y.right = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),    #Distance between the legend box and the plot area
        legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


###################### Allele count statistics for chr22p ######################
telogator2_4 <- read_xlsx("tlens_outlier_filtered.xlsx")

chr22p_count <- subset(telogator2_4, chr == "chr22p")    #Extract chr22p
chr22p_count <- as.data.frame(table(chr22p_count$sample))    #Count the number of chr22p alleles per sample
chr22p_count$Freq <- as.character(chr22p_count$Freq)    #Convert to character

p2_bottom <- ggplot(chr22p_count, aes(x = Freq))+
  geom_bar(width = 0.7, fill = "#82ADD0")+
  theme_bw()+
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = c(10,20,30))+
  labs(x = "Number of ChArmTLs on chr22p", y = "Number of samples")+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

p2_bottom_null <- plot_grid(p2_bottom, NULL, rel_widths = c(1,2))

p2 <- plot_grid(p2_top, p2_bottom_null, ncol = 1, rel_heights = c(1,0.9))
ggsave(filename = "Supplementary Figure 2.pdf", plot = p2, width = 17, height = 9, units = "cm", dpi = 600)


##################### Count telomere-length observations per sample and offspring subtype #####################
tlens_chr_filtered <- subset(tlens, chr != "chr13p" & chr != "chr14p" & chr != "chr15p" & chr != "chr21p" & chr != "chr22p" & chr != "chr15q")    #8263 rows ---> 7794 rows

write_xlsx(tlens_chr_filtered, "tlens_chr_filtered.xlsx")


########Per sample
alleles_count <- as.data.frame(table(tlens_chr_filtered$sample))    # Count telomere-length observations per sample
colnames(alleles_count) <- c("sample", "counts")    #Rename the default column names
alleles_count$clinical_group <- substr(alleles_count$sample,
                                                nchar(as.character(alleles_count$sample))-1,
                                                nchar(as.character(alleles_count$sample)))    #Add the clinical group
alleles_count$clinical_group <- factor(alleles_count$clinical_group, levels = c("s1", "p1", "fa", "mo"))

####Mean values and median
mean(alleles_count$counts)    #62.85484
median(alleles_count$counts)    #63.5

alleles_count %>%
  group_by(clinical_group) %>%
  summarise(
    n = n(),
    mean = mean(counts),
    median = median(counts)
  )  

####Test whether there are differences between groups
alleles_count %>%
  group_by(clinical_group) %>%
  shapiro_test(counts)    # p1 and fa are non-normal distributions
leveneTest(counts ~ clinical_group, data = alleles_count)    #Equal variances
kruskal.test(counts ~ clinical_group, data = alleles_count)    #Test whether there are differences between groups
dunn_test(alleles_count, counts ~ clinical_group, p.adjust.method = "bonferroni")    # Post hoc test

p3 <- ggplot(data = alleles_count, aes(x = clinical_group, y = counts, fill = clinical_group))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#82ADD0","#EB6368","#FCD28c","#A1B5A3"))+
  theme_classic()+
  scale_x_discrete(labels = c("Sibling","ASD","Father","Mother"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(50,75), breaks = c(55,60,65,70))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))+
  geom_signif(y_position = c(70,73), xmin = c(1,3), xmax = c(2,4), annotations = c("ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p3)$data[[1]] # Check the upper whisker endpoint (ymax) to determine geom_signif y_position


########Offspring haplotypes
alleles_count_offsring <- as.data.frame(table(tlens_chr_filtered$sample, tlens_chr_filtered$parent_of_origin))    #Count ChArmTLs observations for each haplotype in offspring
colnames(alleles_count_offsring) <- c("sample", "haplotype","counts")    #Rename the default column names
alleles_count_offsring <- subset(alleles_count_offsring, counts != 0)    #Remove zero results
alleles_count_offsring$clinical_group <- substr(alleles_count_offsring$sample,
                                                nchar(as.character(alleles_count_offsring$sample))-1,
                                                nchar(as.character(alleles_count_offsring$sample)))    #Add the clinical group
alleles_count_offsring$clinical_group[alleles_count_offsring$clinical_group == "p1"] <- "ASD"    #Rename p1 to ASD
alleles_count_offsring$clinical_group[alleles_count_offsring$clinical_group == "s1"] <- "Sib"    #Rename s1 to Sib
alleles_count_offsring$phase_group <- paste(alleles_count_offsring$clinical_group, alleles_count_offsring$haplotype, sep = "-")    #Add phase_group column
alleles_count_offsring$phase_group <- factor(alleles_count_offsring$phase_group, levels = c("Sib-Pat","Sib-Mat","ASD-Pat","ASD-Mat"))

####Mean values and median
mean(alleles_count_offsring$counts)    #30.83065
median(alleles_count_offsring$counts)    #31

alleles_count_offsring %>%
  group_by(phase_group) %>%
  summarise(
    n = n(),
    mean = mean(counts),
    median = median(counts)
  )  

####Test whether there are differences between groups
alleles_count_offsring %>%
  group_by(phase_group) %>%
  shapiro_test(counts)    # ASD-Pat and ASD-Mat are non-normal distributions
leveneTest(counts ~ phase_group, data = alleles_count_offsring)    #Equal variances
kruskal.test(counts ~ phase_group, data = alleles_count_offsring)    #Test whether there are differences between groups
dunn_test(alleles_count_offsring, counts ~ phase_group, p.adjust.method = "bonferroni")    # Post hoc test

p4 <- ggplot(data = alleles_count_offsring, aes(x = phase_group, y = counts, fill = phase_group))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#82ADD0","#82ADD080","#EB6368","#EB636880"))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0) ,limits = c(23,40), breaks = c(25,30,35))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 20),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))+
  geom_signif(y_position = c(36,37,37.5,38.5), xmin = c(1,3,1,2), xmax = c(2,4,3,4), annotations = c("ns","ns","ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p4)$data[[1]] # Check the upper whisker endpoint (ymax) to determine geom_signif y_position

p3_4 <- p3 + p4 + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))

ggsave(filename = "Supplementary Figure 3.pdf", plot = p4_5, width = 10, height = 5, units = "cm", dpi = 600)

