########################## 清空工作空间，设置工作路径 ##########################
rm(list=ls())
setwd("D:\\科研项目\\ASD三代端粒\\3_SCI论文\\第8版\\数据03版\\1_telogator2结果文件处理")
getwd()


#################################### 安装包 ####################################
install.packages("ggpubr")    #安装ggpubr包
install.packages("writexl")    #安装writexl包
install.packages("plyr")    #安装plyr包

library(readxl)    #读取excel
library(ggpubr)    #绘制进阶版Q-Q图
library(ggplot2)
library(writexl)    #保存excel


################################### 载入数据 ###################################
dir()    #查看目录文件

telogator2_3 <- read_excel("3_排除read数小于2.xlsx")    #9599条记录

summary(telogator2_3$TL_p75)    #端粒长度描述分析

hist(telogator2_3$TL_p75, prob = T)    #直方图
lines(density(na.omit(telogator2_3$TL_p75)), col = 'blue')   #密度曲线

ggqqplot(telogator2_3$TL_p75, color = "blue", main="Normal Q-Q Plot")    #Q-Q图


################################## 排除异常值 ##################################
########排除异常值(四分位距)，区分子代和亲代
telogator2_4 <- telogator2_3
telogator2_4$QC_status <- "passed"    #9599条通过
# telogator2_4 <- telogator2_4[,c(1:7,14,8:13)]    #调整QC_status列的位置

telogator2_4$QC_status[telogator2_4$TL_p75 <= 100] <- "outlier"    #107个outlier

# telogator2_4$clinical_group <- substr(telogator2_4$sample,nchar(telogator2_4$sample)-1,nchar(telogator2_4$sample))    #家系中角色
# telogator2_4 <- telogator2_4[,c(1,15,2:14)]    #调整clinical_group列的位置

####offspring
telogator2_4_offspring <- subset(telogator2_4, clinical_group == "p1" | clinical_group == "s1")    #4741条记录

summary(telogator2_4_offspring$TL_p75[telogator2_4_offspring$QC_status == "passed"])    #4702条记录 + 39 outlier
4937 - 1.5 * (7706 - 4937)    #783.5，下限
7706 + 1.5 * (7706 - 4937)    #11859.5，上限

telogator2_4$QC_status[telogator2_4$clinical_group == "p1" & telogator2_4$TL_p75 < 783.5] <- "outlier"    #18 outlier + 27 passed
telogator2_4$QC_status[telogator2_4$clinical_group == "s1" & telogator2_4$TL_p75 < 783.5] <- "outlier"    #21 outlier + 20 passed
telogator2_4$QC_status[telogator2_4$clinical_group == "p1" & telogator2_4$TL_p75 > 11859.5] <- "outlier"    #17 passed
telogator2_4$QC_status[telogator2_4$clinical_group == "s1" & telogator2_4$TL_p75 > 11859.5] <- "outlier"    #24 passed

table(telogator2_4$QC_status)    #195个outlier

####parent
telogator2_4_parent <- subset(telogator2_4, clinical_group == "fa" | clinical_group == "mo")

summary(telogator2_4_parent$TL_p75[telogator2_4_parent$QC_status == "passed"])    #4790条记录 + 68 outlier
3820 - 1.5 * (6360 - 3820)    #10，下限
6360 + 1.5 * (6360 - 3820)    #10170，下限

telogator2_4$QC_status[telogator2_4$clinical_group == "fa" & telogator2_4$TL_p75 < 10] <- "outlier"
telogator2_4$QC_status[telogator2_4$clinical_group == "mo" & telogator2_4$TL_p75 < 10] <- "outlier"
telogator2_4$QC_status[telogator2_4$clinical_group == "fa" & telogator2_4$TL_p75 > 10170] <- "outlier"    #17 passed
telogator2_4$QC_status[telogator2_4$clinical_group == "mo" & telogator2_4$TL_p75 > 10170] <- "outlier"    #30 passed

table(telogator2_4$QC_status)    #242(107 + 88 + 47)个outlier

write_xlsx(telogator2_4, "4_排除异常值(R输出).xlsx")


########排除异常值后的分布与正态性检验
telogator2_4_passed <- subset(telogator2_4, QC_status == "passed")    #9357条记录通过QC

summary(telogator2_4_passed$TL_p75)

hist(telogator2_4_passed$TL_p75, prob = T)  #直方图
lines(density(na.omit(telogator2_4_passed$TL_p75)), col = 'blue')   #密度曲线

ggqqplot(telogator2_4_passed$TL_p75, color = "blue", main="Normal Q-Q Plot")    #近似正态分布

