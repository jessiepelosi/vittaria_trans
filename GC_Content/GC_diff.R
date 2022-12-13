####################################
####       GC DIFFERENCES       ####
####        Jessie Pelosi       ####
#### Last modified Nov 29, 2022 ####
####################################

library(dplyr)
library(ggplot2)
library(BSDA)

#### Compare observed GC content (8334 loci) ####

GC <- read.table("GC_dif.txt")

hist(GC$V2, breaks = 100, main = "Vittaria appalachiana total GC Content")
mean(GC$V2); min(GC$V2); max(GC$V2)
hist(GC$V3, breaks = 100, main = "Vittaria lineata total GC content")
mean(GC$V3); min(GC$V3); max(GC$V3)

ggplot() + geom_boxplot(data = GC, mapping = aes(x = V2), fill = "lightskyblue1") + xlim(20,85) + theme_classic()
ggsave("V.appalachiana_totalGC.png", height = 2, width = 7)
ggplot() + geom_boxplot(data = GC, mapping = aes(x = V3), fill = "indianred2") + xlim(20, 85) + theme_classic()
ggsave("V.lineata_totalGC.png", height = 2, width = 7)

ggplot(data = GC) + geom_histogram(mapping = aes(x = V4), bins = 250, fill = "gray23", alpha = 0.5) +
  xlab(bquote(Delta~'%GC'~ Content [ Asex-Sex])) + ylab("Frequency") +
  theme_classic() + geom_vline(xintercept = mean(GC$V4), color = "red", alpha = 0.7)

ggsave("Total_GC_content_difference_12-18-20.pdf", width = 7, height = 7)

mean(GC$V4) #0.282616

t.test(x= GC$V2, y = GC$V3,paired = T)

# Paired t-test
# data:  GC$V2 and GC$V3
# t= 13.429, df = 8333, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.2413635 0.3238685
# sample estimates:
#  mean of the differences 
# 0.282616

wilcox.test(x = GC$V2, y = GC$V3, paired = T)

# Wilcoxon signed rank test with continuity correction
# data:  GC$V2 and GC$V3
# V = 20307418, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

SIGN.test(GC$V4, mean = 0, alternative = "greater")
# One-sample Sign-Test
# data:  GC$V4
# s = 4723, p-value = 3.331e-16
# alternative hypothesis: true median is greater than 0
# 95 percent confidence interval:
#  0.1788983         Inf
# sample estimates:
#  median of x 
# 0.2006269

## GC Content at the third codon position

hist(GC$V5, breaks = 100, main = "Vittaria appalachiana GC3 Content")
mean(GC$V5); min(GC$V5); max(GC$V5)
hist(GC$V6, breaks = 100, main = "Vittaria lineata total GC3 content")
mean(GC$V6); min(GC$V6); max(GC$V6)

ggplot() + geom_boxplot(data = GC, mapping = aes(x = V5), fill = "lightskyblue1") + xlim(20,85) + theme_classic()
ggsave("V.appalachiana_GC3.png", height = 2, width = 7)
ggplot() + geom_boxplot(data = GC, mapping = aes(x = V6), fill = "indianred2") + xlim(20, 85) + theme_classic()
ggsave("V.lineata_GC3.png", height = 2, width = 7)

ggplot(data = GC) + geom_histogram(mapping = aes(x = V7), bins = 250, fill = "gray23", alpha = 0.5) +
  xlab(bquote(Delta~'%GC'~ Content [ Asex-Sex])) + ylab("Frequency") +
  theme_classic() + geom_vline(xintercept = mean(GC$V7), color = "red", alpha = 0.7)

ggsave("GC3_content_difference_12-18-20.pdf", width = 7, height = 7)

mean(GC$V7) #0.4078634

t.test(x= GC$V5, y = GC$V6,paired = T)

#  Paired t-test
#  data:  GC$V5 and GC$V6
#  t = 10.222, df = 8333, p-value < 2.2e-16
#  alternative hypothesis: true difference in means is not equal to 0
#  95 percent confidence interval:
#  0.3296510 0.4860758
#  sample estimates:
#  mean of the differences 
#  0.4078634 

wilcox.test(x = GC$V5, y = GC$V6, paired = T)
#  Wilcoxon signed rank test with continuity correction
#  data:  GC$V5 and GC$V6
#  V = 19468521, p-value < 2.2e-16
#  alternative hypothesis: true location shift is not equal to 0

SIGN.test(GC$V7, mean = 0, alternative = "greater")
#  One-sample Sign-Test
#  data:  GC$V7
#  s = 4566, p-value = 3.331e-16
#  alternative hypothesis: true median is greater than 0
#  95 percent confidence interval:
#  0.2806429       Inf
#  sample estimates:
#  median of x 
#  0.3325346 

#### Compare equilibrium GC content (Vittaria and Adiantum) #####

EquiGC <- read.delim("../../../TotalGC*.txt",header = T)

GC <- mean(EquiGC$GC.V.appalachiana); GC
GC + sd(EquiGC$GC.V.appalachiana); GC-sd(EquiGC$GC.V.appalachiana)

GC <- mean(EquiGC$GC.V.lineata); GC
GC + sd(EquiGC$GC.V.lineata); GC-sd(EquiGC$GC.V.lineata)

GC3 <- mean(EquiGC$GC3.V.appalachiana); GC3
GC3 + sd(EquiGC$GC3.V.appalachiana); GC3-sd(EquiGC$GC3.V.appalachiana)

GC3 <- mean(EquiGC$GC3.V.lineata); GC3
GC3 + sd(EquiGC$GC3.V.lineata); GC3-sd(EquiGC$GC3.V.lineata)

# Filter to remove inaccurate esimates (less than 1% and greater than 99%)

EquiGC_filt <- EquiGC %>% 
  filter(GC..V.appalachiana > 1 & GC..V.appalachiana < 99 &
           GC..V..lineata < 99 & GC..V..lineata > 1) %>% 
  select(GC..V.appalachiana, GC..V..lineata) %>% 
  mutate(GCEqDiff = GC..V.appalachiana - GC..V..lineata)

EquiGC_filt2 <- EquiGC %>% 
  filter(GC3..V.appalachiana > 1 & GC3..V.appalachiana < 99 &
           GC3..V..lineata < 99 & GC3..V..lineata > 1) %>% 
  select(GC3..V.appalachiana, GC3..V..lineata) %>% 
  mutate(GC3EqDiff = GC3..V.appalachiana - GC3..V..lineata)

GCE3 <- mean(EquiGC_filt2$GC3..V.appalachiana); GCE3
GCE3 + sd(EquiGC_filt2$GC3..V.appalachiana); GCE3-sd(EquiGC_filt2$GC3..V.appalachiana)

GCE3 <- mean(EquiGC_filt2$GC3..V..lineata); GCE3
GCE3 + sd(EquiGC_filt2$GC3..V..lineata); GCE3-sd(EquiGC_filt2$GC3..V..lineata)
  
hist(EquiGC_filt$GC..V.appalachiana, breaks = 100, main = "Vittaria appalachiana total GC Content")
GCE <- mean(EquiGC_filt$GC..V.appalachiana); GCE
GCE + sd(EquiGC_filt$GC..V.appalachiana); GCE-sd(EquiGC_filt$GC..V.appalachiana)
hist(EquiGC_filt$GC..V..lineata, breaks = 100, main = "Vittaria lineata total GC content")
mean(EquiGC_filt$GC..V..lineata); min(EquiGC_filt$GC..V..lineata); max(EquiGC_filt$GC..V..lineata)
GCE <- mean(EquiGC_filt$GC..V..lineata); GCE
GCE + sd(EquiGC_filt$GC..V..lineata); GCE-sd(EquiGC_filt$GC..V..lineata)

ggplot() + geom_boxplot(data = EquiGC_filt, mapping = aes(x = GC..V.appalachiana), fill = "lightskyblue1") + xlim(10,90) + theme_classic()
ggsave("V.appalachiana_totalGC.png", height = 2, width = 7)
ggplot() + geom_boxplot(data = EquiGC_filt, mapping = aes(x = GC..V..lineata), fill = "indianred2") + xlim(20, 85) + theme_classic()
ggsave("V.lineata_totalGC.png", height = 2, width = 7)

ggplot(data = EquiGC_filt) + geom_histogram(mapping = aes(x = GCEqDiff), bins = 250, fill = "gray23", alpha = 0.5) +
  xlab(bquote(Delta~'%GC'~ Content [ Asex-Sex])) + ylab("Frequency") +
  theme_classic() + geom_vline(xintercept = mean(EquiGC_filt$GCEqDiff), color = "red", alpha = 0.7)

ggsave("Total_GC*_content_difference_12-5-22.pdf", width = 7, height = 7)

wilcox.test(x = EquiGC$GC.V.appalachiana, y = EquiGC$GC.V.lineata, paired = T, alternative = "less")
#Wilcoxon signed rank test with continuity correction
#data:  EquiGC$GC.V.appalachiana and EquiGC$GC.V.lineata
#V = 6051743, p-value = 4.409e-09
#alternative hypothesis: true location shift is less than 0
mean(EquiGC$TotalGCDiff)

wilcox.test(x = EquiGC$GC3.V.appalachiana, y = EquiGC$GC3.V.lineata, paired = T, alternative = "less")
#Wilcoxon signed rank test with continuity correction
#data:  EquiGC$GC3.V.appalachiana and EquiGC$GC3.V.lineata
#V = 6150604, p-value = 8.598e-07
#alternative hypothesis: true location shift is less than 0
mean(EquiGC$GC3.Diff)

wilcox.test(x = EquiGC_filt$GC..V.appalachiana, y = EquiGC_filt$GC..V..lineata, paired = T, alternative = "less")
#Wilcoxon signed rank test with continuity correction
#data:  EquiGC_filt$V1 and EquiGC_filt$V2
#V = 2901200, p-value = 0.04283
#alternative hypothesis: true location shift is less than 0
mean(EquiGC_filt$GCEqDiff)

wilcox.test(x = EquiGC_filt2$GC3..V.appalachiana, y = EquiGC_filt2$GC3..V..lineata, paired = T, alternative = "less")
#Wilcoxon signed rank test with continuity correction
#data:  EquiGC_filt2$GC3..V.appalachiana and EquiGC_filt2$GC3..V..lineata
#V = 1238223, p-value = 0.004728
#alternative hypothesis: true location shift is less than 0
mean(EquiGC_filt2$GC3EqDiff)

