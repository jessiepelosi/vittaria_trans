#### GC DIFFERENCES ####
library(dplyr)
library(ggplot2)
library(BSDA)

GC <- read.table("GC_dif.txt")

hist(GC$V2, breaks = 100, main = "Vittaria appalachiana total GC Content")
hist(GC$V3, breaks = 100, main = "Vittaria lineata total GC content")

ggplot(data = GC) + geom_histogram(mapping = aes(x = V4), bins = 250, color = "gray23", alpha = 0.5) +
  xlab(bquote(Delta~'%GC'~ Content [ Asex-Sex])) + ylab("Frequency") +
  theme_classic() + geom_vline(xintercept = mean(GC$V4), color = "red", alpha = 0.5)

ggsave("Total_GC_content_difference_12-3-20.png", dpi = 300, width = 5, height = 4)

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
hist(GC$V6, breaks = 100, main = "Vittaria lineata total GC3 content")

ggplot(data = GC) + geom_histogram(mapping = aes(x = V7), bins = 250, color = "gray23", alpha = 0.5) +
  xlab(bquote(Delta~'%GC'~ Content [ Asex-Sex])) + ylab("Frequency") +
  theme_classic() + geom_vline(xintercept = mean(GC$V7), color = "red", alpha = 0.5)

ggsave("GC3_content_difference_12-3-20.png", dpi = 300, width = 5, height = 4)

mean(GC$V7)

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
