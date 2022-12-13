####################
# dn_ds_diff.R     #
# Jessie Pelosi    #
# 8 Dec 2022       #
####################

library(dplyr)
library(ggplot2)
library(outliers)
library(EnvStats)

# Expression 

NDUV_scaffold_FPKM <- read.delim("NDUV_scaffold_FPKM", sep = " ")
NDUV_FPMK_dNdS <- merge(NDUV_scaffold_FPKM, dn_ds_all, by = 'V1')

NDUV_FPKM_pur <- NDUV_FPMK_dNdS %>% 
  filter(V2 < 1)

ggplot(data = NDUV_FPKM_pur, mapping = aes(x = FPKM, y = V2)) + geom_point() + 
  xlim(0,500) + geom_smooth()

# Revised sheet 
prms <- readxl::read_excel("../../../CodeMLParams/VittariaParams.xlsx")
prmsFilt <- prms %>% 
  filter(appdNdS < 1 & lindNdS < 1) %>% 
  mutate(dNdSDif = appdNdS - lindNdS)

mean(prmsFilt$dNdSDif)

mean(prmsFilt$appdNdS); mean(prmsFilt$appdNdS) + sd(prmsFilt$appdNdS); mean(prmsFilt$appdNdS) - sd(prmsFilt$appdNdS)
mean(prmsFilt$lindNdS); mean(prmsFilt$lindNdS) + sd(prmsFilt$lindNdS); mean(prmsFilt$lindNdS) - sd(prmsFilt$lindNdS)
mean(prmsFilt$appdN); mean(prmsFilt$appdN) + sd(prmsFilt$appdN); mean(prmsFilt$appdN) - sd(prmsFilt$appdN)
mean(prmsFilt$lindN); mean(prmsFilt$lindN) + sd(prmsFilt$lindN); mean(prmsFilt$lindN) - sd(prmsFilt$lindN)
mean(prmsFilt$appS); mean(prmsFilt$appS) + sd(prmsFilt$appS); mean(prmsFilt$appS) - sd(prmsFilt$appS)
mean(prmsFilt$linS); mean(prmsFilt$linS) + sd(prmsFilt$linS); mean(prmsFilt$linS) - sd(prmsFilt$linS)
mean(prmsFilt$appN); mean(prmsFilt$appN) + sd(prmsFilt$appN); mean(prmsFilt$appN) - sd(prmsFilt$appN)
mean(prmsFilt$linN); mean(prmsFilt$linN) + sd(prmsFilt$linN); mean(prmsFilt$linN) - sd(prmsFilt$linN)

wilcox.test(prmsFilt$appdS, prmsFilt$lindS, paired = T)
wilcox.test(prmsFilt$appdN, prmsFilt$lindN, paired = T)
wilcox.test(prmsFilt$appdNdS, prmsFilt$lindNdS, paired = T)

ggplot(data = prmsFilt, mapping = aes(x = dNdSDif)) + geom_histogram( bins = 250, fill = "gray23", alpha = 0.5) + 
  xlab(bquote(Delta~dN/dS [ Asex-Sex])) + ylab("Frequency") + xlim(-0.5,0.5) +
  theme_classic() + geom_vline(xintercept = mean(prmsFilt$dNdSDif), color = "red", alpha = 0.7) 

mean(prmsFilt$dNdSDif)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
} # use function to determine outlier candidates 

sub <- is_outlier(prmsFilt$dNdSDif) 

cand_outliers <- prmsFilt[sub,] #subset dataframe to select candidate outliers

length(cand_outliers$dNdSDif) #709 candidate outliers 

pos_out <- cand_outliers %>% 
  filter(dNdSDif > 0)

write.csv(file = "pos_out.csv", pos_out)

neg_out <- cand_outliers %>% 
  filter(dNdSDif < 0)

write.csv(file= "neg_out.csv", neg_out)


#### Roser Test #####

test <- rosnerTest(prmsFilt$dNdSDif, k=709)
summary(test)
test$all.stats

# 59 observations are outliers based on Rosner test 

dn_ds_pur_test <- prmsFilt

dn_ds_pur_test <- rename(dn_ds_pur_test, Value=dNdSDif)

all_outliers <- inner_join(dn_ds_pur_test, test$all.stats) %>% 
  filter(Outlier == "TRUE")

pos_outliers <- filter(all_outliers, Value > 0)
count(pos_outliers) #48
neg_outliers <- filter(all_outliers, Value < 0)
count(neg_outliers) #11

write.table(x=pos_outliers$Locus,file = "positive_outliers_dec8.txt", quote = F)
write.table(x=neg_outliers$Locus, file = "negative_outliers_dec8.txt", quote = F)
write.table(x=dn_ds_pur_test$Locus, file = "all_dn_ds_dec8.txt", quote = F)

## Concatenated Matrix ## 
cat <- read.csv("dNdSCat.csv", header = T, row.names = 1)
catdf <- as.data.frame(cat)
catNS <- select(catdf, N, S)

fisher.test(catNS, alternative = "greater")
