### DN/DS RATIOS ###
library(dplyr)
library(ggplot2)
library(outliers)
library(EnvStats)

dn_ds_all <- read.delim("../dn_ds_ratios.txt",header = F)

dn_ds_purifying <- dn_ds_all %>% 
  filter(V2 < 1 & V3 < 1) %>% 
  mutate(V4 = V2-V3)

ggplot(data = dn_ds_purifying, mapping = aes(x = V4)) + geom_histogram( bins = 250, color = "gray23", alpha = 0.5) + 
  xlab(bquote(Delta~dn/ds [ Asex-Sex])) + ylab("Frequency") + xlim(-0.2,0.2) +
  theme_classic() + geom_vline(xintercept = mean(dn_ds_purifying$V4), color = "red", alpha = 0.5) 

ggsave("pairwise_dn_ds.png", dpi = 300, height = 4, width = 5)

mean(dn_ds_purifying$V4)

t.test(dn_ds_purifying$V2, dn_ds_purifying$V3, paired = T, alternative = "greater")

# Paired t-test
# data:  dn_ds_purifying$V2 and dn_ds_purifying$V3
# t = 1.7523, df = 4968, p-value = 0.03989
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#  4.868521e-05          Inf
# sample estimates:
#  mean of the differences 
# 0.0007959147

wilcox.test(dn_ds_purifying$V2, dn_ds_purifying$V3,paired = T)

# Wilcoxon signed rank test with continuity correction
# data:  dn_ds_purifying$V2 and dn_ds_purifying$V3
# V = 6479921, p-value = 0.0002455
# alternative hypothesis: true location shift is greater than 0


## Get outer observations in pairwise comparisons

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
} # use function to determine outlier candidates 

sub <- is_outlier(dn_ds_purifying$V4) 

cand_outliers <- dn_ds_purifying[sub,] #subset dataframe to select candidate outliers

length(cand_outliers$V4) #502 candidate outliers 

test <- rosnerTest(dn_ds_purifying$V4, k=502)
summary(test)
test$all.stats

# 99 observations are outliers based on Rosner test 

dn_ds_pur_test <- dn_ds_purifying

dn_ds_pur_test <- rename(dn_ds_pur_test, Value=V4)

all_outliers <- inner_join(dn_ds_pur_test, test$all.stats) %>% 
  filter(Outlier == "TRUE")

pos_outliers <- filter(all_outliers, Value > 0)

neg_outliers <- filter(all_outliers, Value < 0)

write.table(x=pos_outliers$V1,file = "positive_outliers.txt")
write.table(x=neg_outliers$V1, file = "negative_outliers.txt")
write.table(x=dn_ds_pur_test$V1, file = "all_dn_ds.txt")

# DN/DS Hypothesis testing

likelihoods <- read.delim("../likelihoods.txt", header =T)

likelihoods$np <- rep(c(6,7,8,9,9), length= length(likelihoods$Ortholog))
likelihoods$df <- rep(c(0,1,2,3,3), length= length(likelihoods$Ortholog))

likelihoods <- likelihoods %>% 
  group_by(Ortholog) %>% 
  mutate(X2 = -2*(Likelihood[1]-Likelihood[row_number()])) %>% 
  mutate(pval = dchisq(X2[row_number()],df[row_number()])) %>% 
  mutate(AIC = -2*(Likelihood[row_number()])+2*(np[row_number()]))

likelihoods_significant <- likelihoods %>% 
  filter(pval < 0.05 & X2 != 0 & X2 > 1)

ortholog_summary <- likelihoods_significant %>% 
  summarise(count = n()) 

4982-length(ortholog_summary$Ortholog)
#3232 orthologs fit null hypothesis 

orthologs_single_hyp <- merge(ortholog_summary, likelihoods_significant, by = 'Ortholog')
orthologs_single_hyp <- orthologs_single_hyp %>% 
  select(Ortholog, count, Hypothesis) %>% 
  filter(count == 1) %>% 
  group_by(Hypothesis) %>% 
  summarise(no.hyp = n())

likelihoods_significant_AIC <- likelihoods_significant %>% 
  group_by(Ortholog) %>% 
  arrange(AIC) %>% 
  filter(row_number() == 1)

ortholog_summary <- likelihoods_significant_AIC %>%
  group_by(Hypothesis) %>% 
  summarise(count = n()) 

hypothesis <- factor(c(0,1,2,3,4))
number_loci <- c(3232,1254,174,314,8)

overall_summary <- data.frame(hypothesis, number_loci)

ggplot(data = overall_summary, mapping = aes(x = hypothesis, y = number_loci)) + geom_col(fill = "black") +
  theme_classic() + xlab("Hypothesis") + ylab("Number of Loci")
