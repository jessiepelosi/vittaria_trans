library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(mixtools)

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

v_app_fam_ks <- read.delim("../ks_plots/NDUV-cds-ALL.fasta.ks.tsv") #Vittaria appalachiana 
v_lin_fam_ks <- read.delim("../ks_plots/SKYV-all-cds.fasta.ks.tsv") #Vittaria lineata

# Filter out recent tandem duplicates and saturated duplicates (0 < Ks < 4)

v_app_fam_ks_filtered <- v_app_fam_ks %>% 
  select(Ks) %>% 
  filter(Ks > 0) %>% 
  filter(Ks < 4)

v_lin_fam_ks_filtered <- v_lin_fam_ks %>% 
  select(Ks) %>% 
  filter(Ks > 0) %>% 
  filter(Ks < 4)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)

v_app_fam_normalmixEM <- normalmixEM(v_app_fam_ks_filtered$Ks, mu = c(0.1,0.61,2.16)) #values of mu determined from bgmms in wgd 
summary(v_app_fam_normalmixEM)
#summary of normalmixEM object:
#          comp 1   comp 2   comp 3
#lambda 0.0851117 0.346349 0.568539
#mu     0.0515052 0.516681 2.401755
#sigma  0.0400979 0.229875 0.885753
#loglik at estimate:  -8875.475  

v_lin_fam_normalmixEM <- normalmixEM(v_lin_fam_ks_filtered$Ks, k=3) #fit 3 components 
summary(v_lin_fam_normalmixEM)
#summary of normalmixEM object:
#          comp 1   comp 2   comp 3
#lambda 0.0876850 0.374254 0.538061
#mu     0.0562144 0.579358 2.405852
#sigma  0.0423732 0.276331 0.876412
#loglik at estimate:  -7252.284 

# Plot Ks paralog distribution and mixture model in ggplot2
# mapply function from https://dozenoaks.twelvetreeslab.co.uk/2019/06/mixture-models/

ggplot(data = v_app_fam_ks_filtered, mapping = aes(x = Ks)) + geom_histogram(fill = "cornflowerblue", alpha = 0.5, binwidth = 0.025) +
  theme_classic() + ylim(0,130) +
  mapply(
    function(mean, sd, lambda,n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd)) * n * binwidth * lambda
        }
      )
    },
    mean = v_app_fam_normalmixEM[["mu"]],
    sd = v_app_fam_normalmixEM[["sigma"]],
    lambda = v_app_fam_normalmixEM[["lambda"]],
    n = length(v_app_fam_ks_filtered$Ks),
    binwidth = 0.025
  )

ggsave(file = "Vittaria_appalachiana_mix_model_familywise.png", dpi = 300, height = 5, width =7)

ggplot(data = v_lin_fam_ks_filtered, mapping = aes(x = Ks)) + geom_histogram(fill = "tomato", alpha = 0.5, binwidth = 0.025) +
  theme_classic() + ylim(0,130) +
  mapply(
    function(mean, sd, lambda,n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd)) * n * binwidth * lambda
        }
      )
    },
    mean = v_lin_fam_normalmixEM[["mu"]],
    sd = v_lin_fam_normalmixEM[["sigma"]],
    lambda = v_lin_fam_normalmixEM[["lambda"]],
    n = length(v_lin_fam_ks_filtered$Ks),
    binwidth = 0.025
  )

ggsave(file = "Vittaria_lineata_mix_model_family-wise.png", dpi =300, height = 5, width = 7)

# Pair-wise distributions from wgd 
