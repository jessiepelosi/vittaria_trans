## Gene family analysis
## Import results from Count under Wagner Parsimony, gain penalty = 1.2

library(dplyr)

gf_count <- read.delim("count_genefamilies_wagner.txt")

gf_count_losses <- gf_count %>% 
  filter(Losses > 0) %>% 
  filter(NDUV.v.appalachiana < 1) %>% 
  filter(NDUV.v.appalachiana < X3)

OGs_lost_NDUV <- select(gf_count_losses, Family)

gf_count_con <- gf_count %>% #checked against GUI results 
  filter(Contractions > 0) %>% 
  select(Family, NDUV.v.appalachiana, X3) %>% 
  filter(NDUV.v.appalachiana <= (0.5*X3)) %>%
  filter(NDUV.v.appalachiana == 1)
  
OGs_cont_NDUV <- select(gf_count_con, Family)

gf_count_gain <- gf_count %>% 
  filter(Gains > 0) %>%
  filter(X3 == 0) %>% 
  filter(NDUV.v.appalachiana > 0)

OGs_gained_NDUV <- select(gf_count_gain, Family)

gf_count_exp <- gf_count %>% #Checked against GUI results 
  filter(Expansions > 0) %>%
  filter(X3 == 1) %>% 
  filter(NDUV.v.appalachiana >= 1.5*X3)

OGs_expa_NDUV <- select(gf_count_exp, Family)


