#################################
# TE_repeats.R
# Jessie Pelosi 
# Last Modified October 22, 2021 
#################################

### REPEATS 
library(dplyr)
library(ggplot2)
library(chisq.posthoc.test)


repeats <- read.csv("all_repeats.csv")

## Plot repeats

ggplot(data = repeats, mapping = aes(x = Species, y = X._seq, fill = Type)) + geom_col() +
  ylab("% of Transcriptome") + theme_classic() + scale_fill_brewer(palette = "Set3") 

ggsave("Repeats_all_transcriptome.pdf", height = 7, width = 5)

ggplot(data = repeats, mapping = aes(x = Species, y = kb, fill = Type)) + geom_col() +
  ylab("Length (Kb)") + theme_classic() + scale_fill_brewer(palette = "Set3") 

ggsave("Repeats_all_length_transcriptome.pdf", height = 7, width = 5)

retroelements <- repeats %>% 
  filter(Type == "Retroelement")

ggplot(data = retroelements, mapping = aes(x = Species, y = kb, fill = Element)) + geom_col() +
  ylab("Length (Kb)") + theme_classic() + scale_fill_brewer(palette = "Set2") 

ggsave("retroelements_length.pdf", height = 7, width = 5)

ggplot(data = retroelements, mapping = aes(x = Species, y = X._seq, fill = Element)) + geom_col() +
  ylab("% of Transcriptome") + theme_classic() + scale_fill_brewer(palette = "Set2") 

ggsave("retroelements_percent.pdf", height = 7, width = 5)

LINEs <- repeats %>% 
  filter(Element == "LINE")

ggplot(data = LINEs, mapping = aes(x = Species, y = kb, fill = Family)) + geom_col() +
  ylab("Length (Kb)") + theme_classic() + scale_fill_brewer(palette = "Set3")

ggsave("LINEs_length.png", dpi = 300, height = 7, width = 5)

ggplot(data = LINEs, mapping = aes(x = Species, y = X._seq, fill = Family)) + geom_col() +
  ylab("% of Transcriptome") + theme_classic() + scale_fill_brewer(palette = "Set3")

ggsave("LINEs_percent.png", dpi = 300, height = 7, width =5)

LTRs <- repeats %>% 
  filter(Element == "LTR")

ggplot(data = LTRs, mapping = aes(x = Species, y = kb, fill = Family)) + geom_col() +
  ylab("Length (Kb)") + theme_classic() + scale_fill_brewer(palette = "Set2")

ggsave("LTR_length.png", dpi = 300, height = 7, width =5)

ggplot(data = LTRs, mapping = aes(x = Species, y = X._seq, fill = Family)) + geom_col() +
  ylab("% of Transcriptome") + theme_classic() + scale_fill_brewer(palette = "Set2")

ggsave("LTR_percent.png", dpi = 300, height = 7, width =5)

DNA_elements <- repeats %>% 
  filter(Type == "DNA_transposon")

ggplot(data = DNA_elements, mapping = aes(x = Species, y = kb, fill = Family)) + geom_col() +
  ylab("Length(Kb)") + theme_classic() + scale_fill_brewer(palette = "Set3")

ggsave("DNA_elements_length.pdf", height = 7, width = 5)

ggplot(data = DNA_elements, mapping = aes(x = Species, y = X._seq, fill = Family)) + geom_col() +
  ylab("% of Transcriptome") + theme_classic() + scale_fill_brewer(palette = "Set3")

ggsave("DNA_elements_percent.pdf", height = 7, width = 5)

#### Contigency Table Tests ####

Number <- repeats %>% 
  dplyr::select(Species, Family, Number_of_elements)

Table = xtabs(Number_of_elements ~ Species + Family,
              data=Number)

chisq.test(Table)

chisq.posthoc.test(Table, method = "bonferroni", round = 7)


Prop <- repeats %>% 
  dplyr::select(Species, Family, X._seq)

Table = xtabs(X._seq ~ Species + Family,
              data=Prop)

chisq.test(Table)

chisq.posthoc.test(Table, method = "bonferroni", round = 7)


Length <- repeats %>% 
  dplyr::select(Species, Family, Length)

Table = xtabs(Length ~ Species + Family,
              data=Length)

chisq.test(Table)

chisq.posthoc.test(Table, method = "bonferroni", round = 7)
