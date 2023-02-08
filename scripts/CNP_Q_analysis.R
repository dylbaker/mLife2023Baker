library(kableExtra)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggtext)
library(ggpubr)

CNP_Quotas <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                         sheet = "final_table") 
matrix <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                     sheet = "corr_matrix") |>
  mutate(state = case_when(str_detect(Culture, "Axenic") ~ "Axenic",
                           T ~ "Xenic"))
Q_lm <- lm(Qn~Qp, matrix)
mu_lm <- lm(mumax_N~ mumax_P, matrix)
ks_lm <- lm(ks_N~ks_P, matrix)

p1 <- ggplot(matrix, aes(x = Qn, y = Qp)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  labs(y = "Q<sub>P</sub> fmol • cell<sup>-1</sup>",
       x = "Q<sub>N</sub> fmol • cell<sup>-1</sup>") +
  #facet_wrap(~state) +
  stat_cor(method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_pubr() +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12))

p2 <- ggplot(matrix, aes(x = mumax_N, y = mumax_P)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  labs(y = "μmax day<sup>-1</sup> P",
       x = "μmax day<sup>-1</sup> N") +
  #facet_wrap(~state) +
  stat_cor(method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_pubr() +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12))

p3 <- ggplot(matrix, aes(x = ks_N, y = ks_P)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  #facet_wrap(~state) +
  labs(y = "K<sub>s</sub> μg/L P",
       x = "K<sub>s</sub> μg/L N") +
  stat_cor(method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) + 
  theme_pubr() +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12))

png(filename = "./figures/supplemental_figures/regressions_figS4.png",
    res = 300,
    type = "cairo",
    units = "in",
    width = 10,
    height = 3)

(p3 | p2 | p1) + plot_annotation(tag_levels = "A")

dev.off()

kable(CNP_Quotas,escape = F, digits = 4, align = c("l") ,col.names = c("Culture", "μmax day<sup>-1</sup> (SE)", 
                                                                                                            "Ks (μg/L) (SE)",
                                                                                                            "Initial Slope day<sup>-1</sup> • μM<sup>-1</sup>",
                                                                                                            "Q<sub>P</sub> (SD) fmol • cell<sup>-1</sup>",
                                                                                                            "P<sup>*</sup>, d = 0.2 day<sup>-1</sup>",
                                                                                                            "μmax day<sup>-1</sup> (SE)","Ks (μg/L) (SE)", 
                                                                                                            "Initial Slope day<sup>-1</sup> • μM<sup>-1</sup>",
                                                                                                            "Q<sub>N</sub> (SD) fmol • cell<sup>-1</sup>", 
                                                                                                            "N<sup>*</sup>, d = 0.2 day<sup>-1</sup>",
                                                                                                            "Q<sub>C</sub> (SD) fmol • cell<sup>-1</sup>",
                                                                                                            "N:P","C:N","C:P" ))  %>%
  kable_styling("striped", full_width = T) %>%
  column_spec(c(4,9),width_min = "1in") %>%
  add_header_above(c("", "Phosphorus" = 5, "Nitrogen" = 5, "Carbon" = 1, "Elemental Ratios" = 3)) %>% save_kable(file = "./tables/table2.png", zoom = 1.5)
  

