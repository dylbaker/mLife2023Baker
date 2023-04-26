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
  pivot_longer(cols = c(mumax_P, mumax_N, ks_P, ks_N, Q_N, Q_P), names_to = c('Parameter','Element'), names_sep = '_' ) |>
  mutate(Element = case_when( Element == 'N' ~ 'Nitrogen',
                              T ~ 'Phosphorus'))
matrix$facets = factor(matrix$Parameter, 
  labels = c(
  "K[s]", 
  "μ[max]", 
  "Quota"), 
  levels = c(
    "ks", 
    "mumax", 
    "Q"))

barplot_ksN <- ggplot(data = filter(matrix, Parameter == 'ks', Element == "Nitrogen"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0, 400), expand = F) +
  facet_wrap(~ Element + facets, nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(K[s]~(μgL^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  ) 
barplot_mumaxN <- ggplot(data = filter(matrix, Parameter == 'mumax', Element == "Nitrogen"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0, 0.65), expand = F) +
  facet_wrap(~ Element + facets,nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(μ[max]~(day^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  ) 
barplot_QN <- ggplot(data = filter(matrix, Parameter == 'Q', Element == "Nitrogen"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,160), expand = F) +
  facet_wrap(~ Element + facets,nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(Quota~(fmol~cell^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  ) 
N_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = "AUTO", 
                     barplot_ksN, barplot_mumaxN, barplot_QN)
ggsave('./figures/figure3.png',N_plots,width = 180, height = 152, units = "mm", dpi = 600 )

barplot_ksP <- ggplot(data = filter(matrix, Parameter == 'ks', Element == "Phosphorus"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,8), expand = F) +
  facet_wrap(~ Element + facets, nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(K[s]~(μgL^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  ) 
barplot_mumaxP <- ggplot(data = filter(matrix, Parameter == 'mumax', Element == "Phosphorus"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,0.65), expand = F) +
  facet_wrap(~ Element + facets,nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(μ[max]~(day^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  )  

barplot_QP <- ggplot(data = filter(matrix, Parameter == 'Q', Element == "Phosphorus"), aes(x = Culture, y = value)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,4), expand = F) +
  facet_wrap(~ Element + facets,nrow = 1, scales = 'free_y', labeller = label_parsed) +
  scale_fill_manual(values = c("Axenic" = "gray20",
                               "Xenic" = "gray85")) +
  ylab(expression(Quota~(fmol~cell^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        strip.text = element_text(size = 16),
        legend.position = "right"
  ) 

P_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = "AUTO", 
                     barplot_ksP, barplot_mumaxP, barplot_QP)
ggsave('./figures/figure4.png',P_plots,width = 180, height = 152, units = "mm", dpi = 600)


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
  add_header_above(c("", "Phosphorus" = 5, "Nitrogen" = 5, "Carbon" = 1, "Elemental Ratios" = 3)) %>% save_kable(file = "./tables/table2.doc", zoom = 1.5)
  

