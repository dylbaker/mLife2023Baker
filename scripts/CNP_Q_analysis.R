library(kableExtra)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggtext)
library(ggpubr)

CNP_Quotas <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                         sheet = "final_table") 
matrix <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                     sheet = "corr_matrix")

barplot_ksN <- ggplot(data = matrix, aes(x = Culture, y = ks_N, 
                                         ymin = ks_N - ks_N_SE,
                                         ymax = ks_N + ks_N_SE,
                                         fill = State)) +
  geom_bar(position = "dodge", stat = 'identity', aes(fill = State)) +
  coord_cartesian(ylim = c(0, 499), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(K[s]~(μgL^-1))) +
  rremove("xlab") +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 
barplot_mumaxN <- ggplot(data = matrix, aes(x = Culture, y = mumax_N,
                                            ymin = mumax_N - mumax_N_SE,
                                            ymax = mumax_N + mumax_N_SE,
                                            fill = State)) +
  geom_bar(position = "dodge", stat = 'identity', aes(fill = State)) +
  coord_cartesian(ylim = c(0, 0.65), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(μ[max]~(day^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 

barplot_QN <- ggplot(data = matrix, aes(x = Culture, y = Q_N)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,165), expand = F) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(Quota~(fmol~cell^-1))) +
  rremove("xlab") +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 
N_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = "AUTO", 
                     barplot_ksN, barplot_mumaxN, barplot_QN)

ggsave('./figures/figure_pieces/figure5ABC.png',N_plots,width = 180, height = 152, units = "mm", dpi = 600 )
ggsave('./figures/figure_pieces/figure5ABC.svg',N_plots,width = 180, height = 152, units = "mm", dpi = 600 )

barplot_ksP <- ggplot(data = matrix, aes(x = Culture, y = ks_P, 
                                         ymin = ks_P - ks_P_SE,
                                         ymax = ks_P + ks_P_SE,
                                         fill = State)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,9), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(K[s]~(μgL^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        legend.position = "right"
  ) 
barplot_mumaxP <- ggplot(data = matrix, aes(x = Culture, y = mumax_P,
                                            ymin = mumax_P - mumax_P_SE,
                                            ymax = mumax_P + mumax_P_SE,
                                            fill = State)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,0.65), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(μ[max]~(day^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        legend.position = "right"
  )  

barplot_QP <- ggplot(data = matrix, aes(x = Culture, y = Q_P)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,3.5), expand = F) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 ΔmcyB","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(Quota~(fmol~cell^-1))) +
  theme_pubr(margin = F, border = T, base_size = 16) +
  theme(axis.text.x = element_text(hjust = 1, angle = 70),
        legend.position = "right"
  ) 

P_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = "AUTO", 
                     barplot_ksP, barplot_mumaxP, barplot_QP)

ggsave('./figures/figure_pieces/figure5DEF.png',P_plots,width = 180, height = 152, units = "mm", dpi = 600)
ggsave('./figures/figure_pieces/figure5DEF.svg',P_plots,width = 180, height = 152, units = "mm", dpi = 600)

all_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 2,ncol = 3, labels = "AUTO", 
                       barplot_ksN, barplot_mumaxN, barplot_QN,
                       barplot_ksP, barplot_mumaxP, barplot_QP)
ggsave('./figures/Finalized Figs for Publication/figure5_wcolor.png',all_plots,width = 180, height = 300, units = "mm", dpi = 600)
ggsave('./figures/Finalized Figs for Publication/figure5_wcolor.svg',all_plots,width = 180, height = 300, units = "mm", dpi = 600)

Q_lm <- lm(Q_N~Q_P, matrix)
mu_lm <- lm(mumax_N~ mumax_P, matrix)
ks_lm <- lm(ks_N~ks_P, matrix)

p1 <- ggplot(matrix, aes(x = Q_N, y = Q_P)) +
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

png(filename = "./figures/Finalized Figs for Publication/FigureS3_final.png",
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
  add_header_above(c("", "Phosphorus" = 5, "Nitrogen" = 5, "Carbon" = 1, "Elemental Ratios" = 3)) %>% save_kable(file = "./tables/new_table1.doc", zoom = 1.5)
  

