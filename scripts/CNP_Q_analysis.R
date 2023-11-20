library(kableExtra)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggtext)
library(ggpubr)

CNP_Quotas <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                         sheet = "final_table") |>
  mutate(Culture = str_replace(Culture, "mcyB", "<i>mcyB</i>"))
matrix <- read_excel("./data/Microcystis Quota Estimates.xlsx", 
                     sheet = "corr_matrix") |>
  mutate(Culture = str_replace(Culture, 'mcyB', '*mcyB*'))

barplot_ksN <- ggplot(data = matrix, aes(x = Culture, y = ks_N, 
                                         ymin = ks_N - ks_N_SE,
                                         ymax = ks_N + ks_N_SE,
                                         fill = State)) +
  geom_bar(position = "dodge", stat = 'identity', aes(fill = State)) +
  geom_text(inherit.aes = F,
            aes(label = sig_ksN, x = Culture, group = State, y = ks_N + ks_N_SE + ks_N_SE*.1), 
            position = position_dodge(.9), size = 8, show.legend = F, check_overlap = T,
            data = filter(matrix, ! sig_ksN == "ns")) +
  coord_cartesian(ylim = c(0, 499), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(K)[s]~(μgL^-1))) +
  rremove("xlab") +
  theme_pubr(margin = F, border = T,base_family = 'arial', base_size = 8) +
  theme(axis.text.x = element_markdown(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 
barplot_mumaxN <- ggplot(data = matrix, aes(x = Culture, y = mumax_N,
                                            ymin = mumax_N - mumax_N_SE,
                                            ymax = mumax_N + mumax_N_SE,
                                            fill = State)) +
  geom_bar(position = "dodge", stat = 'identity', aes(fill = State)) +
  geom_text(inherit.aes = F,
            aes(label = sig_mumaxN, x = Culture, group = State, y = mumax_N + mumax_N_SE + mumax_N_SE*.1), 
            position = position_dodge(.9), size = 8, show.legend = F, check_overlap = T,
            data = filter(matrix, ! sig_mumaxN == "ns")) +
  coord_cartesian(ylim = c(0, 0.65), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(μ)[max]~(day^-1))) +
  theme_pubr(margin = F, border = T,base_family = 'arial', base_size = 8) +
  theme(axis.text.x = element_markdown(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 

barplot_QN <- ggplot(data = matrix, aes(x = Culture, y = Q_N)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  coord_cartesian(ylim = c(0,165), expand = F) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(Q)[N]~(fmol~cell^-1))) +
  rremove("xlab") +
  theme_pubr(margin = F, border = T, base_family = 'arial',base_size = 8) +
  theme(axis.text.x = element_markdown(hjust = 1, angle = 70),
        axis.title.x = element_blank(),
        legend.position = "right"
  ) 
N_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = c('[A]','[B]','[C]'), 
                     barplot_ksN, barplot_mumaxN, barplot_QN,font.label = list(size = 9.5, color = "black"))

ggsave('./figures/figure_pieces/figure5ABC.png',N_plots,width = 160, height = 152, units = "mm", dpi = 600 )
ggsave('./figures/figure_pieces/figure5ABC.svg',N_plots,width = 160, height = 152, units = "mm", dpi = 600 )

barplot_ksP <- ggplot(data = matrix, aes(x = Culture, y = ks_P, 
                                         ymin = ks_P - ks_P_SE,
                                         ymax = ks_P + ks_P_SE,
                                         fill = State)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge()) +
  geom_text(inherit.aes = F,
            aes(label = sig_ksP, x = Culture, group = State, y = ks_P + ks_P_SE + ks_P_SE*.1), 
            position = position_dodge(.9), size = 8, show.legend = F, check_overlap = T,
            data = filter(matrix, ! sig_ksP == "ns")) +
  coord_cartesian(ylim = c(0,9), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(K)[s]~(μgL^-1))) +
  theme_pubr(margin = F, border = T,base_family = 'arial', base_size = 8) +
  theme(axis.text.x = element_markdown(hjust = 1, angle = 70),
        legend.position = "right"
  ) 
barplot_mumaxP <- ggplot(data = matrix, aes(x = Culture, y = mumax_P,
                                            ymin = mumax_P - mumax_P_SE,
                                            ymax = mumax_P + mumax_P_SE,
                                            fill = State)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge(.9)) +
  geom_text(inherit.aes = F,
            aes(label = sig_mumaxP, x = Culture, group = State, y = mumax_P + mumax_P_SE + mumax_P_SE*.1), 
            position = position_dodge(.9), size = 8, show.legend = F, check_overlap = T,
            data = filter(matrix, ! sig_mumaxP == "ns")) +
  coord_cartesian(ylim = c(0,0.65), expand = F) +
  geom_errorbar(position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(μ)[max]~(day^-1))) +
  theme_pubr(margin = F, border = T, base_size = 8) +
  theme(axis.text.x = element_markdown(hjust = 1, angle = 70),
        legend.position = "right"
  )  

barplot_QP <- ggplot(data = matrix, aes(x = Culture, y = Q_P)) +
  geom_bar(stat = 'identity', aes(fill = State), position = position_dodge(.9)) +
  coord_cartesian(ylim = c(0,3.5), expand = F) +
  scale_fill_manual(values = c("Axenic" = "darkseagreen1",
                               "Xenic" = "forestgreen")) +
  scale_x_discrete(limits = c("PCC 7806 Δ*mcyB*","PCC 7806","PCC 9701", "NIES-843")) +
  ylab(expression(italic(Q)[P]~(fmol~cell^-1))) +
  theme_pubr(margin = F, border = T,base_family = 'arial', base_size = 8) +
  theme(axis.text.x = element_markdown(angle = 70, hjust = 1),
        legend.position = "right"
  ) 

P_plots <- ggarrange(common.legend = T, legend = "bottom",nrow = 1, labels = c('[D]','[E]','[F]'), 
                     barplot_ksP, barplot_mumaxP, barplot_QP, font.label = list(size = 9.5, color = "black"))

ggsave('./figures/figure_pieces/figure5DEF.png',P_plots,width = 160, height = 152, units = "mm", dpi = 600)
ggsave('./figures/figure_pieces/figure5DEF.svg',P_plots,width = 160, height = 152, units = "mm", dpi = 600)

all_plots <- ggarrange(common.legend = T, legend = "bottom",
                       nrow = 2,ncol = 3, hjust = 0, labels = c('[A]','[B]','[C]',
                                                     '[D]','[E]','[F]'),
                       font.label = list(size = 9.5, color = "black"),
                       barplot_ksN, barplot_mumaxN, barplot_QN,
                       barplot_ksP, barplot_mumaxP, barplot_QP)
ggsave('./figures/Finalized Figs for Publication/figure5_wcolor.png',all_plots,width = 160, height = 250, units = "mm", dpi = 600)
ggsave('./figures/Finalized Figs for Publication/figure5_wcolor.svg',device = svglite::svglite,
       all_plots,width = 160, height = 250, units = "mm", dpi = 600)

Q_lm <- lm(Q_N~Q_P, matrix)
mu_lm <- lm(mumax_N~ mumax_P, matrix)
ks_lm <- lm(ks_N~ks_P, matrix)

#Tradeoff plots
p1 <- ggplot(matrix, aes(x = Q_N, y = Q_P)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  labs(y = "*Q*<sub>P</sub> fmol • cell<sup>-1</sup>",
       x = "*Q*<sub>N</sub> fmol • cell<sup>-1</sup>") +
  #facet_wrap(~state) +
  stat_cor(method = "pearson", size = 8/.pt,
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_pubr(base_size = 8, base_family = 'arial', border = T) +
  theme(axis.title.x = ggtext::element_textbox(),
        axis.title.y = ggtext::element_textbox(orientation = 'left-rotated'))

p2 <- ggplot(matrix, aes(x = mumax_N, y = mumax_P)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  labs(y = "*μ*<sub>max</sub> day<sup>-1</sup> P",
       x = "*μ*<sub>max</sub> day<sup>-1</sup> N") +
  #facet_wrap(~state) +
  stat_cor(method = "pearson", size = 8/.pt,
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_pubr(base_size = 8, base_family = 'arial', border = T) +
  theme(axis.title.x = ggtext::element_textbox(),
        axis.title.y = ggtext::element_textbox(orientation = 'left-rotated'))

p3 <- ggplot(matrix, aes(x = ks_N, y = ks_P)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y~x, se = F, color = "black") +
  #facet_wrap(~state) +
  labs(y = "*K*<sub>s</sub> μg/L P",
       x = "*K*<sub>s</sub> μg/L N") +
  stat_cor(method = "pearson", size = 8/.pt,
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) + 
  theme_pubr(base_size = 8, base_family = 'arial', border = T) +
  theme(axis.title.x = ggtext::element_textbox(),
        axis.title.y = ggtext::element_textbox(orientation = 'left-rotated'))

tradeoff_plots <- ggarrange(hjust = 0, nrow = 1, labels = c('[A]','[B]','[C]'), 
          p3, p2, p1, font.label = list(size = 9.5, color = "black"))

ggsave("./figures/Finalized Figs for Publication/FigureS3_final.png",tradeoff_plots,width = 160, height = 100, units = "mm", dpi = 600)

kable(CNP_Quotas,escape = F, digits = 4, align = c("l") ,col.names = c("Culture", "<i>μ</i><sub>max</sub> day<sup>-1</sup> (SE)", 
                                                                                                            "<i>K</i><sub>s</sub> (μg/L) (SE)",
                                                                                                            "Initial Slope day<sup>-1</sup> • μM<sup>-1</sup>",
                                                                                                            "<i>Q</i><sub>P</sub> (SD) fmol • cell<sup>-1</sup>",
                                                                                                            "P<sup>*</sup>, d = 0.2 day<sup>-1</sup>",
                                                                                                            "<i>μ</i><sub>max</sub> day<sup>-1</sup> (SE)","<i>K</i><sub>s</sub> (μg/L) (SE)", 
                                                                                                            "Initial Slope day<sup>-1</sup> • μM<sup>-1</sup>",
                                                                                                            "<i>Q</i><sub>N</sub> (SD) fmol • cell<sup>-1</sup>", 
                                                                                                            "N<sup>*</sup>, d = 0.2 day<sup>-1</sup>",
                                                                                                            "<i>Q</i><sub>C</sub> (SD) fmol • cell<sup>-1</sup>",
                                                                                                            "N:P","C:N","C:P" ))  %>%
  kable_styling("striped", full_width = T) %>%
  column_spec(c(4,9),width_min = "1in") %>%
  add_header_above(c("", "Phosphorus" = 5, "Nitrogen" = 5, "Carbon" = 1, "Elemental Ratios" = 3)) %>% save_kable(file = "./tables/table1_final.doc", zoom = 1.5) 

  

