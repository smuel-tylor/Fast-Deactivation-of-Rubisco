####Comparing Method Type
#Suppl Fig - bicarbonate tests
#Modfied from scripts provided by Emmanuel Gonzalez-Escobar by Samuel Taylor
# to fit with GitHub submission

library(here)
library(ggplot2)
library(patchwork)
library(emmeans)

#Data for Supplementary Fig 3
SDF3 <- read.csv(here("data/Data_FigS3.csv"))

head(SDF3)
SDF3$Condition

#Statistical tests
#The interest here, is in establishing whether incubation time
# affects activation state
lmSDF3_AS <- lm(Activation.State ~ Condition, data=SDF3)
anova(lmSDF3_AS)
emmeans(lmSDF3_AS, pairwise ~ Condition)
#this is sig diff between the 20 min and 1 h incubation

lmSDF3_Vi <- lm(Vi_umol_m2_s ~ Condition, data=SDF3)
emmeans(lmSDF3_Vi, pairwise ~ Condition)
#ns

lmSDF3_Vt <- lm(Vt3_umol_m2_s ~ Condition, data=SDF3)
emmeans(lmSDF3_Vt, pairwise ~ Condition)
#ns


level_order_3 = c("20min", "40min", "60min")
SDF3$Condition <- factor(SDF3$Condition, levels = level_order_3)

#Light System Incubation Times Compared----
#Activation state
FigS3_ActS_Incubation <- ggplot(SDF3, aes(x=Condition, y=Activation.State)) + 
  geom_boxplot(fill = c("#A6CEE3"),
               color = "black",
               width = 0.50) +
  geom_point() +
  xlab("Incubation Time") +
  ylab("Activation State (%)") +
  ylim(0,100) +
  annotate("text",
           x = c(1:3), y = 97, label = c("a", "ab", "b"),
           fontface = "bold", size = 7
           ) +
  theme_bw() +
  theme(
    plot.title = element_text(color = "black", size  =20, face = "bold.italic"),
    axis.title = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "white", size = 20),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = 'black', size = 20),
    axis.ticks.length = unit(-0.25, "cm"),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

FigS3_ActS_Incubation

#Vi
FigS3_Vi_Incubation = ggplot(SDF3, aes(x=Condition, y=Vi_umol_m2_s)) + 
  geom_boxplot(fill = c("#A6CEE3"),
               color = "black", width = 0.50) +
  geom_point() +
  xlab("Incubation Time (min)") +
  ylab(expression("Initial Activity "*(mu*mol~~m^-2~~s^-1))) +
  ylim(0, 80) +
  annotate("text",
           x = c(1:3), y = 77, label = "a",
           fontface = "bold", size = 7
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold.italic"),
    axis.title = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = 'white', size = 20),
    axis.text.y = element_text(color = 'black', size = 20),
    axis.title.x = element_text(color = 'white', size = 20),
    axis.title.y = element_text(color = 'black', size = 20),
    axis.ticks.length = unit(-0.25, "cm"),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

FigS3_Vi_Incubation

#Vt
FigS3_Vt_Incubation = ggplot(SDF3, aes(x = Condition, y = Vt3_umol_m2_s)) + 
  geom_boxplot(fill = c("#A6CEE3"),
               color = "black",
               width = 0.50) +
  geom_point() +
  ylab(expression("Total Activity "*(mu*mol~~m^-2~~s^-1))) +
  xlab("Incubation Time") +
  ylim(0, 80) +
  annotate("text",
           x = c(1:3), y = 77, label = "a",
           fontface = "bold", size = 7
  ) +
  theme_bw() +
  scale_x_discrete(labels=c("20min" = "20 min\n(N = 6)",
                            "40min" = "40 min\n(N = 5)",
                            "60min" = "60 min\n(N = 6)")
                   ) + 
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold.italic"),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = 'black', size = 20),
    axis.title.x = element_text(color = "black", size = 20, vjust = -1),
    axis.title.y = element_text(color = 'black', size = 20),
    axis.ticks.length = unit(-0.25, "cm"),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

FigS3_Vt_Incubation

#compile multipanel figure
Layout_RubiscoPaper_FigS3 = FigS3_ActS_Incubation /
                              FigS3_Vi_Incubation /
                                FigS3_Vt_Incubation +
  plot_layout(widths = c(1, 8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20, face =  "bold")
        )

Layout_RubiscoPaper_FigS3

ggsave(here("output/FigureS3.pdf"),
       Layout_RubiscoPaper_FigS3,
       width=8, height=15
)
