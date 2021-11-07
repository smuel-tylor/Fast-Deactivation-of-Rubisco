####Comparing Method Type
#Suppl Fig 2 only Frozen and GH Lights ----
#modfied from scripts provided by Emmanuel Gonzalez-Escobar by Samuel Taylor
# to fit with Github submission

library('here')
library('ggplot2')
library('patchwork')
library('readxl')
library('lme4')
library('emmeans')
library('dplyr')
library('grid')


#BT4 DF----
# names indicates data from experiment Bicarbonate Trials 4
BT4 <- read.csv(here("data/Data_FigS4.csv"))

head(BT4)

#Statistical tests
#The interest here, is in establishing whether incubation time and NaHCO3
# affect activation state within each genotype
#There is repeated measures because "Replicate" represents
# technical repeats using the same material
lmerBT4 <- lmer(ActivationState ~ Bicarbonate_mM * Duration_min * Genotype +
                  (1|Replicate),
                data=BT4
)
emmeans(lmerBT4, pairwise ~ Bicarbonate_mM*Duration_min|Genotype)

#plots

#full colours
blck = rgb(0,0,0, maxColorValue=255)
dgry = rgb(0.5*255,0.5*255,0.5*255, maxColorValue=255)
dgrn = rgb(70,190,110, maxColorValue=255)
lgrn = rgb(135,240,115, maxColorValue=255)

#plot for IT86D-1010
A <- BT4 %>%
  filter(Genotype == "IT86D-1010") %>%
  ggplot(aes(x=interaction(Bicarbonate_mM, Duration_min), y=ActivationState)) + 
  geom_boxplot(fill=rep(dgrn, 4),
               colour="black"
  ) +
  geom_point() +
  xlab("") +
  ylab("Activation State (%)") +
  ylim(0,100) +
  theme_bw() +
  #annotate("text",
  #         x=c(1:4),
  #         y=-10,
  #         label=c("MES","NaHCO[3]","MES","NaHCO[3]"),
  #         parse=TRUE
  #) +
  #annotate("text", x=c(1.5, 3.5), y=-15, label=c("40 min", "80 min")) +
  annotate("text", x=3.5, y=30, label="IT86D-1010", fontface="bold", size=6) +
  annotate("text", x=c(1:4), y=97, label="a", fontface="bold", size=5) +
  coord_cartesian(ylim=c(0, 100), clip="off") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    text = element_text(colour = 'black', size=20),
    axis.text.x = element_blank(),#element_text(colour = 'white', size=14),
    axis.title.x = element_blank(),#element_text(vjust=-1),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.title.y = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA)
  )

A

#plot for IT82E-16
B <- BT4 %>%
  filter(Genotype == "IT82E-16") %>%
  ggplot(aes(x=interaction(Bicarbonate_mM, Duration_min), y=ActivationState)) + 
  geom_boxplot(fill=rep(lgrn, 4),
               colour="black"
  ) +
  geom_point() +
  xlab("") +
  ylab("Activation State (%)") +
  ylim(0,100) +
  theme_bw() +
  #annotate("text",
  #         x=c(1:4),
  #         y=-10,
  #         label=c("MES","NaHCO[3]","MES","NaHCO[3]"),
  #         parse=TRUE
  #) +
  #annotate("text", x=c(1.5, 3.5), y=-15, label=c("40 min", "80 min")) +
  annotate("text", x=3.5, y=30, label="IT82E-16", fontface="bold", size=6) +
  annotate("text", x=c(1:4), y=97, label="a", fontface="bold", size=5) +
  coord_cartesian(ylim=c(0, 100), clip="off") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    text = element_text(colour = 'black', size=20),
    axis.text.x = element_blank(),#element_text(colour = 'white', size=14),
    axis.title.x = element_blank(),#element_text(vjust=-1),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.title.y = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA)
  )

B

#plot for TVNu-1948
C <- BT4 %>%
  filter(Genotype == "TVNu-1948") %>%
  ggplot(aes(x=interaction(Bicarbonate_mM, Duration_min), y=ActivationState)) + 
  geom_boxplot(fill=rep(blck, 4),
               colour="black"
  ) +
  geom_point(colour = "gray") +
  xlab("") +
  ylab("Activation State (%)") +
  ylim(0,100) +
  theme_bw() +
  #annotate("text",
  #         x=c(1:4),
  #         y=-10,
  #         label=c("MES","NaHCO[3]","MES","NaHCO[3]"),
  #         parse=TRUE
  #) +
  #annotate("text", x=c(1.5, 3.5), y=-15, label=c("40 min", "80 min")) +
  annotate("text", x=3.5, y=30, label="V. sp. Savi", fontface="bold.italic", size=6) +
  annotate("text", x=c(1:4), y=97, label="a", fontface="bold", size=5) +
  coord_cartesian(ylim=c(0, 100), clip="off") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    text = element_text(colour = 'black', size=20),
    axis.text.x = element_blank(),#element_text(colour = 'white', size=14),
    axis.title.x = element_blank(),#element_text(vjust=-1),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.title.y = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA)
  )

C

#plot for V.adenantha
D <- BT4 %>%
  filter(Genotype == "Adenantha") %>%
  ggplot(aes(x=interaction(Bicarbonate_mM, Duration_min), y=ActivationState)) + 
  geom_boxplot(fill=rep(dgry, 4),
               colour="black"
  ) +
  geom_point() +
  xlab("") +
  ylab("Activation State (%)") +
  theme_bw() +
  annotate("text",
           x=c(1:4),
           y=-12,
           label=c("MES","NaHCO[3]","MES","NaHCO[3]"),
           parse=TRUE,
           size = 5.5
  ) +
  annotate("text", x=c(1.5, 3.5), y=-22, label=c("40 min", "80 min"), size = 5.5) +
  annotate("text", x=3.5, y=30, label="V. adenantha", fontface="bold.italic", size=6) +
  annotate("text", x=c(1:4), y=97, label=c("b", rep("a",3)), fontface="bold", size=5) +
  coord_cartesian(ylim=c(0, 100), clip="off") +
  theme(
    plot.margin = unit(c(0.2, 1, 2.5, 1), "lines"),
    text = element_text(colour = 'black', size=20),
    axis.text.x = element_blank(),#element_text(colour = 'white', size=14),
    axis.title.x = element_blank(),#element_text(vjust=-1),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.title.y = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA)
    )

D

Layout_RubiscoPaper_FigureS4 = wrap_elements(
  textGrob("Activation State (%)",
                 rot=90,
           gp=gpar(fontsize=20) 
                 )) +
  A / B / C / D +
  plot_layout(widths = c(1, 8)) +
  plot_annotation(tag_levels = list(c("", "a", "b", "c", "d")),
                  caption = "N = 4",
                  ) &
  theme(plot.tag = element_text(size = 20, face =  "bold"),
        plot.caption = element_text(size = 20)
        )

Layout_RubiscoPaper_FigureS4

ggsave(here("output/FigureS4.pdf"),
       Layout_RubiscoPaper_FigureS4,
       width=8, height=15
) ##Wait to confirm measurements


