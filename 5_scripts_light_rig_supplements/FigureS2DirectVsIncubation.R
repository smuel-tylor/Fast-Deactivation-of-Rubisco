####Comparing Method Type
#Suppl Fig 2 only Frozen and GH Lights ----
library(here)
library(ggplot2)

#This data from Leaf Disc Trial 2
LDT2 <- read.csv("data/Data_FigS2.csv")

#level_order_2 = c("Frozen", "Greenhouse Lights", "60min")
LDT2$Sample <- factor(LDT2$Sample)

lmLDT2_AS <- lm(Act.S ~ Sample, data=LDT2)
anova(lmLDT2_AS)
emmeans(lmLDT2_AS, pairwise ~ Sample)
#this is sig diff between the 20 min and 1 h incubation

lmLDT2_Vi <- lm(Vi ~ Sample, data=LDT2)
emmeans(lmLDT2_Vi, pairwise ~ Sample)
#ns

lmLDT2_Vt <- lm(Vt ~ Sample, data=LDT2)
emmeans(lmLDT2_Vt, pairwise ~ Sample)


fcols <- c("lightgreen", "paleturquoise2", "lavender")

##Title: Methods Compared
#SupplFig2A_ActS Methods Compared----
#Activation state figure 2A for Supplementary data
FigS2_ActS_buffer = ggplot(LDT2, aes(x=Sample, y=Act.S)) + 
  geom_boxplot(fill = fcols, color="black", width = 0.50) +
  geom_point() +
  xlab("Method") +
  ylab("Activation State (%)") +
  ylim(0,100) +
  scale_x_discrete(labels=c("FROZEN" = "Directly Frozen",
                            "MES" = "MES",
                            "WATER" = expression(H[2]*O)
                            )
                   ) + 
  annotate("text",
           x = c(1:3), y = 100, label="a",
           fontface = "bold", size = 7
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = 'black', size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

FigS2_ActS_buffer

#SupplFig2A2_Vi Methods Compared----
#Activation state figure 2A for Supplementary data
FigS2_Vi_buffer = ggplot(LDT2, aes(x=Sample, y=Vi)) + 
  geom_boxplot(fill = fcols, color = "black", width = 0.50) +
  geom_point() +
  xlab("Method") +
  ylab("Initial Activity"~~(mu*mol~~m^-2~~s^-1)) +
  ylim(0, 80) +
  theme_bw() +
  scale_x_discrete(labels = c("FROZEN" = "Directly Frozen",
                              "MES" = "MES",
                              "WATER" = expression(H[2]*O)
  )
  ) +
  annotate("text",
           x = c(1:3), y = 80, label="a",
           fontface = "bold", size = 7
  ) +
  theme(
    axis.title.y = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = 'black', size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

FigS2_Vi_buffer


#SupplFig2A3_Vt Methods Compared----
#Activation state figure 2A for Supplementary data
FigS2_Vt_buffer = ggplot(LDT2, aes(x=Sample, y=Vt)) + 
  geom_boxplot(fill = fcols, color = "black", width = 0.50) +
  geom_point() +
  xlab("60 min greenhouse incubation") +
  ylab("Total Activity"~~(mu*mol~~m^-2~~s^-1)) +
  ylim(0, 80) +
  theme_bw() +
  scale_x_discrete(labels = c("FROZEN" = "Directly Frozen",
                              "MES" = "MES",
                              "WATER" = expression(H[2]*O)
  )
  ) +
  annotate("text",
           x = c(1:3), y = 80, label="a",
           fontface = "bold", size = 7
  ) +
  theme(
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = 'black', size = 20),
    #This one is at the bottom of the layout, so it keeps its x-axis labels
    axis.title.x = element_text(vjust = -1, hjust = 0.85), 
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

FigS2_Vt_buffer

Layout_RubiscoPaper_FigureS2 <- FigS2_ActS_buffer /
                                  FigS2_Vi_buffer /
                                    FigS2_Vt_buffer +
  plot_annotation(tag_levels = "a",
                  caption = "N = 4",
  ) &
  theme(plot.tag = element_text(size = 20, face =  "bold"),
        plot.caption = element_text(size = 20)
  )

Layout_RubiscoPaper_FigureS2

ggsave(here("output/FigureS2.pdf"),
       Layout_RubiscoPaper_FigureS2,
       width = 8, height = 15
)

