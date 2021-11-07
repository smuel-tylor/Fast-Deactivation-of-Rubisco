#2020-01-15 (Original File)
#Graphing individual light response curves
#Dragons Den - Light Response Workflow
#Note that DF changed to newer version including new 1948 Data

#Note: DragonsDen_RAN_Original and DragonsDen_RAN_4 are different DFs
#Note 2: The means DF is based on Original DF values

#Update 2020-07-16 (RAN_6, built on top final version RAN_5)
#Formatting Light Response figures to match ST colour scheme0

#Update 2021-08-16 Code for Publication figures from Original 20200716 DragonsDen_RAN_6 file extracted.
#Sheets that are used here:  DragonsDen_D_Processed, DragonsDen_E_Processed, DragonsDen_N_Processed, DragonsDen_T_Processed

#Update 2021-08-18 Samuel Taylor
# modified to fit with GitHub submission

#Update 2021-11-06 Samuel Taylor
# updated data files to be .csv to keep GitHub data submission tidy

#install.packages("ggplot2")
#install.packages("readxl")
#install.packages("dplyr")

library('ggplot2')
#library('readxl')
library('gridExtra')
library('dplyr')
library('patchwork')
library('here')
library('Rcpp')

#full colours
blck = rgb(0, 0, 0, maxColorValue = 255)
dgry = rgb(0.5*255, 0.5*255, 0.5*255, maxColorValue = 255)
dgrn = rgb(70, 190, 110, maxColorValue = 255)
lgrn = rgb(135, 240, 115, maxColorValue = 255)


#transparency colours

t_blck = rgb(0, 0, 0, maxColorValue = 255, alpha = 0.5*255)
t_dgry = rgb(0.5*255, 0.5*255, 0.5*255, maxColorValue = 255, alpha = 0.5*255)
t_dgrn = rgb(70, 190, 110, maxColorValue = 255, alpha = 0.5*255)
t_lgrn = rgb(135, 240, 115, maxColorValue = 255, alpha = 0.5*255)

#IT86D-1010 DF----
#Simplified Nov 2021 by Sam Taylor - supplying ony the relevant data as .csv
#Dragon_1010_4_Processed_DF = read_excel(
#  here("data/20200115_DragonsDen_RAN_4.xlsx"),
#  sheet = "DragonsDen_D_Processed"
#) #Processed DF
DF1010 = read.csv(here("data/Data_FigS5_IT86D1010.csv"))

##Figure 1A Annotated Final----
FigureS5A_Annotated = ggplot(DF1010, aes(x = Treatment, y = ActS)) +
  annotate(geom = "text",
           x = 887,
           y = 12.5,
           hjust = 0.445,
           size = 5,
           label = "IT86D-1010"
  ) +
  ggtitle("IT86D-1010 Original Light Response") +
  geom_point(color = t_dgrn, shape = 19, size = 4) +
  xlab("PPFD") + ylab("Activation State (%)") + theme_bw() +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200),
                     limits=c(0, 1250),
                     labels = c("0","200","400","600","800","1000","1200")
  ) +
  ylim(0,100) +
  stat_smooth(method='loess',
              color = dgrn,
              size = 2,
              span = 0.85,
              alpha = 0.2
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(color="black", size=14),
    axis.text = element_text(color = 'black', size=14),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = NA), #0.7
    legend.position = "none"
  )

FigureS5A_Annotated


FigureS5A_Annotated_NCS = FigureS5A_Annotated +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(vjust=-0.25, margin=margin(t=10)),
        axis.text.y = element_text(margin=margin(r=10))
  )
FigureS5A_Annotated_NCS

#ggsave(here("output/20200716_FigureS5A_NCS.tiff"),
#       FigureS5A_Annotated_NCS, width=6, height=4, dpi=600)

#IT82E-16 Activation Figure 1B
#Simplified Nov 2021 by Sam Taylor - supplying ony the relevant data as .csv
#Dragon_E16_4_Processed_DF = read_excel(here("data/20200115_DragonsDen_RAN_4.xlsx"),
#                                       sheet = "DragonsDen_E_Processed") #Processed DF

DF16 <- read.csv(here("data/Data_FigS5_IT82E16.csv"))
###Figure 1B

FigureS5B_Annotated = ggplot(DF16, aes(x = Treatment, y = ActS)) +
  annotate(geom = "text",
           x = 887,
           y = 12.5,
           hjust = 0.445,
           size = 5,
           label = "IT82E-16") +
  ggtitle("IT82E-16 Light Response") +
  geom_point(color = t_lgrn, shape = 19, size = 4) +
  xlab("PPFD") +
  ylab("Activation State (%)") +
  theme_bw() +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200),
                     limits=c(0, 1250),
                     labels = c("0","200","400","600","800","1000","1200")
  ) +
  ylim(0,100) +
  stat_smooth(method='loess',
              color = lgrn,
              size = 2,
              span = 0.945,
              alpha = 0.3) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(color="black", size=14),
    axis.text = element_text(color = 'black', size=14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = NA),
    legend.position = "none"
  )
FigureS5B_Annotated


FigureS5B_Annotated_NCS = FigureS5B_Annotated +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(vjust=-0.25, margin=margin(t=10)),
        axis.text.y = element_text(margin=margin(r=10))
  )

FigureS5B_Annotated_NCS

#ggsave("20200716_FigureS5B_NCS.tiff",
#       FigureS5B_Annotated_NCS, width=6, height=4, dpi=600) ##Wait to confirm measurements



##Figure 1C TVNu-1948
#Simplified Nov 2021 by Sam Taylor - supplying ony the relevant data as .csv
#Dragon_1948_4_Processed_DF = read_excel(here("data/20200115_DragonsDen_RAN_4.xlsx"),
#                                        sheet = "DragonsDen_T_Processed") #Processed DF
DF1948 <- read.csv(here("data/Data_FigS5_TVNu1948.csv"))


FigureS5C_Annotated = ggplot(DF1948, aes(x = Treatment, y = ActivationState)) +
  annotate(geom = "text",
           x = 887,
           y = 12.5,
           hjust = 0.445,
           size = 5,
           label = "italic('V. sp. Savi.')",
           parse = TRUE) + ##Species name change
  ggtitle("TVNu-1948 Light Response") +
  geom_point(color = t_blck, shape = 19, size = 4) +
  xlab("PPFD"~~(mu*mol~~m^-2~~s^-1)) +
  ylab("Activation State (%)") + theme_bw() +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200),
                     limits=c(0, 1250),
                     labels = c("0","200","400","600","800","1000","1200")
  ) +
  ylim(0,100) +
  stat_smooth(method='loess',
              color = blck,
              size = 2,
              span = 1.1628,
              alpha = 0.3) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(color="black", size=14),
    axis.text = element_text(color = 'black', size=14),
    axis.title.x = element_text(vjust=-1),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = NA),
    legend.position = "none"
  )
FigureS5C_Annotated


FigureS5C_Annotated_NCS = FigureS5C_Annotated +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(vjust=-0.25, margin=margin(t=10)),
        axis.text.y = element_text(margin=margin(r=10))
  )
FigureS5C_Annotated_NCS


#ggsave("20200716_FigureS5C_NCS.tiff",
#       FigureS5C_Annotated_NCS, width=6, height=4, dpi=600) ##Wait to confirm measurements

##Figure 1D

#Simplified Nov 2021 by Sam Taylor - supplying ony the relevant data as .csv
#Dragon_Aden_4_Processed_DF = read_excel(here("data/20200115_DragonsDen_RAN_4.xlsx"),
#                                        sheet = "DragonsDen_N_Processed") #Processed DF
DFAden <- read.csv(here("data/Data_FigS5_VAdenantha.csv"))

FigureS5D_Annotated = ggplot(DFAden, aes(x = Treatment, y = ActS)) +
  annotate(geom = "text",
           x = 887,
           y = 12.5,
           hjust = 0.44,
           size = 5,
           label = "italic('V. adenantha')",
           parse = TRUE) +
  ggtitle("Adenantha Light Response") +
  geom_point(color = t_dgry, shape = 19, size = 4) +
  xlab("PPFD"~~(mu*mol~~m^-2~~s^-1)) +
  ylab("Activation State (%)") + theme_bw() +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200),
                     limits = c(0, 1250),
                     labels = c("0","200","400","600","800","1000","1200")
  ) +
  ylim(0,100) +
  stat_smooth(method='loess',
              color = dgry,
              size = 2,
              span = 0.85,
              alpha = 0.3) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(color="black", size=14),
    axis.text = element_text(color = 'black', size=14),
    axis.title.x = element_text(vjust=-1),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.minor = element_line(size = NA),
    panel.grid.major = element_line(size = NA),
    legend.position = "none")
FigureS5D_Annotated

FigureS5D_Annotated_NCS = FigureS5D_Annotated +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(vjust=-0.25, margin=margin(t=10)),
        axis.text.y = element_text(margin=margin(r=10))
  )

FigureS5D_Annotated_NCS

#ggsave("20200716_FigureS5D_NCS.tiff",
#       FigureS5D_Annotated, width=6, height=4, dpi=600) ##Wait to confirm measurements


##Figure 1 Layout with patchwork

##Final Figure 1
FigureS5A_Annotated_NCS #1010
FigureS5B_Annotated_NCS #E16
FigureS5C_Annotated_NCS #1948
FigureS5D_Annotated_NCS #Adenantha


Layout_FigureS5_Annotated_NCS  = (FigureS5A_Annotated_NCS |
                                    FigureS5B_Annotated_NCS) /
  (FigureS5C_Annotated_NCS |
     FigureS5D_Annotated_NCS) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 14, face =  "bold"))

Layout_FigureS5_Annotated_NCS


ggsave(here("output/FigureS5.pdf"),
       Layout_FigureS5_Annotated_NCS,
       width=10, height=8
) ##Wait to confirm measurements
