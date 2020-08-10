
package_lst <- c("ape", "reshape", "gridExtra", "ggplot2", "scales", "grid")

not_installed <- package_lst[!(package_lst %in%
                               installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

library(ape)
library(ggplot2)
library(reshape)
library(scales)
library(grid)
library(gridExtra)

# colors
alpha <- .7
c_yellow <-          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
c_blue <-            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
c_orange <-          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
c_green <-           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_very_dark_green <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_sea_green <-       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
c_black <-           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
c_grey <-            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_red <-        rgb(255 / 255, 130 / 255,   0 / 255, alpha)
c_very_light_red <- rgb(  255 / 255, 105 / 255, 108 / 255, alpha)
c_cyan3 <-          rgb(0 / 255,  205 / 255,   205 / 255, alpha)

# ggplot2 theme
main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line.x=element_line(color="black"),
              axis.line.y=element_line(color="black"),
              axis.ticks=element_line(color="black"),
              axis.text=element_text(colour="black", size=10),
              legend.position="top",
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(family="sans"))

