
library(ggplot2)

dat.plot.fixed <- read.csv(“SGTF_growths_4b.csv”)

ggplot(dat.plot.fixed, aes(x = Vaccinated, y = Rate, group = Type, color = Type)) + geom_point() + xlab("2020-21") + ylab("Growth rate") + scale_color_manual(values = c("#4d4d4d", "#b2182b", "#2166ac", "#e08214")) + theme(legend.position = c(0.15,0.15), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + labs(color = "Variant")+scale_y_continuous() + geom_hline(yintercept = 0, linetype = "dashed")
