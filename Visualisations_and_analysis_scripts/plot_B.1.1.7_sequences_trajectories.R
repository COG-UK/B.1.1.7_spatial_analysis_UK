#Plots the trajectories of B.1.1.7 through time
#Run on adm2_weekly_sequences.csv

args <- commandArgs(TRUE)

library(ggplot2)
library(gridExtra)
library(cowplot)

#Import the lineage trajectories
trajectories <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE)
#Convert dates to date format
trajectories$Date <- as.Date(trajectories$week_start)

#Remove the final week
trajectories <- trajectories[which(trajectories$week_start <= as.Date("2021-01-10")),]

#Import adm2 colours and extract
adm2Colours <- read.csv(args[2], header = TRUE, stringsAsFactors = FALSE)
colours <- c()
for (i in adm2Colours$colour) {
  colours <- c(colours, i)}
names(colours) <- adm2Colours$adm2

#Extract the unique adm2s
adm2s <- unique(trajectories$adm2)
#Extract and sort the unique weeks
weeks <- unique(trajectories$week_start)

#Will be filled with the cumulative number of B.1.1.7 sequences in each adm2
adm2Trajectory <- data.frame(adm2 = as.character(), Week = as.character(), Sequences = as.numeric(), stringsAsFactors = FALSE)

#Incremented with each new row for adm2Trajectory
iterator <- 1

#Iterate through the weeks in each adm2, calculate the cumulative number of sequences and add to adm2Trajectory
for (adm2 in adm2s) {
  adm2Sequences <- trajectories[which(trajectories$adm2 == adm2),]
  
  for (week in weeks) {
    wDate <- as.Date(week)
    adm2Trajectory[iterator,1] <- adm2
    adm2Trajectory[iterator,2] <- week
    adm2Trajectory[iterator,3] <- sum(adm2Sequences$B.1.1.7_sequences[which(adm2Sequences$Date <= wDate)])
    iterator <- iterator + 1}}

  #pdf(paste(adm2, ".pdf", sep = ""))
  #proportionPlot <- ggplot(adm2Sequences, aes(x = Date, y = Proportion_sequences_B.1.1.7)) + theme_bw() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) + geom_bar(stat = "identity")
  #print(proportionPlot)
  #dev.off()}

#Convert dates to date format
adm2Trajectory$Date <- as.Date(adm2Trajectory$Week)

pdf("adm2_trajectories.pdf")
adm2Plot <- ggplot(adm2Trajectory, aes(x = Date, y = Sequences, group = adm2, colour = adm2)) + theme_bw() + geom_line() + scale_colour_manual(values = colours) + guides(colour = FALSE) + ylab("Cumulative B.1.1.7 sequences") + theme(axis.title.y = element_text(hjust = -0.2))
print(adm2Plot)
dev.off()

#Plot proportions through time for specific areas
kSequences <- trajectories[which(trajectories$adm2 == "KENT"),]
kPlot <- ggplot(kSequences, aes(x = Date, y = proportion_sequences_B.1.1.7)) + theme_bw() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + geom_bar(stat = "identity", fill = "dodgerblue") + ylab("Proportion of sequences B.1.1.7") + ggtitle("Kent") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank())
lSequences <- trajectories[which(trajectories$adm2 == "GREATER_LONDON"),]
lPlot <- ggplot(lSequences, aes(x = Date, y = proportion_sequences_B.1.1.7)) + theme_bw() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + geom_bar(stat = "identity", fill = "red") + ylab("Proportion of sequences B.1.1.7") + ggtitle("Greater London") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank())
eSequences <- trajectories[which(trajectories$adm2 == "ESSEX"),]
ePlot <- ggplot(eSequences, aes(x = Date, y = proportion_sequences_B.1.1.7)) + theme_bw() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + geom_bar(stat = "identity", fill = "orange") + ylab("Proportion of sequences B.1.1.7") + ggtitle("Essex") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank())
nSequences <- trajectories[which(trajectories$adm2 == "NORFOLK"),]
nPlot <- ggplot(nSequences, aes(x = Date, y = proportion_sequences_B.1.1.7)) + theme_bw() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + geom_bar(stat = "identity", fill = "purple") + ylab("Proportion of sequences B.1.1.7") + ggtitle("Norfolk") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank())

regionPlot <- grid.arrange(arrangeGrob(kPlot, lPlot, ePlot, nPlot, ncol = 2, left = "Proportion of sequences B.1.1.7", bottom = "Date"))

pdf("combined_plot.pdf")
combinedPlot <- plot_grid(adm2Plot, regionPlot, ncol = 1, rel_heights = c(0.4, 0.6), labels = c("a.", "b."))
print(combinedPlot)
dev.off()


#print(kSequences)
