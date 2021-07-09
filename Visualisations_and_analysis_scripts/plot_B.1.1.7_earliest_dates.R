#Plots the cumulative distribution of adm2s with B.1.1.7 and calculates the number of adm2s added each week
#Run on adm2_earliest_sequence.csv from adm2_earliest_date.py

args <- commandArgs(TRUE)

library(ggplot2)

#Import the adm2 data
adm2s <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE)
#Convert dates to date format
adm2s$Date <- as.Date(adm2s$earliest_sequence)

outFile <- args[2]

#Calculate the number of UTLAs added each week
ewAdm2s <- data.frame(week_start = as.Date(character()), adm2s = numeric())
#Start date of the first week
startDate <- as.Date("2020-09-20")
endDate <- as.Date("2021-01-10")
#Incremented with each epi week
iterator <- 1
while (startDate < endDate) {
  ewAdm2s[iterator,1] <- startDate
  ewAdm2s[iterator,2] <- length(which(adm2s$Date <= (startDate + 6)))
  startDate <- startDate + 7
  iterator <- iterator + 1}

#Convert week starts to weeks to fit linear model
ewAdm2s$week <- 1:length(ewAdm2s[,1])

#Fit lm of epi week vs number of regions and calculate average regions added each week
linearMod <- lm(adm2s ~ week, data = ewAdm2s)
print(linearMod)

pdf(outFile)
uPlot <- ggplot(ewAdm2s, aes(x = week_start, y = adm2s)) + theme_bw() + scale_y_continuous(limits = c(0, 100), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + geom_rect(aes(xmin = as.Date("2020-11-05"), xmax = as.Date("2020-12-02"), ymin = 0, ymax = 100), fill = "turquoise4", colour = NA, alpha = 0.1) + geom_bar(stat = "identity", fill = "grey70") + xlab("Week start") + ylab("Number of adm2s with B.1.1.7 sequenced") + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))
print(uPlot)
dev.off()