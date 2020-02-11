#' ---
#' title: "Class 5: Data Visualization and graphs in R"
#' author: "Cindy Tran"
#' date: "January 21st, 2020"
#' ---

# Class 5
# Data Visualization and graphs in R
plot(1:10, col="blue", typ= "o")

# Section 2 
# Now we will import a file
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

# Plot the weight data
plot(weight$Age, weight$Weight, typ="o")
# change the character to be a filled square
plot(weight$Age, weight$Weight, typ= "o", pch=15)
# change the point size to be 1.5x normal size
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5)
# change thickness to be twice the default size
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2)
# changes the y-axis limits to scale between 2 and 10 kg
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10))
# change the x-axis label to be Age (months)
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab ="Age (months)")
# change y axis label to Weight (kg)
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab ="Age (months)", ylab= "Weight (kg)")
# add a suitable title to the top of the plot
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab ="Age (months)", ylab= "Weight (kg)", main= "Baby weight with age")

# complete code of 2A
plot(weight$Age, weight$Weight, typ= "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab ="Age (months)", ylab= "Weight (kg)", main= "Baby weight with age")

# Enter in feature data for 2B
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header= TRUE, sep = "\t")

# make a barplot
barplot(mouse$Count)
# barplot with bars horizontal
barplot(mouse$Count, horiz = TRUE)
# set feature names 
barplot(mouse$Count, horiz = TRUE, names.arg= mouse$Feature, ylab ="", las=1)
# give it a title
barplot(mouse$Count, horiz = TRUE, names.arg= mouse$Feature, ylab ="", las=1, main="Numbers of features in the mouse GRCm38 genome")
# change x axis limits to 0 to 80000
barplot(mouse$Count, horiz = TRUE, names.arg= mouse$Feature, ylab ="", las=1, main="Numbers of features in the mouse GRCm38 genome", xlim = c(0,80000))
# adjust the margins so that the labels show in the plot on the right side
par(mar=c(3.1,11,4,2))
barplot(mouse$Count, horiz = TRUE, names.arg= mouse$Feature, ylab ="", las=1, main="Numbers of features in the mouse GRCm38 genome", xlim = c(0,80000))
# for the margins, it goes bottom, left, top, and right

# 2C
x <- c(rnorm(10000), rnorm(10000)+4)
par(mar=c(5,4,2,2))
hist(x, breaks=10)

# Section 3: Using color in plots
# call male female counts file
mfcounts <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
# plot as barplot
barplot(mfcounts$Count, names.arg = mfcounts$Sample, las=2, col=rainbow(10))
# using nrow to determine number of colors
barplot(mfcounts$Count, names.arg = mfcounts$Sample, las=2, col=rainbow(nrow(mfcounts)))
# replot to make the males and females different colors
# final graph
barplot(mfcounts$Count, names.arg = mfcounts$Sample, las=2, col=c("purple2", "pink2"), ylab = "Counts")

# Section 3B
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
# how many rows
nrow(genes)
# gives you a count of how many values in wach condition
table(genes$State)
# plot condition1 vs. condition2
plot(genes$Condition1, genes$Condition2, col=genes$State, ylab="Expression condition2", xlab="Expression condition 1")
# Make the colors red, blue and gray
palette(c("blue", "gray", "red"))
plot(genes$Condition1, genes$Condition2, col=genes$State, ylab="Expression condition2", xlab="Expression condition 1")

#Create a lab report
# Press the notepad at the top and then compile