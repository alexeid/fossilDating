require(coda)
source("~/Subversion/compevol/research/fossilDating/fossilDatingFunctions.r")

setwd("~/Subversion/compevol/research/fossilDating/run1_2")
fossilTable1 <- read.table("1_fossilTable.txt", sep="\t", header=T)
regexp = "\\-[0-9]+\\.log"
seq <- read.table("~/Subversion/compevol/research/fossilDating/sorted_sequences.txt", header=T, sep="\t")
seq$Sequence <- as.character(seq$Sequence)

result1 <- loadFossilFiles(".", regexp, fossilTable1, 10001)
height1 <- result1$heights
df1 <- createSummaryTable(result1$names, height1, fossilTable1, seq)

caption <- "Summary of results for 36 fossil penguins under Model 1. {\\em post} is the posterior probability that the phylogenetic age is within the paleaontological age range. {\\em BF} is the bayes factor in support of the palaeontogical age. {\\em phylo age} is the phylogenetic estimate of the age, along with the upper and lower of the corresponding 95\\% HPD credible interval. {\\em error} is the difference in millions of years between the phylogenetic point estimate of the fossil's age and the mean of it's paleaontological age range. {\\em ESS} is the estimated effective sample size for the phylogenetic age estimate."
label <- "fossilTable1"

writeSummaryTable("1_summaryTable.tex",caption, label, df1)

ind <- 1:nrow(df8)
ind <- ind[df8$post < 0.05]
orient <- rep(4,length(ind))
orient <- c(4,2,4,4,2)

createPhyloAgeVsGeoAgePDF("1_phyloAgeVsGeoAge.pdf", df1, fossilTable1, ind, orient)
createPrecisionVsKnownCharactersFigure("1_precisionVsKnownCharacters.pdf", df1, ind, orient);

createHistPlotFigures("1_fossilDatingHist", result1$names, height1, fossilTable1, combinedPlots=TRUE)