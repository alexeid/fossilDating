require(coda)
source("~/Subversion/compevol/research/fossilDating/fossilDatingFunctions.r")

setwd("~/Subversion/compevol/research/fossilDating/run8_5")
fossilTable8 <- read.table("~/Subversion/compevol/research/fossilDating/run8_2/fossilTable.txt", sep="\t", header=T)
regexp = "\\-[0-9]+\\.log"
seq <- read.table("~/Subversion/compevol/research/fossilDating/sorted_sequences.txt", header=T, sep="\t")
seq$Sequence <- as.character(seq$Sequence)

result8 <- loadFossilFiles(".", regexp, fossilTable8, 10001)
height8 <- result8$heights
df8 <- createSummaryTable(result8$names, height8, fossilTable8, seq)

caption <- "Summary of results for 36 fossil penguins under Model 8. {\\em post} is the posterior probability that the phylogenetic age is within the paleaontological age range. {\\em BF} is the bayes factor in support of the palaeontogical age. {\\em phylo age} is the phylogenetic estimate of the age, along with the upper and lower of the corresponding 95\\% HPD credible interval. {\\em error} is the difference in millions of years between the phylogenetic point estimate of the fossil's age and the mean of it's paleaontological age range. {\\em ESS} is the estimated effective sample size for the phylogenetic age estimate."
label <- "fossilTable8"

writeSummaryTable("8_summaryTable.tex",caption, label, df8)

ind <- 1:nrow(df8)
ind <- ind[df8$post < 0.05]
orient <- rep(4,length(ind))
orient <- c(4,2,4,4,2)

createPhyloAgeVsGeoAgePDF("8_phyloAgeVsGeoAge.pdf", df8, fossilTable8, ind, orient)
createPrecisionVsKnownCharactersFigure("8_precisionVsKnownCharacters.pdf", df8, ind, orient);

# Need to split younger and older ones

createHistPlotFigures("8_fossilDatingHist_younger", result8$names[df8$geo < 32], height8[df8$geo < 32], fossilTable8[df8$geo < 32,], combinedPlots=TRUE, xmin=0, xmax=45)

createHistPlotFigures("8_fossilDatingHist_older", result8$names[df8$geo > 32], height8[df8$geo > 32], fossilTable8[df8$geo > 32,], combinedPlots=TRUE, xmin=10, xmax=70)