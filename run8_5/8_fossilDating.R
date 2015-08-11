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

caption <- "Table of 36 fossil penguins. $post$ is the posterior probability that the phylogenetic age is within the paleaontological age range. Error is the difference in millions of years between the phylogenetic point estimate of the fossil's age and the mean of it's paleaontological age range. {\\em Age} is the phylogenetic estimate of the age, along with the upper and lower of the corresponding 95\\% HPD credible interval."
label <- "fossilTable"

writeSummaryTable("8_summaryTable.tex",caption, label, df8)

ind <- 1:nrow(df8)
ind <- ind[df8$post < 0.05]
orient <- rep(4,length(ind))

createPhyloAgeVsGeoAgePDF("8_phyloAgeVsGeoAge.pdf", df8, fossilTable8, ind, orient)
createPrecisionVsKnownCharactersFigure("8_precisionVsKnownCharacters.pdf", df8, ind, orient);

createDensityPlotFigures("fossilDatingDensities", result8$names, height8, fossilTable8, combinedPlots=TRUE)