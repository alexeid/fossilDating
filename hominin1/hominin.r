require(coda)
source("~/Subversion/compevol/research/fossilDating/fossilDatingFunctions.r")

setwd("~/Subversion/compevol/research/fossilDating/hominin1")
fossilTable_h <- read.table("homininFossilTable.txt", sep="\t", header=T)
regexp = "\\-[0-9]+\\.log"
seq <- read.table("~/Subversion/compevol/research/fossilDating/hominin1/hominin_sorted_sequences.txt", header=T, sep="\t")
seq$Sequence <- as.character(seq$Sequence)

result_h <- loadFossilFiles(".", regexp, fossilTable_h, 5001)
height_h <- result_h$heights
df_h <- createSummaryTable(result_h$names, height_h, fossilTable_h, seq, hpdProb=0.99)

caption <- "Summary of results for 19 fossil hominins under Model 1. {\\em post} is the posterior probability that the phylogenetic age is within the paleaontological age range. {\\em BF} is the bayes factor in support of the palaeontogical age. {\\em phylo age} is the phylogenetic estimate of the age, along with the upper and lower of the corresponding 95\\% HPD credible interval. {\\em error} is the difference in millions of years between the phylogenetic point estimate of the fossil's age and the mean of it's paleaontological age range. {\\em ESS} is the estimated effective sample size for the phylogenetic age estimate."
label <- "fossilTable_h"

writeSummaryTable("1h_summaryTable.tex",caption, label, df_h)

ind <- 1:nrow(df_h)
ind <- ind[df_h$probInHPD == 00]
orient <- rep(4,length(ind))

createPhyloAgeVsGeoAgePDF("1h_phyloAgeVsGeoAge.pdf", df_h, fossilTable_h, ind, orient, min=0, max=8)
createPrecisionVsKnownCharactersFigure("1h_precisionVsKnownCharacters.pdf", df_h, ind, orient);

createHistPlotFigures("1h_fossilDatingHist", result_h$names, height_h, fossilTable_h, combinedPlots=TRUE, xmax=8, binsPerMY=2)