require(coda)
source("~/Git/fossilDating/fossilDatingFunctions.r")

setwd("~/Git/fossilDating/8_penguins_eNprior")
fossilTable8 <- read.table("fossilTable.txt", sep="\t", header=T)
regexp = "\\-[0-9]+\\.log"
seq <- read.table("~/Git/fossilDating/sorted_sequences.txt", header=T, sep="\t")
seq$Sequence <- as.character(seq$Sequence)

h1prob <- read.table("~/Git/fossilDating/1_penguins_eNprior/H1prob-final.txt", sep=" ", header=T)
fossilTable8$prior <- h1prob$fossilH1

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

createHistPlotFigures("8_fossilDatingHist_younger", result8$names[df8$geo < 32], height8[df8$geo < 32], fossilTable8[df8$geo < 32,], combinedPlots=TRUE, xmin=0, xmax=45, numrow=6, numcol=3)

createHistPlotFigures("8_fossilDatingHist_older", result8$names[df8$geo > 32], height8[df8$geo > 32], fossilTable8[df8$geo > 32,], combinedPlots=TRUE, xmin=10, xmax=70, numrow=6, numcol=3)

createHistPlotFigures("8_fossilDatingHist_younger_wide", result8$names[df8$geo < 32], height8[df8$geo < 32], fossilTable8[df8$geo < 32,], combinedPlots=TRUE, xmin=0, xmax=45, numrow=3, numcol=6, pdfWidth=12, pdfHeight=6, histCol=rgb(0.1, 0.1, 0.9, 1.0))

createHistPlotFigures("8_fossilDatingHist_older_wide", result8$names[df8$geo > 32], height8[df8$geo > 32], fossilTable8[df8$geo > 32,], combinedPlots=TRUE, xmin=10, xmax=70, numrow=3, numcol=6, pdfWidth=12, pdfHeight=6, histCol=rgb(0.1, 0.1, 0.9, 1.0))