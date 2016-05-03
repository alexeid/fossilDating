require(coda)
source("~/Git/fossilDating/fossilDatingFunctions.r")

setwd("~/Git/fossilDating")
fossilTable_c <- read.table("~/Git/fossilDating/canids/stratRanges/canidFossilTable_stratRanges.txt", sep="\t", header=T)
regexp = "\\-[0-9]+\\.log"
seq <- read.table("~/Git/fossilDating/canids/stratRanges/canids-sorted-sequences.txt", header=T, sep="\t")
seq$Sequence <- as.character(seq$Sequence)

result_c <- loadFossilFiles(".", regexp, fossilTable_c, 10001)
height_c <- result_c$heights
df_c <- createSummaryTable(result_c$names, height_c, fossilTable_c, seq, hpdProb=0.95)

caption <- "Summary of results for 123 fossil canids under Model 1. {\\em post} is the posterior probability that the phylogenetic age is within the paleaontological age range. {\\em BF} is the bayes factor in support of the palaeontogical age. {\\em phylo age} is the phylogenetic estimate of the age, along with the upper and lower of the corresponding 95\\% HPD credible interval. {\\em error} is the difference in millions of years between the phylogenetic point estimate of the fossil's age and the mean of it's paleaontological age range. {\\em ESS} is the estimated effective sample size for the phylogenetic age estimate."
label <- "fossilTable_c"

writeSummaryTable("1c_summaryTable.tex",caption, label, df_c)

ind <- 1:nrow(df_c)
ind <- ind[df_c$probInHPD == 0]
orient <- rep(4,length(ind))

createPhyloAgeVsGeoAgePDF("1c_phyloAgeVsGeoAge.pdf", df_c, fossilTable_c, ind, orient, min=0, max=38,xlabel="stratigraphic age range (Myr)")
createPrecisionVsKnownCharactersFigure("1c_precisionVsKnownCharacters.pdf", df_c, ind, orient, relativeHPD=T);

createHistPlotFigures("1c_fossilDatingHist_0to3", result_c$names[df_c$geo < 3], height_c[df_c$geo < 3], fossilTable_c[df_c$geo < 3,], combinedPlots=TRUE, xmin=0, xmax=10)

cond <- (df_c$geo >= 3 & df_c$geo < 10)

createHistPlotFigures("1c_fossilDatingHist_3to10", result_c$names[cond], height_c[cond], fossilTable_c[cond,], combinedPlots=TRUE, xmin=0, xmax=20)