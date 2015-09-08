setwd("~/Git/fossilDating")
createPhyloAgeVsPhyloAgePDF("compareAgeM1M8.pdf",df8,df1,"Model 8", "Model 1",c(),c())

pdf(file="comparePostM1M8.pdf", width=5, height=5, pointsize=10)
par(mar = c(4,4,1,1) + 0.1)
plot(df8$post, df1$post, xlab="Posterior probability in Model 8", ylab="Posterior probability in Model 1", pch=16, col="blue", xlim=c(0,1), ylim=c(0,1))
lm <- lm(df1$post ~ df8$post)
abline(lm, col="red")
r2 <- format(summary(lm)$r.squared,digits=3)
eq <- bquote(bold(R^2 == .(r2)))
    
text(0.05,0.925, labels=eq,pos=4)
dev.off()

logbf1 <- log(df1$bf)
logbf8 <- log(df8$bf)

pdf(file="compareBFM1M8.pdf", width=5, height=5, pointsize=10)
par(mar = c(4,4,1,1) + 0.1)
plot(logbf8, logbf1, xlab="log(BF) in Model 8", ylab="log(BF) in Model 1", pch=16, col="blue")
lm <- lm(logbf1 ~ logbf8)
abline(lm, col="red")
r2 <- format(summary(lm)$r.squared,digits=3)
eq <- bquote(bold(R^2 == .(r2)))
    
text(0.5,4.1, labels=eq,pos=4)
dev.off()

pdf(file="compareErrorM1M8.pdf", width=5, height=5, pointsize=10)
par(mar = c(4,4,1,1) + 0.1)
plot(df8$err, df1$err, xlab="Error in Model 8 age (My)", ylab="Error in Model 1 age (My)", pch=16, col="blue")
lm <- lm(df1$err ~ df8$err)
abline(lm, col="red")
r2 <- format(summary(lm)$r.squared,digits=3)
eq <- bquote(bold(R^2 == .(r2)))
    
text(0.5,15.5, labels=eq,pos=4)
dev.off()
