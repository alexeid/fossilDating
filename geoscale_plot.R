setwd("~/Git/fossilDating")

library("phytools")
library("phyloch")
library("strap")
library("coda")

## user epochs
Name <- c("Ps","Po","Miocene","Oligocene","Eocene","Paleocene")
Start <- c(2.58,5.333,23.03,33.9,56.0,66)
End <- c(0.0117,2.58,5.333,23.03,33.9,56.0)
user_ep <- data.frame(Start, End, Name)

t <- read.beast("hominin1r_cog/1r_hominin_ranges_cog-1440423074150-mcc.tree")
t$root.time <- t$height[1]
#t$posterior[50] <- 0.0958
#t$height[50] <- 37.5857
#t$"height_95%_HPD_MIN"[50] <- 34.2724
#t$"height_95%_HPD_MAX"[50] <- 40.0015
#t$posterior[21] <- 0.154
#t$height[21] <- 15.6909
#t$"height_95%_HPD_MIN"[21] <- 12.824
#t$"height_95%_HPD_MAX"[21] <- 18.6343

#t$posterior[22] <- 0.481
#t$height[22] <- 17.6519
#t$"height_95%_HPD_MIN"[22] <- 16.0011
#t$"height_95%_HPD_MAX"[22] <- 19.7024

tre <- read.nexus("hominin1r_cog/1r_hominin_ranges_cog-1440423074150-mcc.tree")

num_taxa <- length(t$tip.label)

display_all_node_bars <- TRUE

names_list <-vector()
for (name in t$tip){
  v <- strsplit(name,"_")[[1]]
  if(display_all_node_bars){
    names_list = c(names_list,name)
  }
  else if(v[length(v)]=="0"){
    names_list = c(names_list,name)
  }
}

nids <- vector()
pos <- 1
len_nl <- length(names_list)
for(n in names_list){
  for(nn in names_list[pos:len_nl]){
    if(n != nn){
      m <- getMRCA(t,c(n,nn))
      if(m %in% nids == FALSE){
        nids <- c(nids,m)
      }
    }
  }
  pos<-pos+1
}

root_max <- t$"height_95%_HPD_MAX"[1]
x_max <- root_max #origin_HPD[2] * 0.1 + origin_HPD[2]

pdf("geoscaled.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(t,right=FALSE), boxes="Epoch", units=c("Period", "Epoch"), cex.tip=0.95,cex.age=1, quat.rm=TRUE,
             cex.ts=1,label.offset=0.25,x.lim=c(-1,x_max),width=1.5, erotate=0,urotate=0)
mtext("Ps",side=1,line=0.9,adj=0.77,cex=0.8)
mtext("Quat.",side=1,line=3,adj=0.78,cex=1)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.3),lwd=5)
}

t$node.label<-t$posterior
p <- character(length(t$node.label))
p[is.na(t$node.label)] <- "red"
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=1, bg=p)
dev.off()