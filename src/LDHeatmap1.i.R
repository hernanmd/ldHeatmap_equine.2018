# Description: LDheatmap interactive version
# Author : <hernan.morales@gmail.com>
# Input files:
#   .ld (inter-variant correlation table or matrix). Produced by PLINK "--r2 square".
# Output:
#   LD heat map plot

source("https://bioconductor.org/biocLite.R")
install.packages("LDheatmap")
library(LDheatmap)
# Install latest dev version:
# devtools::install_github("mcneney/LDheatmap")
# Test LDheatmap: LDheatmap(matrix(runif(100,0,1),nrow = 10))

###########################################################
#
# Input settings 
#
############################################################

getwd()
# setwd("c:\\")
ldSqFilename <- "EQN65kCAR-GENO01-chr32sq_10kb-90kb.ld"
mapFilename <- "EQN65kCAR-GENO01-chr32_10kb-90kb.map"

###########################################################
#
# Begin code 
#
############################################################

# Read Square LD matrix
ldDataFrame <- read.csv(ldSqFilename, header = FALSE, sep = "\t")
# Print data frame dimensions
dim(ldDataFrame)
# View(ldDataFrame)

# Read SNP Names
mapDataFrame <- read.csv(mapFilename, header = FALSE, sep = "\t")
View(mapDataFrame)
snpNames <- as.vector.factor(mapDataFrame[[2]])

ldMatrix <- as.matrix(ldDataFrame)
rm(ldDataFrame)
# View(ldMatrix)
class(ldMatrix)

###########################################################
#
# Plot LDHeatmap
#
############################################################

LDheatmap(ldMatrix, SNP.name = snpNames, color = "blueToRed")
?LDheatmap


###########################################################
#
# Edit grob 
#
############################################################

require(grid)
ldHMGrob <- LDheatmap(ldMatrix, SNP.name = snpNames, color = "blueToRed")

getNames()
# Find the names of the component grobs of "ldheatmap"
childNames(grid.get("ldheatmap"))
#[1] "heatMap" "geneMap" "Key"
#Find the names of the component grobs of heatMap
childNames(grid.get("heatMap"))
#[1] "heatmap" "title"
#Find the names of the component grobs of geneMap
childNames(grid.get("geneMap"))
#[1] "diagonal" "segments" "title"    "symbols"  "SNPnames"
#Find the names of the component grobs of Key
childNames(grid.get("Key"))
#[1] "colorKey" "title"    "labels"   "ticks"    "box"

###########################################################
#
# Edit grob : Not tested
#
############################################################

#Change the plotting symbols that identify SNPs rs2283092 and rs6979287
?grid.edit
#on the plot to bullets
grid.edit("symbols", pch = 20, gp = gpar(cex = 1))
#Change the color of the main title
grid.edit(gPath("ldheatmap", "heatMap", "title"), gp = gpar(col = "red"))
#Change size of SNP labels
grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=1.5))
#Add a grid of white lines to the plot to separate pairwise LD measures
grid.edit(gPath("ldheatmap", "heatMap", "heatmap"), gp = gpar(col = "white",
                                                              lwd = 2))
#### Modify a heat map using 'editGrob' function ####
MyHeatmap <- LDheatmap(MyHeatmap, color = grey.colors(20))
new.grob <- editGrob(MyHeatmap$LDheatmapGrob, gPath("geneMap", "segments"),
                     gp=gpar(col="orange"))
##Clear the old graphics object from the display before drawing the modified heat map:
grid.newpage()
grid.draw(new.grob)
# now the colour of line segments connecting the SNP
# positions to the LD heat map has been changed from black to orange.
#### Draw a resized heat map (in a 'blue-to-red' color scale ####
grid.newpage()
pushViewport(viewport(width=0.5, height=0.5))
LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"), newpage=FALSE,
          color="blueToRed")
popViewport()
#### Draw and modify two heat maps on one plot ####
grid.newpage()
##Draw and the first heat map on the left half of the graphics device
pushViewport(viewport(x=0, width=0.5, just="left"))
LD1<-LDheatmap(MyHeatmap, color=grey.colors(20), newpage=FALSE,
               title="Pairwise LD in grey.colors(20)",
               SNP.name="rs6979572", geneMapLabelX=0.6,
               geneMapLabelY=0.4, name="ld1")
upViewport()
##Draw the second heat map on the right half of the graphics device
pushViewport(viewport(x=1,width=0.5,just="right"))
LD2<-LDheatmap(MyHeatmap, newpage=FALSE, title="Pairwise LD in heat.colors(20)",
               SNP.name="rs6979572", geneMapLabelX=0.6, geneMapLabelY=0.4, name="ld2")
upViewport()
##Modify the text size of main title of the first heat map.
grid.edit(gPath("ld1", "heatMap","title"), gp=gpar(cex=1.5))
##Modify the text size and color of the SNP label of the second heat map.
grid.edit(gPath("ld2", "geneMap","SNPnames"), gp=gpar(cex=1.5, col="DarkRed"))
#### Draw a lattice-like plot with heat maps in panels ####
# Load CHBJPTSNP and CHBJPTDist
data(CHBJPTData)
# Make a variable which indicates Chinese vs. Japanese
pop <- factor(c(rep("chinese",45), rep("japanese",45)))
require(lattice)
xyplot(1:nrow(CHBJPTSNP) ~ 1:nrow(CHBJPTSNP) | pop,
       type="n", scales=list(draw=FALSE), xlab="", ylab="",
       panel=function(x, y, subscripts,...) {
         LDheatmap(CHBJPTSNP[subscripts,], CHBJPTDist, newpage=FALSE) })
data(GIMAP5)
require(lattice)
n<-nrow(GIMAP5$snp.data)
xyplot(1:n ~ 1:n | GIMAP5$subject.support$pop,
       type="n", scales=list(draw=FALSE), xlab="", ylab="",
       panel=function(x, y, subscripts,...) {
         LDheatmap(GIMAP5$snp.data[subscripts,],
                   GIMAP5$snp.support$Position, SNP.name="rs6598", newpage=FALSE) })
#Reset the user's setting for prompting on the graphics output
#to the original value before running these example commands.
devAskNewPage(old.prompt)

