#!/usr/bin/env Rscript
library(MSCquartets)

fileName = "../data/tree_ord.txt"
gene_trees = read.tree(fileName)
QT = quartetTable(gene_trees, taxonnames = NULL, epsilon = 0, random = 0)
QTR = quartetTableResolved(QT, omit = F)

MSC = quartetTreeTestInd(QTR, model = "T3")
Star = quartetStarTestInd(QTR)

CFpVals = cbind(MSC, Star[, "p_star"])
outFileName = "msquartets"

write.csv(CFpVals, paste(outFileName, sep = "")) 

sum = 0
for (i in 6:100){
  sum =sum+ choose(99, x-1)*(i/100)^(100-k)*(1-i/100)^k
}