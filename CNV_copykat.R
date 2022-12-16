library(data.table)
library(copykat)

LFILES <- grep(".counts.txt", list.files(), value=T)

for (i in 1:length(LFILES)) {
    
    print(paste0("Analyzing ", LFILES[i]))
    
    exp.rawdata <- fread(LFILES[i], header=T, data.table=F)
    rownames(exp.rawdata) <- exp.rawdata[, 1]
    exp.rawdata <- as.matrix(exp.rawdata[, -1])
    sample <- sapply(strsplit(LFILES[i], "\\."), "[", 1)
    
    copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=10, win.size=25, KS.cut=0.1, 
                        sam.name=sample, distance="euclidean", norm.cell.names="", output.seg="FLASE", 
                        plot.genes="FALSE", genome="hg20", n.cores=16)
    
}
