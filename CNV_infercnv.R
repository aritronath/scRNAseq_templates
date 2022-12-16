library(infercnv)
library(data.table)

args = commandArgs(trailingOnly = TRUE) 
sample <- as.character(args[1])
Ncores <- as.numeric(args[2])


dat <- fread(paste(sample, ".ic.counts.txt", sep=""), header = T, data.table=F, nThread=Ncores)
annot.file <- paste(sample, ".ic.annot.txt", sep="")

rownames(dat) <- dat[,1]
dat <- as.matrix(dat[,2:ncol(dat)])


# Rename non-malignant cells as normal and save 
annot <- read.table(annot.file, header=F, sep="\t", stringsAsFactors = FALSE)
x <- which(annot$V2 != "malignant")
annot$V2[x] <- "normal"
fwrite(annot, file=annot.file, quote=F, row.names=F, col.names=F, sep="\t") 


out.dir <- paste(sample, "output", sep="_")
dir.create(out.dir)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = dat, 
                                    annotations_file = annot.file, 
                                    delim = "\t", 
                                    ref_group_names = c("normal"),
                                    gene_order_file = "hg19.RefSeq.NM_pos_unique_sort.txt")

infercnv_obj = infercnv::run(infercnv_obj, num_threads=Ncores, resume_mode=TRUE, cutoff=0.1,
                             
                             #### Clustering paramaters
                             analysis_mode="subclusters",
                             tumor_subcluster_partition_method="qnorm",
                             tumor_subcluster_pval = 0.5,
                             cluster_by_groups=FALSE, #FALSE if submitting individual samples 
                             cluster_references=FALSE,
			     #k_obs_groups=2, #hopefully, splits malignant and non malignant cells	

                             #### Denoising
                             denoise=T,
                             sd_amplifier=1.5,  # sets midpoint for logistic
                             #noise_logistic=TRUE, # turns gradient filtering on
                             
                             #### CNV calls
                             HMM=TRUE,
                             HMM_type="i6",
                             BayesMaxPNormal=0.2, #Default is 0.5. Setting lower could result in cleaner profiles 

                             #### Output
                             out_dir=out.dir,
			     plot_steps=FALSE,
                             no_plot=FALSE)
                             
