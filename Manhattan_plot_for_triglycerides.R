rm(list = ls() )
library(topr)
library(dplyr)
library(data.table)
library(fastman)

# Load the GWAS dataset
mydata <- fread("Triglycerides_GWAS.txt")
dim(mydata)

# Preview data
mydata[1:3, ]

# Rename columns to match fastman format
mydata <- mydata %>% rename(POS=BP, CHROM=CHR)
max(-log10(mydata$P))
mydata <- mydata %>% dplyr::filter(P< 0.05)
dim(mydata)
mydata <- mydata %>% dplyr::filter(P> 0)
# graphics.off()

# Manhattan Plot with adjusted annotation for TG GWAS
mydata <- mydata %>% dplyr::filter(P > 2.2e-308)
manhattan(df = mydata,
          annotate = c(1e-150),  # Show fewer genes for TG
          sign_thresh = c(5e-8, 1e-5),
          sign_thresh_color = "red",
          sign_thresh_label_size = 5, gene_label_fontface = "bold",
          ymax = 330, ymin = 1, scale = 0.8, max.overlaps = 50,
          build = 37,
          title = "Triglycerides GWAS Manhattan Plot")

# my Q-Q plot need a complete set of p-value
mydata <- fread("Triglycerides_GWAS.txt")
fastqq(mydata$P, 
   main="Triglycerides GWAS Q-Q Plot")

