rm(list = ls() )
library(topr)
library(data.table)
library(dplyr)
library(fastman)

mydata <- fread("pQTL_PCSK9_GWAS.txt")
dim(mydata)
mydata[1:3, ]

mydata1 <- mydata %>% dplyr::filter(P < 0.05)
dim(mydata1)

max(-log10(mydata1$P))

mydata2 <- mydata1 %>% dplyr::filter(P > 0)
max(-log10(mydata2$P))

mydata3 <- mydata2 %>% rename(BP=POS, chr=CHR)
#png("ManhattanPlot_pQTL_PCSK9_GWAS.png", width = 10, height = 6, units = "in", res = 300)
manhattan(df = mydata3, 
          annotate= c(5e-12), 
          sign_thresh = c(5e-8, 1e-5),
          sign_thresh_color = "red", 
          sign_thresh_label_size=5, gene_label_fontface = "bold",
          ymax = 90, ymin = 1, scale = 0.8, max.overlaps = 50,
          build=37,
          title = c("pQTL_PCSK9 GWAS Manhattan Plot")) 
dev.off()

mydata$P[mydata$P == 0] 
# my Q-Q plot need a complete set of p-value
mydata <- fread("pQTL_PCSK9_GWAS.txt")
fastqq(mydata$P, 
       main="pQTL PCSK9 GWAS QQ Plot")

#GeneCard
mychr=1
mystart=55505221
myend=55530525

#cis- pQTL
#within 1MB of the gene
#with p_value < 5e-8.
mystart_with_buffer = mystart - 1000000 
myend_with_buffer = myend + 1000000

cis_data <- mydata3 %>% dplyr::filter(chr==1) %>% 
  dplyr::filter(BP > mystart_with_buffer) %>%
  dplyr::filter(BP < myend_with_buffer) %>% 
  dplyr::filter(P < 5e-8) %>%
  arrange(P)


dim(cis_data)

#the top 10 cis-pQTL for PCSK9 is
cis_data[1:10, ]


#trans pQTL

trans_data <- mydata3 %>% dplyr::filter(P < 1.7e-11) %>% 
  dplyr::filter(! SNP %in% cis_data$SNP) %>%
  arrange(P)
trans_data[1:10, ]

dim(trans_data)
#157 trans
#243 cis 

write.table(mydata2, "pQTL_PCSK9_GWAS_p005.txt", row.names=F, sep="\t", quote=F)
