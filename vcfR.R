### vcfR: helping visualize, manipulate and quality filter data in VCF file

library(vcfR)

# read
vcf <- read.vcfR("D:/Rice-NJU/nitrate/data/3k-SNP/01_Nipponbare_B001.snp.vcf")

# access and process
vcf@meta
vcf@fix
vcf@gt

queryMETA(vcf, element = "AD")
getFIX(vcf)[1:2, 1:7]
vcf@gt[1:5, 1:7]
gt <- extract.gt(vcf, element = "GT")
gt[1:2, 1:7]



## Add Genome Info and Visualize
# Note: only one chromosome can be used!
library(vcfR)
library(ape)

# read data
vcf <- read.vcfR("D:/Rice-NJU/nitrate/data/3k-SNP/01_Nipponbare_B001.snp.vcf")
fa <- ape::read.dna("D:/Rice-NJU/data/Ref/IRGSP-1.0_genome.fasta", format = "fasta")
gff <-
  read.table(
    "D:/Rice-NJU/data/Ref/Oryza_sativa.IRGSP-1.0.49.gff3",
    sep = "\t",
    quote = ""
  )

# get info
chrom <- create.chromR(
  name = "Oryza sativa japonica WGS",
  vcf = vcf,
  seq = fa,
  ann = gff,
  verbose = TRUE
)

# plot
chrom
plot(chrom)
chromoqc(chrom)

# process and plot
pchrom <- proc.chromR(chrom, verbose = TRUE)
plot(pchrom)
chromoqc(pchrom)

