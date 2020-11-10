## Fetch UTR or promoter sequences (in R)

library("GenomicFeatures")
gtf_file="/home/test/motif/sequence_motif/genome/gencode.v27.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# get 3'UTR & 5'UTR site range
utr5p = fiveUTRsByTranscript(txdb, use.names=T)
utr3p = threeUTRsByTranscript(txdb, use.names=T)

utr3p.df=as.data.frame(utr3p)
utr5p.df=as.data.frame(utr5p)

write.table(utr3p.df, "utr3p.info", row.names=FALSE, sep='\t',quote=FALSE )
write.table(utr5p.df, "utr5p.info", row.names=FALSE, sep='\t' ,quote=FALSE)
# view them with cat, head, less, column, notepad, vscode, and so on

# get promoter site range
promoter=promoters(txdb)

promoter.df=as.data.frame(promoter)

write.table(promoter.df, "promoter.info", row.names=FALSE, sep='\t' ,quote=FALSE)

## Intersect with interested genes (in linux)

## Convert to bed format (in linux)

## Get genome sequence (firstly in linux then in R)

# concatenate sequences of the same 3’ UTR and the same promoter
library(dplyr)
concatenate_seq <- function(fasta_file) {
    biozhuoer::read_fasta(fasta_file) %>%
        dplyr::mutate(name = stringr::str_extract(name, 'ENST[\\d\\.]+')) %>% 
        dplyr::group_by(name) %>% dplyr::summarise(seq = paste0(seq, collapse = '')) %>% 
        biozhuoer::write_fasta(fasta_file)
}
concatenate_seq('interested_three_prime_UTR.fa')
concatenate_seq('interested_promoter.fa')

## Generate random sequence as background sequence (in linux except specified)

# downstream 1000bp as background (in R)
library(dplyr)

slide <- function(input_bed, output_bed, n = 1000) {
    col_names <- c('chr', 'start', 'end', 'name', 'score', 'strand');

    original <- readr::read_tsv(input_bed, col_names) %>% 
        dplyr::group_by_at(-2:-3) %>% 
        dplyr::summarise(length = sum(end - start), end = max(end)) %>% 
        dplyr::ungroup()

    if (n > 0) {
       slide <- original %>% dplyr::mutate(start = end + n, end = start + length)
    } else {
       slide <- original %>% dplyr::mutate(end = start + n, start = end - length)
    }

    slide %>% dplyr::select(chr, start, end, name, score, strand) %>% 
        readr::write_tsv(output_bed, col_names = F)
}
slide('interested_three_prime_UTR.bed', 'interested_three_prime_UTR_downstream.bed')
slide('interested_promoter.bed', 'interested_promoter_downstream.bed')

## Motif enrichment (in linux)