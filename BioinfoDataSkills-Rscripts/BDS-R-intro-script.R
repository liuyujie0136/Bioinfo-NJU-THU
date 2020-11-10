# Selected basic R scripts in book Bioinformatics Data Skills (chapter8).

## exploratory data analysis (EDA)
print(sqrt(3.5), digits=10) 
# overall 10 digits

round(sqrt(3.5), digits=3) 
# after . have 3 digits

seq(3,5)	#3:5
# sequence from 3 to 5

c(1, 2) + c(0, 0, 0, 0)
# R adds a shorter vector c(1, 2) to a longer vector c(0, 0, 0, 0) by recycling the shorter values.

x[1]
# R's vectors are 1-indexed, meaning that the index 1 corresponds to the first element in a list (in contrast to 0-indexed languages like Python). 

b <- c(a=3.4, b=5.4, c=0.4)
# Vectors can also have names
names(b) <- c("x", "y", "z")
# change these names 

z[c(-4, -5)]
# exclude fourth and fifth elements

z[order(z)]
# order() returns a vector of indexes that indicate the (ascending) order of the elements
z[order(z, decreasing=TRUE)] 

v[v > 2 & v < 4] 
# comparision for vectors

q <- c(2, 3.5, -1.1, 3.8)
typeof(q)
# what is the variable's type?

chr_hits <- c("chr2", "chr2", "chr3", "chrX", "chr2", "chr3", "chr3") 
hits <- factor(chr_hits, levels=c("chrX", "chrY", "chr2", "chr3", "chr4")) 
# create factor with levels

levels(hits) <- list(chrX="chrX", chrY="chrY", chr2="chr2", chr3="chr3", chr4="chr4") 
# for existing vector, & we use a named character vector to provide a mapping between the original names and the new names
table(hits) 
summary(hits)
# show levels
class(hits) 
# what is the variable's class?

nums <- c(0.97, -0.7, 0.44, 0.25, -1.38, 0.08)
summary(nums)
# show basic statistic values

d$percent.GC[d$Pi > 16]
# part of data.frame extracted as a vector can also do check on other cols in origin

which(d$Pi > 3) 
# return position
d[which.min(d$total.Bases),]
# find specified row according to returned position

subset(d, Pi > 16 & percent.GC > 80)
# create a subset


## ggplot2
ggplot(d) + geom_point(aes(x=position, y=diversity, color=cent))  
# position, diversity, and cent are columns in data.frame d

ggplot(d) + geom_density(aes(x=diversity), fill="black")
ggplot(d) + geom_density(aes(x=diversity, fill=cent), alpha=0.4)
# calculates a density from column diversity. cent column's content is true or false and can be used to classify into two density lines

ggplot(d, aes(x=depth, y=total.SNPs)) + geom_point() + geom_smooth()
# add a smoothing line to plots and look for an unexpected trend

d$GC.binned <- cut(d$percent.GC, 5) 
# divides the data into 5 equally sized bins.
ggplot(d) + geom_density(aes(x=depth, linetype=GC.binned), alpha=0.5) 
#  cut() is useful in grouping data

ggplot(d) + geom_bar(aes(x=percent.GC))
# bar plot or histogram. can also x=GC.binned


## data handling
c(3, 4, -1) %in% c(1, 3, 4, 8)
table(mtfs$pos %in% rpts$pos) 
# whether some of a vector's values are in another vector
match(c("A", "C", "E", "A"), c("A", "B", "A", "E"))
# return matched first position in vector 2

mtfs$pos <- paste(mtfs$chr, mtfs$motif_start, sep="-")
# merge two strings together with separator "-"

mtfs$repeat_name <- rpts$name[match(mtfs$pos, rpts$pos)] 
recm <- merge(mtfs, rpts, by.x="pos", by.y="pos") 
# merge two dataframes according to same position

unique(mtfs$motif) 
# show unique items or levels

ggplot(mtfs, aes(x=dist, y=recom)) + geom_point(size=1, color="grey") + geom_smooth(method='loess', se=FALSE, span=1/16) + facet_wrap(~ motif)
# facet_wrap() takes a factor column and creates a panel for each level and wraps around horizontally

ggplot(mtfs, aes(x=dist, y=recom)) + geom_point(size=1, color="grey") + geom_smooth(method='loess', se=FALSE, span=1/16) + facet_grid(repeat_name ~ motif)
# facet_grid() allows finer control of facets by allowing you to specify the columns to use for vertical and horizontal facets, i.e. there are many plots in one picture with different colnames in different corner

adh <- list(chr="2L", start=14615555L, end=14618902L, name="Adh") 
# list has different types of data
adh[1:2] #sublist
adh[[2]] #extract data
adh$chr #extract data (also in dataframes)
adh$chr <- NULL #remove data

ll <- list(a=rnorm(6, mean=1), b=rnorm(6, mean=4), c=rnorm(6, mean=6)) 
lapply(ll, mean) 
# calculate mean for each vector stored in this list; lapply(): the 'l' is for list, as lapply() returns the result as a list
ll$a[3] <- NA
lapply(ll, mean, na.rm=TRUE) 
# ignore NA values

library(parallel)
results <- mclapply(my_samples, slowFunction) 
# mclapply runs the function 'slowFunction' on each of the elements of 'my_samples' in parallel to enhance effciency. but often, efficient R code will lead to sizable performance gains without requiring parallelization.

sapply(ll, function(x) mean(x, na.rm=TRUE))
# similar to lapply(), except that it simplifies the results into a vector, array, or matrix 
mapply(function(a, b) length(intersect(a, b)), ind_1, ind_2, SIMPLIFY=FALSE) 
# a multivariate version of sapply(): the function you pass to mapply() can take in and use multiple arguments.

fun_name <- function(args) {
   # body, containing R expressions
   return(value) 
}
# create R functions


## split-apply-combine strategy
d_split <- split(d$depth, d$GC.binned)
# We split a dataframe or vector using split(x, f), where x is a dataframe/vector and f is a factor(grouping levels, can be more than one), returns a list(can use lapply/sapply to use functions on each item).
lapply(d_split, mean) 
# then apply function mean to each splitted group

cbind() #column bind
rbind() #row bind
# combine a list of vectors by binding each element together into a matrix or dataframe.
do.call(rbind, lapply(split(d$depth, d$GC.binned), summary))
# do function rbind until no data from lapply

tapply(d$depth, d$GC.binned, mean) 
aggregate(d$depth, list(gc=d$GC.binned), mean)
# advanced options, combines split and lapply together


## dplyr: arrange(), filter(), mutate(), select(), and summarize()
d_df <- tbl_df(d) 
# like head()

select(d_df, start, end, Pi, Recombination, depth)  #equivalent to d[, c("start", "end", "Pi", "Recombination", "depth")]
select(d_df, start:total.Bases)
select(d_df, -(start:cent)) 
# selecting some columns

filter(d_df, Pi > 16, percent.GC > 80) 
# select specific rows

arrange(d_df, depth)  #like d[order(d$depth), ]
arrange(d_df, desc(total.SNPs), desc(depth))  #in descending order 
# sorting by columns

d_df <- mutate(d_df, diversity = Pi/(10*1000)) 
# add new columns to our dataframe
# note: dplyr can use %>% (known as pipe) from the magrittr package (like | in linux)

mtfs_df %>% group_by(chr) %>% summarize(max_recom = max(recom), mean_recom = mean(recom), num=n())
# group and apply functions to each group


## string
nchar(c("AGCTAG", "ATA", "GATCTGAG", ""))  #returns: 6 3 8 0
# length("AGCTAG") will return 1

re_sites <- c("CTGCAG", "CGATCG", "CAGCTG", "CCCACA") 
grep("CAG", re_sites) 
grep("CT[CG]", re_sites) 
# grep function with POSIX extended regular expressions

chrs <- c("chrom6", "chr2", "chr6", "chr4", "chr1", "chr16", " chrom8")
regexpr("[^\\d]6", chrs, perl=TRUE)  #returns: 5 -1  3 -1 -1 -1 -1 ...(start from 0)
# returns where in each element of x it matched pattern; perl=TRUE means Perl Compatible Regular Expressions (PCRE) is supported

substr(x, start, stop) 
sub(pattern, replacement, x) 
# can be used with regular expression (further study! or perl)

paste("chr", c(1:22, "X", "Y"), sep="")
paste0("chr","om",1:22)  # same as paste("chr","om",1:22,sep="")
# paste and paste0

region <- "chr10:158395-172881"
chunks <- sub("(chr[\\d+MYX]+):(\\d+)-(\\d+)", "\\1;;\\2;;\\3", region, perl=TRUE)
strsplit(chunks, ";;") 
# extracting multiple values (see more in linux part)


## control flow: if, for, and while 
if (x == some_value) {
  # do some stuff in here
} else {
  # else is optional
}

for (element in some_vector) {  # e.g. i in 1:100
  # iteration happens here
}

while (something_is_true) {
  # do some stuff
} 
# basic syntax

ifelse(x < 0, -1, 1)  # ifelse(test, yes, no) 
# vector version of if


## loading and exporting data
read.table(x, header=FALSE, col.names=bedcols)
write.table(mtfs, file="hotspot_motifs.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
save(tmp, file="example.Rdata")
load("example.Rdata") 

