#The file gene_expression.tsv downloaded from the github repository.
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "try.tsv")

#Read the file with gene accession numbers as row numbers
a <- read.table("try.tsv", header = TRUE, row.names = 1)

#Displaying the values of the first six genes
head(a)

#Making a new column which contains the mean of other columns
a$Mean <- rowMeans(a)

#Displaying the values of the first six genes
head(a)

#Listing the 10 genes with highest mean expression
list <- a[order(-a$Mean),]
head(list,10)

#Number of genes with a mean>10
nrow( subset(a, a$Mean<10))

#histogram plot of the means
a$Mean <- as.matrix(a)
range(a$Mean)
hist(a$Mean)
hist(as.matrix(a$Mean),10, xlab = "Mean", breaks = 50, col = "blue", xlim = c(0,75000))


#download file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",
              destfile = "try2.csv")
y <- read.table("try2.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
colnames(y)

#subset
y[1:50,]
northeast <- y[1:50, ]
#calculate mean and sd
mean(northeast$Circumf_2004_cm)
sd(northeast$Circumf_2004_cm)
mean(northeast$Circumf_2019_cm)
sd(northeast$Circumf_2019_cm)
#subset
y[51:100,]
southwest <- y[51:100, ]
#calculate mean and sd
mean(southwest$Circumf_2004_cm)
sd(southwest$Circumf_2004_cm)
mean(southwest$Circumf_2019_cm)
sd(southwest$Circumf_2019_cm)

#boxplot
boxplot(southwest$Circumf_2004_cm,southwest$Circumf_2019_cm,northeast$Circumf_2004_cm,northeast$Circumf_2019_cm,
        names = c("SW2004","SW2019","NE2004","NE2019"),ylab="Cirumference (cm)",
        main="Growth at Two Plantation Sites")

#mean growth for 10 years

GrowthSW <- (southwest$Circumf_2019_cm-southwest$Circumf_2009_cm)     
GrowthNE <- (northeast$Circumf_2019_cm-northeast$Circumf_2009_cm)
mean(GrowthSW)
mean(GrowthNE)

head(y)

#t.test
res <- t.test(GrowthSW,GrowthNE, var.equal = FALSE)
res


#part 2

#libraries that are required
library("seqinr")
library("rBLAST")
library("R.utils")

#Download the E.coli CDS sequence from the Ensembl FTP page
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "ecoli.fa.gz")

# uncompress the file
gunzip("ecoli.fa.gz")
# create the blast DB
makeblastdb("ecoli.fa",dbtype="nucl", "-parse_seqids")

#Download the sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = "sample.fa")

#Read the sample file into R
d <- read.fasta("sample.fa")
mygene <- d[[3]]
mygene

#Length of sequence in bp
str(mygene)
length(mygene)

#Proportion of GC content
seqinr::GC(mygene)

#function to create blast
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast.R")
source("mutblast.R")

#Blast search for E. coli genes that matches best
res <- myblastn_tab(myseq = mygene, db = "ecoli.fa")

#Blast results
res
View(res)

str(res)
head(res)

#making mutations to mygene
mygene_mutation <- mutator(mygene,20)
res <- myblastn_tab(myseq = mygene_mutation, db = "ecoli.fa")
res

#first need to write a blast index
write.fasta(mygene,names="mygene",file.out = "mygene.fa")
makeblastdb(file = "mygene.fa",dbtype = "nucl")

# test with 100 mismatches
mygene_mutation <- mutator(myseq=mygene,100)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#test with 150 mismatches
mygene_mutation <- mutator(myseq=mygene,150)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#test with 200 mismatches
mygene_mutation <- mutator(myseq=mygene,200)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#test with 250 mismatches
mygene_mutation <- mutator(myseq=mygene,250)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#test with 300 mismatches
mygene_mutation <- mutator(myseq=mygene,450)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

plot(res$bitscore)

