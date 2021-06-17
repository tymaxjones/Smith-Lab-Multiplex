install.packages("devtools")
library(devtools)
find_rtools() 

install.packages("BiocManager")
BiocManager::install()

BiocManager::install("GenomeInfoDb")
library(GenomeInfoDb)

BiocManager::install("openPrimeRui")

library(openPrimeRui)

openPrimeRui::startApp()

library(openPrimeR)

fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                          "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
# Load the template sequences from 'fasta.file'
seq.df.simple <- read_templates(fasta.file)


BiocManager::install("Biostrings")
BiocManager::install("TmCalculator")
library(TmCalculator)

library(Biostrings)
dna <- readDNAStringSet("C:/Users/tymax/Desktop/Sars-CoV-2.fasta")

#64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)

melt_temp <- 64.9 + (41 * (4 - 16.4)/(18))



#primer should be 18 30
dna2 <- dna$`NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome`

current_start_position <- 1
current_end_position <- 18


current_seq <- dna2[seq(current_start_position,current_end_position)]

matchPattern(current_seq, dna2)

sum(letterFrequency(current_seq, c("C","G")))

current_seq


library(TmCalculator)



TmCalculator::Tm_GC(as.character(current_seq))
Tm_NN(as.character(current_seq))
Tm_Wallace(as.character(current_seq))


Tm_GC("TTGGACTAATTGACTCTGACTGCCG" )
Tm_NN("TTGGACTAATTGACTCTGACTGCCG" )
Tm_Wallace("TTGGACTAATTGACTCTGACTGCCG" )


Tm_GC("TTGGACTAATTGACTCTGACTGCCG", Na = 50, Mg = 0 )

#This is for analysis of the different melting temperatures stuff

random <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(c("A", "C", "T", "G")), n, TRUE), FALSE))
  paste0(a, sprintf(sample(c("A", "C", "T", "G"), n, TRUE))
}


seqs <- rep("", 100)
id <- seq(1,100)
tm1 <- rep(0, 100)
tm2 <- rep(0, 100)
tm3 <- rep(0, 100)
for (i in 1:1000) {
  seq <- paste0(sample(c("A", "C", "T", "G"), replace = TRUE, size = sample(18:30, size = 1)), collapse = "")
  seqs[i] <- seq
  tm1[i] <- Tm_GC(seq)
  tm2[i] <- Tm_NN(seq)
  tm3[i] <- Tm_Wallace(seq)
}

library(ggplot2)


df <- data.frame(id,tm1,tm2,tm3)

install.packages("tidyverse")
library(tidyverse)

ggplot(df) +
  geom_line(aes(x = id, y = tm1), color = "red") +
  geom_line(aes(x = id, y = tm2), color = "blue") +
  geom_line(aes(x = id, y = tm3), color = "orange") +
  geom_smooth(aes(x = id, y = tm1), color = "red") +
  geom_smooth(aes(x = id, y = tm2), color = "blue") +
  geom_smooth(aes(x = id, y = tm3), color = "orange")



df$dif12 <- df$tm2 - df$tm1  
df$dif13 <- df$tm3 - df$tm1  
df$dif23 <- df$tm3 - df$tm2  


ggplot(df) +
  geom_line(aes(tm1, dif12))
ggplot(df) +
  geom_line(aes(tm1, dif13))

ggplot(df) +
  geom_line(aes(tm2, dif23))





########################################################################################################


########################################################################################################

library(TmCalculator)
library(Biostrings)

#Get the sequence
dna <- readDNAStringSet("C:/Users/tymax/Desktop/Sars-CoV-2.fasta")
dna2 <- dna$`NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome`

#Sequence Stuff
current_start_position <- 1
current_end_position <- 18

#This subsets our sequence to just the portion that we want
current_seq <- dna2[seq(current_start_position,current_end_position)]

matchPattern(current_seq, dna2) #Get the specificity
sum(letterFrequency(current_seq, c("C","G"))) #the number of G and C



#remove the poly a tail

index <- length(dna2)
while (as.character(dna2[index]) == "A"){
  index = index - 1
}
dna3 <- dna2[1:index]


#Create some vectors that we will store info in
start_pos <- rep(1:length(dna3), each = length(18:30))
dif <- rep(17:29, times = length(unique(start_pos)))
end_pos <- start_pos + dif
length <- dif + 1


#turn them into a df
primer_df <- data.frame(start_pos,end_pos,length)
primer_df <- primer_df[primer_df$end_pos < length(dna3),]

matchPattern(dna3[1:18], dna3)
#now we are going to get the specificity

n <- nrow(primer_df)
specificity <- rep(0,n)
for (i in 1:10) {
#  specificity[i] = length(matchPattern(dna3[seq(primer_df$start_pos[i], primer_df$end_pos[i])], dna3))
  
  print(i)
  
}
i = 1
while (i <= n) {
  start <-  primer_df$start_pos[i]
  stop <- primer_df$end_pos[i]
  cur_seq <- dna3[seq(start,stop)]
  
  if (length(matchPattern(cur_seq, dna3)) == 1 ) {
    #This gives us the index of all the values for our current start pos
    all_start <- which(primer_df$start_pos == start)
    #make all those indices 1
    specificity[all_start] <- 1
    
  }else{
    #Get all of the indices of our current start pos
    all_start <- which(primer_df$start_pos == start)
    
    #loop through all of the values in our current start pos
    for (j in all_start) {
      specificity[j] = length(matchPattern(dna3[seq(primer_df$start_pos[j], primer_df$end_pos[j])], dna3))
      
    }
  }
  #increase i to go past all the indices of our current start pos
  i <- i + length(all_start)
  
  
  
  if (i %% 50 == 0) {
    print(i)
  }
  
  
}

primer_df$specificity <- specificity

View(primer_df)


sum(letterFrequency(dna3[1:18], c("G", "C")))


gc_content <- rep(0, n)
for (i in 1:n) {
  gc_content[i] <- sum(letterFrequency(dna3[ seq(primer_df$start_pos[i], primer_df$end_pos[i])], c("G", "C")))
  
  if (i %% 1000 == 0 ) {
    print(i)
    
  }
  
}

letterFrequency(dna3[ seq(primer_df$start_pos[i], primer_df$end_pos[i])], c("G", "C", "A", "T"))

primer_df$gc_content <- gc_content
library(readr)
write_csv(primer_df, "C:/Users/tymax/Desktop/Sars-CoV-2_tbl.csv")


dna

#Do primers that are 18-22/23
#Avoid long stretch of A,C,T,G no longer than a stretch of 4
#avoid a 3 prime T
#GC content between 40 and 60 percent
#online primer design tips
#melting temperature closer to 60
#specificity:
#inclusive for innfluenza A substrains
#exclusive ot other types of unfluenza

##########################################################################
library(RJSONIO)
report <- RJSONIO::fromJSON("C:/Users/tymax/Downloads/scheme/scheme.report.json")

View(report)

primers <- read_tsv("C:/Users/tymax/Downloads/scheme/scheme.primer.tsv")

View(primers)
dna3
