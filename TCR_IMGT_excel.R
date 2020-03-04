args<-commandArgs(TRUE)
infile=args[1]
outfile=args[2]
infile2 = args[3]

library(readxl)
library(plyr)
library(stringr)

#load in the Excel file from IMGT
imgt <- read_excel(infile,sheet = 1)

if (lengths(strsplit(as.character(imgt$`Sequence ID` ),"-")[1]) == 1) {
  print("Simple ID input")
  achain <- imgt[as.numeric(str_extract(imgt$`Sequence ID`, "[0-9]+")) %% 2 != 0 ,]
  bchain <- imgt[as.numeric(str_extract(imgt$`Sequence ID`, "[0-9]+")) %% 2 == 0 ,]
} else {
  print("Other ID input")
  achain <- imgt[grep("JMWA",imgt$`Sequence ID`),]
  bchain <- imgt[grep("JMWB",imgt$`Sequence ID`),]
}


if (length(args)==2) {
  print("No Second File Passed")
} else if (length(args)==3) {
  imgt2 <- read_excel(infile2,sheet = 1)
  imgt <- rbind(imgt,imgt2)
} else { 
  print("Wrong number of arguments, check command")
  stop()
}

#get only these columns from the excel file
interestA <- c("Sequence ID","V-GENE and allele","J-GENE and allele","AA JUNCTION","V-DOMAIN Functionality")
interestB <- c("Sequence ID","V-GENE and allele","J-GENE and allele","D-GENE and allele","AA JUNCTION","V-DOMAIN Functionality")

#make a dataframe from columns of interest
dfa<- data.frame(achain[,interestA])
dfb<- data.frame(bchain[,interestB])

#Get the well number
if (lengths(strsplit(as.character(imgt$`Sequence ID` ),"-")[1]) == 1) {
  dfa$number <- dfa[,1]
  dfb$number <- dfb[,1]
} else {
  dfa$number <- sapply(strsplit(dfa[,1],"-"), "[", 3)
  dfb$number <- sapply(strsplit(dfb[,1],"-"), "[", 3)
}


#order by well number for pairing
dfa <- dfa[order(as.numeric(gsub("([0-9]+)([A-Z]+)", "\\1", dfa$number)),gsub("([0-9]+)([A-Z]+)", "\\2", dfa$number)),]
dfb <- dfb[order(as.numeric(gsub("([0-9]+)([A-Z]+)", "\\1", dfb$number)),gsub("([0-9]+)([A-Z]+)", "\\2", dfb$number)),]


#give order for later sorting 
dfa$order <- 1:nrow(dfa)
dfb$order <- 1:nrow(dfb)

#retrieve only the chain form
clean <- function(region){
  t1 <- gsub("\\s*\\([^\\)]+\\)","",as.character(region))
  t2 <- str_replace(t1, "F", "")
  t3 <- gsub("Homsap", "",t2)
  t4 <- gsub(" |, ", "",t3)
  return(t4)
}

#Get only the amino acids(remove any notes/comments)
cleanAA <- function(region){
  t1 <- gsub("\\s*\\([^\\)]+\\)","",as.character(region))
  t6 <- gsub(" |, ", "",t1)
  return(t6)
}


#get coordinates for matching
dfa$genes <- paste(clean(dfa$V.GENE.and.allele),clean(dfa$J.GENE.and.allele),cleanAA(dfa$AA.JUNCTION), sep = ":")
dfb$genes <- paste(clean(dfb$V.GENE.and.allele),clean(dfb$J.GENE.and.allele),clean(dfb$D.GENE.and.allele),cleanAA(dfb$AA.JUNCTION), sep = ":")

#get the frequency of each of the genes 
freqa <- count(dfa$genes)
freqb <- count(dfb$genes)

#number for each unique chain
freqa$group <- 1:nrow(freqa)
freqb$group <- 1:nrow(freqb)

#now each repeated chain will be grouped 
dfa <- merge(dfa,freqa,by.x='genes', by.y=1)
dfb <- merge(dfb,freqb,by.x='genes', by.y=1)

#makes list of each chain that is repeated
groupings <- function(df){
  grouped <- list()
  number <- df$number 
  names(number) <- df$group
  for(i in 1:length(unique(df$group))){
    grouped[[i]]<- paste(number[names(number) == i], collapse=",")
  }
  df <- data.frame(row.names  = 1:length(unique(df$group)) ,unlist(grouped))
  return(df)
}

#makes list of each pair of a and b that are repeated 
pair_groupings <- function(df){
  grouped <- list()
  numbera <- df$number.x
  names(numbera) <- df$group
  numberb <- df$number.y 
  names(numberb) <- df$group
  for(i in 1:length(unique(df$group))){
    grouped[[i]]<- paste(paste(numbera[names(numbera) == i], collapse=","),paste(numberb[names(numberb) == i], collapse=","),collapse = ",")
  }
  df <- data.frame(row.names  = 1:length(unique(df$group)) ,unlist(grouped))
  return(df)
}

#Outputs unproducive chains
productive <- function(df){
  grouped <- list()
  number <- df$number 
  group <- df$group
  prod <- df$V.DOMAIN.Functionality
  df <- data.frame(number,group,prod)
  for(i in 1:length(unique(df$group))){
    grouped[[i]]<- paste(df[df$group == i & df$prod != "productive",]$number, collapse=",")
  }
  df <- data.frame(row.names  = 1:length(unique(df$group)) ,unlist(grouped))
  return(df)
}

#add to table for each chain 
groupinga <- groupings(dfa)
groupingb <- groupings(dfb)
proda <- productive(dfa)
prodb <- productive(dfb)
colnames(proda) <- "Unproductive_A_Chains"
colnames(prodb) <- "Unproductive_B_Chains"


dfa <- merge(dfa,groupinga,by.x='group', by.y=0)
dfb <- merge(dfb,groupingb,by.x='group', by.y=0)
dfa <- merge(dfa,proda,by.x='group', by.y=0)
dfb <- merge(dfb,prodb,by.x='group', by.y=0)

#merge a and b chain tables in order for pairing 
merged <- merge(dfa,dfb, by ="order",all = T)

#make coordinate for every pair
merged$pair <- paste(merged$genes.x,merged$genes.y, sep = ":")

#count duplicate pairs
freqab <- count(merged$pair)

#number for each unique pair
freqab$group <- 1:nrow(freqab)

#each pair is now grouped
merged <- merge(merged,freqab,by.x='pair', by.y=1)

#Print a and b chains that are repeated
groupingab <- pair_groupings(merged)

#add to merge dataframe 
merged <- merge(merged,groupingab,by.x='group', by.y=0)

#add unproductive chain pairs 
merged$Unproductive_Chains <- paste(merged$Unproductive_A_Chains,merged$Unproductive_B_Chains)

#tidy up dataframe
final <- merged[!duplicated(merged$group),-(grep("genes|order|pair|number|group.x|group.y|V.DOMAIN.Functionality|Unproductive_B_Chains|Unproductive_A_Chains",colnames(merged)))]
final <- final[order(final$freq,decreasing = T),-1]
colnames(final) <- c("Sequence_ID_Alpha", "Alpha_V","Alpha_J","Alpha_CDR3","Frequency_of_A","Same_A_Wells","Sequence_ID_Beta", "Beta_V","Beta_J","Beta_D","Beta_CDR3","Frequency_of_B","Same_B_Wells","Frequency_of_Pair","Same_Pair_Wells", "Unproductive Wells")

print(paste("There were",length(rownames(achain)),"Alpha chains in the input file"))
print(paste("There were",length(rownames(bchain)),"Beta chains in the input file"))
print(paste("There were",length(rownames(final)),"Unique Pairs"))

nums <- data.frame(c(paste("There were",length(rownames(achain)),"Alpha chains in the input file"),
                     paste("There were",length(rownames(bchain)),"Beta chains in the input file"),
                     paste("There were",length(rownames(final)),"Unique Pairs")))

#write output
write.table(final,file = outfile,row.names = F, quote = F, sep ="\t")
write.table(nums,file = outfile,row.names = F, quote = F, sep ="\t",append = T,col.names = F)
