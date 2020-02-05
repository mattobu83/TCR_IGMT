#specify full file path of IGMT output in first argument
#specify full file name for output og script in second argument
#specify the number of sequences in the fastq file in the third argument

args<-commandArgs(TRUE)
infile=args[1]
outfile=args[2]
nseq = args[3]

#load required packages
library(data.table)

#read in output from IGMT
test <- scan(infile, what = "list", sep = "\n")

#create empty list
info <- list()

#seperate information for each sequence
for (n in 1:nseq){
  index1 <- grep("Sequence number", test)[n]
  index2 <- ifelse(is.na(grep("Sequence number", test)[n+1]),length(test), grep("Sequence number", test)[n+1])
  info[[n]] <- test[index1:(index2-1)]
}

#seperate a chain and bchain info
achain <- info[grep("JMWA",info)]
bchain <- info[grep("JMWB",info)]

#outputs percentage for v/d/j
matchper <- function(chain){
  identity <- ifelse(is.na(sapply(strsplit(chain,";"),"[",4)[1]),"NA",sapply(strsplit(chain,";"),"[",4))
  percent <- ifelse(is.na(sapply(strsplit(chain,";"),"[",4)[1]),"NA",sapply(strsplit(identity,"="),"[",2))
  return(percent)
}

#Gets info for A chain, places into dataframe 
checkA <- function(infoA){
  name <- sapply(strsplit(infoA[1],":"), "[", 2)
  number <- sapply(strsplit(infoA[1],"-"), "[", 3)
  v <- ifelse(is.na(grep("V-GENE and allele;Homsap TRA",infoA)[1]),"NA",infoA[grep("V-GENE and allele;Homsap TRA",infoA)])
  j <- ifelse(is.na(grep("J-GENE and allele;Homsap TRA",infoA)[1]),"NA",infoA[grep("J-GENE and allele;Homsap TRA",infoA)])
  cdr <- ifelse(is.na(grep("CDR-IMGT lengths",infoA)[1]),"NA",infoA[grep("CDR-IMGT lengths",infoA)])
  res <- ifelse(is.na(grep("IMGT/V-QUEST results:",infoA)[1]),"NA",infoA[(grep("IMGT/V-QUEST results:",infoA))+1])
  prod <- ifelse(is.na(sapply(strsplit(res," "),"[",2)[1]),"NA",sapply(strsplit(res," "),"[",2))
  prod2 <- ifelse(prod == "Productive", "Yes", "No")
  vchain <- ifelse(is.na(sapply(strsplit(v,";"),"[",2)[1]),"NA",sapply(strsplit(v,";"),"[",2))
  jchain <- ifelse(is.na(sapply(strsplit(j,";"),"[",2)[1]),"NA",sapply(strsplit(j,";"),"[",2))
  cdr2 <- ifelse(is.na(sapply(strsplit(cdr,";"),"[",4)[1]),"NA",sapply(strsplit(cdr,";"),"[",4))
  percentv <- matchper(v)
  percentj <- matchper(j)
  df <- data.frame(nameA = name,numberA = number,vchain,percentv,jchain,percentj,cdr2,Productive = prod2)
  return(df)
}


#Gets info for B chain, places into dataframe 
checkB <- function(infoB){
  name <- sapply(strsplit(infoB[1],":"), "[", 2)
  number <- sapply(strsplit(infoB[1],"-"), "[", 3)
  v <- ifelse(is.na(grep("V-GENE and allele;Homsap TRB",infoB)[1]),"NA",infoB[grep("V-GENE and allele;Homsap TRB",infoB)])
  j <- ifelse(is.na(grep("J-GENE and allele;Homsap TRB",infoB)[1]),"NA",infoB[grep("J-GENE and allele;Homsap TRB",infoB)])
  d <- infoB[ifelse(is.na(grep("D-GENE and allele",infoB)[1]),"NA",grep("D-GENE and allele",infoB))]
  cdr <- ifelse(is.na(grep("CDR-IMGT lengths",infoB)[1]),"NA",infoB[grep("CDR-IMGT lengths",infoB)])
  res <- ifelse(is.na(grep("IMGT/V-QUEST results:",infoB)[1]),"NA",infoB[(grep("IMGT/V-QUEST results:",infoB))+1])
  prod <- ifelse(is.na(sapply(strsplit(res," "),"[",2)[1]),"NA",sapply(strsplit(res," "),"[",2))
  prod2 <- ifelse(prod == "Productive", "Yes", "No")
  vchain <- ifelse(is.na(sapply(strsplit(v,";"),"[",2)[1]),"NA",sapply(strsplit(v,";"),"[",2))
  jchain <- ifelse(is.na(sapply(strsplit(j,";"),"[",2)[1]),"NA",sapply(strsplit(j,";"),"[",2))
  dchain <- ifelse(is.na(sapply(strsplit(d,";"),"[",2)[1]),"NA",sapply(strsplit(d,";"),"[",2))
  cdr2 <- ifelse(is.na(sapply(strsplit(cdr,";"),"[",4)[1]),"NA",sapply(strsplit(cdr,";"),"[",4))
  percentv <- matchper(v)
  percentj <- matchper(j)
  df <- data.frame(nameB = name,numberB = number,vchain,percentv,jchain,percentj,dchain,cdr2, Productive = prod2)
  return(df)
}

#applies functions above on all sequences in list 
test1 <- lapply(achain, checkA)

test2 <- lapply(bchain, checkB)

#make dataframe from list of dataframes, orders numerically, then by letter
dfa <- rbindlist(test1)
dfa <- dfa[order(as.numeric(gsub("([0-9]+)([A-Z]+)", "\\1", dfa$number)),gsub("([0-9]+)([A-Z]+)", "\\2", dfa$number))]

dfb <- rbindlist(test2)
dfb <- dfb[order(as.numeric(gsub("([0-9]+)([A-Z]+)", "\\1", dfb$number)),gsub("([0-9]+)([A-Z]+)", "\\2", dfb$number))]


#match A and B 
final <-cbind(dfa,dfb)

#write output
write.table(final, file = outfile, sep = "\t", row.names = F, quote = F)
