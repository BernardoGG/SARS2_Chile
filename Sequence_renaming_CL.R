################################################################################
########### SARS-CoV-2 Chile sequence re-namer #################################
################################################################################

############## Rhys Inward, Bernardo Gutierrez #################################

# Load data
data1 <- as.data.frame(fread("SARS2_Chile/Data/SARS2_CL_Alpha_ISP.csv"))
data2 <- as.data.frame(fread("SARS2_Chile/Data/SARS2_CL_Alpha_GISAID.csv"))

# Generate .txt file for sequence naming
df1 <- data%>% unite (hk,'Accession ID': 'Virus name',sep=",")
df2 <- data %>% unite (hk2,'Accession ID': 'Collection date',sep=",")
df3 <- df2 %>% unite (hk3,'hk2': 'Location',sep=",")
df4 <- df3 %>% unite (hk4,'hk3':'Pango lineage',sep=",")
df5 <- df4 %>% unite (hk5,'hk4':'surveillance',sep=",")
df6 <- bind_cols (df1$hk,df5$hk5)
df6 <- df6 %>% unite (hk,'...1':'...2' ,sep="|")
colnames(df6) <- "EpiID,SeqName,Date,Location,lineage,surveillance"
df7 <- df6[,1, drop=FALSE]

# To export as text
write.table(df7, file = "Data/21jun2022_EPI_ISL_summary.txt", row.names = FALSE,
            sep = "\t", quote = FALSE)

#Use the csv to copy and paste Accession ID into GISAID to get required sequences

#make sure to point to the correct directory
path <- "~/rhys/Chile/Data/"
fnames<- dir(path)
fnames<- fnames[grep(".fasta",fnames)]

#remember to rename the file with all the sequences BetaCoV_Wuhan_177seqs_22dec2020 
# create one file with all of the sequences
# need to rename file
combinedFname <- paste("gisaid_hcov-19_2022_03_13_16.fasta", sep = "")
pos <- grep(combinedFname,dir(),fixed=TRUE)

#read the sequence information file
info <- read.csv(paste(path,"08mar2022_EPI_ISL_summary.txt",sep=""))
print(info)
#read file and alter the sequence names 
rootname <-	paste(gsub("\\.fas","",combinedFname),"_beastNames",sep="")
pos <- grep(paste(rootname,".fas",sep=""),dir(),fixed=TRUE)
if (length(pos)==0) {
  
  seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
  taxa <- as.matrix(attributes(seqs)$names)
  isName <- apply(taxa, 1, getEl, ind=3, sep="/")
  epiISL <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  minds  <- match(epiISL, info$EpiID)
  all( epiISL==info$EpiID[minds] )
  dateTxt <- as.matrix(info$Date[minds])
  lineage <- as.matrix(info$lineage[minds])
  surveillance <- as.matrix(info$surveillance[minds])
  decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  location<- as.matrix(info$Location[minds])
  country <- apply(as.matrix(location), 1, getEl, ind=1, sep=" / ")
  state   <- apply(as.matrix(location), 1, getEl, ind=2, sep=" / ")
  place   <- apply(as.matrix(location), 1, getEl, ind=3, sep=" / ")
  
  pos     <- which((state=="Kanagawa Prefecture") & is.na(place))
  place[pos] <- "Yokohama"
  print(paste("Changed place of",info$SeqName[pos],"from NA to Yokohama (capital city)"))
  
  newTaxa <- paste(state,place,isName,epiISL,lineage,surveillance,dateTxt,decDate,sep="|")
  newTaxa <- gsub(" ","_",newTaxa)
  attributes(seqs)$names <- newTaxa
  write.dna(seqs, file=paste(path,rootname,".fas",sep=""), format="fasta", nbcol=-1, colsep="")
  
  newInfo <- cbind(newTaxa,epiISL,isName,dateTxt,country,state,place,decDate)
  colnames(newInfo) <- c("SeqName","EPI_ISL","IsolateName","CollectionDate","Country","State","Place","decDate")
  write.table(newInfo,file=paste(path,rootname,"_infoTbl.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
} else {
  print("Already done renaming")
}

#load results of QC

nextclade <- as.data.frame(fread("Data/nextclade.tsv"))

#select for sequences that weren't deemed good 

aligned_seq <- read.fasta("Data/sequences.aln.fasta")

taxa <- as.data.frame(as.matrix(attributes(aligned_seq)$names))

taxa <- as.data.frame(taxa[!duplicated(taxa$V1), ])

nextclade_remove <- filter(nextclade, qc.overallStatus != 'good')

species.to.remove <- nextclade_remove$seqName

vec.names<-unlist(lapply(strsplit(names(aligned_seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which(! vec.names %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=aligned_seq[vec.tokeep], names=names(aligned_seq)[vec.tokeep], file.out="Data/chile_data.fas")

#plot of data

#split location of data 

treemmer_seq <- read.fasta("Data/chile_data.fas")
treemmer_seq_name <- as.data.frame(as.matrix(attributes(treemmer_seq)$names))
taxa_split_spatial_temporal <- data.frame(do.call('rbind',strsplit(as.character(treemmer_seq_name$V1),'|',fixed = TRUE)))

table(taxa_split_spatial_temporal$X6)
#plot sequences over time and compare 
#want to sum up up all those on the same day from the same country 

by_location_spatial_temporal_sum <- taxa_split_spatial_temporal %>%
  group_by(X1, X7) %>%
  dplyr :: summarize(count = n())
colnames(by_location_spatial_temporal_sum) <- c("Location", "Date", "count")

by_location_spatial_temporal_sum$Date <- as.Date (by_location_spatial_temporal_sum$Date )

by_location_spatial_temporal_sum <- by_location_spatial_temporal_sum %>% group_by(week = cut(Date, "week", start.on.monday = FALSE),Location) %>% 
  dplyr :: summarise(count = sum(count))

by_location_spatial_temporal_sum$week <- as.Date(by_location_spatial_temporal_sum$week)

by_location_spatial_temporal_sum <- by_location_spatial_temporal_sum %>%
  filter(count >= 1)

ggplot(by_location_spatial_temporal_sum, aes(x= week, y = count, fill = Location
)) + geom_bar(stat = 'identity', position = 'fill') +
  labs(x="Date of Collection",y="Proportion of Sequences from Each Country") + theme_bw() +
  scale_x_date(date_breaks = "2 month", date_labels = "%Y %b %d")

ggplot(by_location_spatial_temporal_sum, aes(x= week, y = count, fill = Location
)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Date of Collection",y="Total Number of Sequences from Each Country") + theme_bw() +
  scale_x_date(date_breaks = "2 month", date_labels = "%Y %b %d")
