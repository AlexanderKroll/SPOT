# ==============================================================================
# Need a database with the following entries:
#
# Protein
# - Name - name of protein
# - TCNumber - classification of a protein
# - UniProt - UniProt id for a protein
# - Sequence - sequence of the protein
#
# Substrate
# - Substrate - name of substrate
# - CHEBI - identification of a substrate
# - InChI - inchi code of substrate
# ==============================================================================

# ==============================================================================
# 1. Read data from "tcdb.fasta" and "getSubstrates.tsv" from
# "https://www.tcdb.org/" and process them for further use
# ==============================================================================

library(Biostrings) # To read .fasta files
tcdb <- readAAStringSet("original_data/tcdb.fasta")
map <- read.csv("original_data/getSubstrates.tsv",
                header=FALSE, 
                sep="\t", 
                stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1.1 Process tcdb.fasta
# ------------------------------------------------------------------------------
# Get AA sequences
tcdb_seq <- as.character(tcdb, use.names = FALSE) 

# Split rest of information into 4 parts 
# ("gnl", "TC-DB", UniProt, TCNumber + Name) - last two needed
x <- strsplit(names(tcdb), split="\\|") 

# Some problems when splitting; 3 wrong splits
x_n <- unlist(lapply(x, length))
levels(as.factor(x_n))
sum(x_n == 5)

# Get entries with length 5 and correct them manually
bl <- x[x_n == 5]; bl

v1 <- bl[1][[1]]
v1 <- c(v1[1:3], paste(v1[4], "|", v1[5], sep=""))
bl[1][[1]] <- v1

v2 <- bl[2][[1]]
v2 <- c(v2[1:3], paste(v2[4], "|", v2[5], sep=""))
bl[2][[1]] <- v2

v3 <- bl[3][[1]]
v3 <- c(v3[1:3], paste(v3[4], "|", v3[5], sep=""))
bl[3][[1]] <- v3

x[x_n == 5] <- bl; x[x_n == 5]

# Get only relevant entries and split into names and TCNumber
all_cols <- matrix(unlist(x), ncol=4, byrow=TRUE)[,3:4]
uniprot <- all_cols[,1]
tcnumber_names <- all_cols[,2]

f <- function(x) {
  c(x[1], paste(x[-1], collapse=""))
}

tcnumber_names <- matrix(unlist(lapply(strsplit(tcnumber_names, split="\\ "),
                                       f)),
                         ncol=2, byrow=TRUE)

# Create tcdb dataframe and save it in a .csv file for further use
tcdb <- data.frame(Name = tcnumber_names[,2],
                   TCNumber = tcnumber_names[,1],
                   UniProt = uniprot,
                   Sequence = tcdb_seq)
write.csv(tcdb, "component_data/tcdb.csv", row.names = F)

# ------------------------------------------------------------------------------
# 1.2 Process getSubstrates.tsv
# ------------------------------------------------------------------------------
# Split data in second row - more than one entries for one TCNumber allowed
sepdata1 <- strsplit(map[,2], split="\\;|\\|")
m <- matrix(unlist(sepdata1), ncol=2, byrow=TRUE)

# Create map dataframe and save it in a .csv file for further use
map <- data.frame(TCNumber = rep(map[,1], lengths(sepdata1)/2),
                  ChEBI = m[,1],
                  Substrate = m[,2])
write.csv(map, "component_data/map_tcnumber_to_chebi.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 1.3 Merge tcdb and map dataframe for further use - no InChI
# ------------------------------------------------------------------------------
data <- merge(tcdb, map, by="TCNumber")
write.csv(data, "component_data/database_without_inchi.csv", row.names = FALSE)

# ==============================================================================
# 2. Get InChI from CHEBI
# ==============================================================================
rm(list = ls())
library(webchem)
map <- read.csv("component_data/map_tcnumber_to_chebi.csv")

# ------------------------------------------------------------------------------
# 2.1 Get InChI codes from https://www.ebi.ac.uk/
# ------------------------------------------------------------------------------
chebiids <- levels(as.factor(map$ChEBI))
chebi_list <- chebi_comp_entity(chebiids)

prop <- lapply(chebi_list, "[[", "properties")
inchi <- sapply(prop, "[[", "inchi")

conv <- data.frame(InChI = inchi)
conv$ChEBI <- rownames(conv)
rownames(conv) <- NULL
conv <- conv[c("ChEBI", "InChI")]

merged_map <- merge(map, conv, by="ChEBI")

# ------------------------------------------------------------------------------
# 2.2 Get InChI codes from "https://pubchem.ncbi.nlm.nih.gov/"
# ------------------------------------------------------------------------------
missing_map <- merged_map[is.na(merged_map$InChI),]
unique_names <- missing_map[!duplicated(missing_map[,"ChEBI"]),]$Substrate

# Get Pubchem IDs - NA if not found - and split them
cids <- get_cid(unique_names)

cids_without_na <- na.omit(cids)
cids_with_na <- cids[is.na(cids$cid),]
names_na <- unique(cids_with_na$query)

# Get information about rows frequency (more than one value equal?)
x <- data.frame(table(cids_without_na$query))

# Get duplicates and their names
dups <- cids_without_na[cids_without_na$query %in%
                          x$Var1[x$Freq > 1],]
names_dups <- unique(dups$query)
View(dups)
# Remove duplicates
ucids <- cids_without_na[cids_without_na$query %in%
                           x$Var1[x$Freq == 1],]

# Add InChI-Codes
ucids$inchi <- pc_sect(ucids$cid, "inchi")$Result
dups$inchi <- pc_sect(dups$cid, "inchi")$Result

# Merge 
x <- ucids[,c("query", "inchi")]
colnames(x) <- c("Substrate", "InChI")

View(x)
View(merged_map)

merged_map2 <- merge(merged_map, x, by="Substrate", all=T)
View(merged_map2)
merged_map2$InChI <- merged_map2$InChI.x
merged_map2$InChI[!is.na(merged_map2$InChI.y)] <- 
  merged_map2$InChI.y[!is.na(merged_map2$InChI.y)]
merged_map2$InChI.x <- NULL
merged_map2$InChI.y <- NULL



# ------------------------------------------------------------------------------
# 2.3 Save data
# ------------------------------------------------------------------------------
write.csv(na.omit(merged_map2), "component_data/chebi_to_inchi.csv", row.names = F)
write.csv(merged_map2, "component_data/chebi_to_inchi_with_na.csv", row.names = F)
write.csv(merged_map2[is.na(merged_map2$InChI),], "component_data/missing_inchi.csv", row.names = F)
write.csv(dups, "component_data/no_unique_inchi.csv", row.names = F)

# ==============================================================================
# 3. Merge data and analyze it
# ==============================================================================
rm(list = ls())
data <- read.csv("component_data/database_without_inchi.csv")
inchi <- read.csv("component_data/chebi_to_inchi.csv")
View(inchi)
View(data)

# ------------------------------------------------------------------------------
# 3.1 Merge data with InChI mapping by CHEBI
# ------------------------------------------------------------------------------
inchi$Substrate <- NULL
inchi$TCNumber <- NULL

data2 <- unique(merge(inchi, data, by = "ChEBI"))
View(data2)


write.csv(data2, "database.csv", row.names = F)
test <- read.csv("database.csv")
all.equal(data2, test)
View(test)
