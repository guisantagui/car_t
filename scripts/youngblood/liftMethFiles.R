################################################################################
# CAR-T: lift coordinates of methylation bed files that have been built with   #
# one reference genome to a different reference genome.                        #
################################################################################

if(!require("argparser", quietly = T)){
        install.packages("argparser",
                         repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("BiocManager", quietly = T)){
        install.packages("BiocManager",
                         repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("rtracklayer", quietly = T)) BiocManager::install("rtracklayer",
                                                              update = F)
if(!require("GenomicRanges", quietly = T)) BiocManager::install("GenomicRanges",
                                                                update = F)

library(rtracklayer)
library(GenomicRanges)
library(argparser)

# Parser
################################################################################
parser <- arg_parser("lift coordinates of bed files that have been built with one reference #
# genome to a different reference genome.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--chain",
                               "--reorder",
                               "--compress",
                               "--outTag",
                               "--outDir"),
                       help = c("Input bed file.",
                                "Chain file with mappings from one start reference to target reference (Obtained from http://hgdownload.cse.ucsc.edu/downloads.html#human).",
                                "If chromosomes and positions should be reordered in output",
                                "If output file should be compressed.",
                                "Tag to add to new files.",
                                "Output directory."),
                       flag = c(F, F, T, T, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# read.csv but faster
readCsvFst <- function(pth, header = T, row.names = T){
        df <- data.frame(data.table::fread(pth, header = header))
        if(row.names){
                rownames(df) <- df$V1
                df <- df[, colnames(df) != "V1"]
        }
        return(df)
}

# write.csv, but faster
writeCsvFst <- function(df, file, rowNames = T, colNames = T, sep = ","){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames, sep = sep)
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Extract numeric part from chromosome names
extract_number <- function(x){
        num_part <- gsub("chr", "", x)
        if(is.na(suppressWarnings(as.numeric(num_part)))){
                return(Inf) # Assign infinity to non-numeric parts so they appear last
        }else{
                return(as.numeric(num_part))
        }
}

# Directory stuff
################################################################################
inFile <- parsed$input
chainFile <- parsed$chain
reord <- parsed$reorder
compress <- parsed$compress
outTag <- parsed$outTag
outDir <- parsed$outDir

inFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/data/youngblood/GSM5677818_2-1245637.txt.gz"
chainFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/data/hg19ToHg38.over.chain.gz"
reord <- T
compress <- T
outTag <- "hg38"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/data/youngblood/hg38Conv"

outDir <- addSlashIfNot(outDir)
createIfNot(outDir)


outName <- basename(inFile)
outName <- gsub(".gz", "", outName)
outName <- gsub(".txt", "", outName)

outName <- sprintf("%s%s_%s.txt", outDir, outName, outTag)

# Load data
################################################################################
df <- readCsvFst(inFile, header = F, row.names = F)

gr <- GRanges(seqnames = Rle(df$V1),
              ranges = IRanges(start = df$V2, end = df$V2),
              coverage = df$V3,
              numCs = df$V4,
              numTs = df$V3 - df$V4)

if(grepl(".gz", chainFile)){
        system(sprintf("gzip -d -c %s > %s", chainFile,
                       gsub(".gz", "", chainFile, fixed = T)))
        chainFile <- gsub(".gz", "", chainFile, fixed = T)
}

chain <- import.chain(chainFile, exclude = "_")

# Lift the coordinates
################################################################################
seqlevelsStyle(gr) = "UCSC"
gr_lifted <- liftOver(gr,  chain)

gr_lifted <- unlist(gr_lifted)
genome(gr_lifted) <- outTag

# Parse result to output
################################################################################

# Convert to dataframe
gr_lifted <- as.data.frame(gr_lifted)
gr_lifted$seqnames <- as.character(gr_lifted$seqnames)

if(reord){
        # Order rows according to chromosomes
        print("Reordering chromosomes and positions...")
        ordChroms <- sapply(unique(gr_lifted$seqnames), extract_number)
        ordChroms_nonNum <- sort(names(ordChroms)[is.infinite(ordChroms)])
        ordChroms <- names(sort(ordChroms[!is.infinite(ordChroms)]))
        ordChroms <- c(ordChroms, ordChroms_nonNum)
        
        gr_lifted_ord <- data.frame(matrix(nrow = 0,
                                           ncol = ncol(gr_lifted),
                                           dimnames = list(NULL,
                                                           colnames(gr_lifted))))
        for(chr in ordChroms){
                chr_df <- gr_lifted[gr_lifted$seqnames == chr, ]
                chr_df <- chr_df[order(chr_df$start), ]
                gr_lifted_ord <- rbind.data.frame(gr_lifted_ord, chr_df)
        }
        gr_lifted <- gr_lifted_ord
}

# Save result
################################################################################

writeCsvFst(gr_lifted, rowNames = F, colNames = F, file = outName, sep = "\t")
if(compress){
        system(sprintf("gzip %s", outName))
        outName <- paste0(outName, ".gz")
}
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))
