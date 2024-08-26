################################################################################
# CAR-T: Given a bismark wgbs file and a LMR file obtained by hmr (DNMTools),  #
# obtain a dataframe indicating methylation percentage, and other information, #
# for the introduced sample on each LMR.                                       #
################################################################################

if(!require("BiocManager", quietly = T)){
        install.packages("BiocManager",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("argparser", quietly = T)){
        install.packages("argparser",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("methylKit", quietly = T)) BiocManager::install("methylKit",
                                                            update = F)

library(parallel)
library(methylKit)
library(argparser)
library(dplyr)

# Parser
################################################################################
parser <- arg_parser("Given a bismark wgbs file and a LMR file obtained by hmr (DNMTools),
obtain a dataframe indicating methylation percentage, and other information,
for the introduced sample on each LMR.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--lmrFile",
                               "--minCov",
                               "--outDir"),
                       help = c("Bismark file.",
                                "LMR hmr file.",
                                "Minimum coverage for considering a CpG.",
                                "Output directory."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# write.csv, but faster
writeCsvFst <- function(df, file, rowNames = T, colNames = T){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames)
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

# Directory stuff
################################################################################

bismFile <- parsed$input
lmrFile <- parsed$lmrFile
minCov <- as.numeric(parsed$minCov)
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

# Load data
################################################################################

lmrs <- read.table(lmrFile)
colnames(lmrs) <- c("chr", "start", "end", "id", "n_CpGs", "strand")

samp <- basename(bismFile)

if(grepl("_", samp)){
        samp <- strsplit(samp, split = "_")[[1]]

        samp <- paste(samp[1:min(grep("bismark", samp))], collapse = "_")
}else{
        samp <- strsplit(samp, split = ".", fixed = T)[[1]][1]
}


bism <- methRead(bismFile,
                 sample.id = samp,
                 assembly = "h38",
                 pipeline = "bismarkCoverage",
                 context = "CpG",
                 mincov = minCov)


# Get LMRs for the introduced sample
################################################################################
print(sprintf("Parsing LMRs in %s for sample %s...", basename(lmrFile), samp))

lmrDF <- list()

chrs <- unique(lmrs$chr)
chrs <- gsub("chr", "", chrs)
pb <- txtProgressBar(min = 0, max = length(chrs), initial = 0, style = 3)
for(i in seq_along(chrs)){
        setTxtProgressBar(pb, i)
        chr <- chrs[i]
        bism_chr <- bism[bism$chr == chr, ]
        bism_chr <- getData(bism_chr)
        lmrs_chr <- lmrs[lmrs$chr == paste0("chr", chr), ]
        dfList <- apply(lmrs_chr,
                        1,
                        function(x) bism_chr[bism_chr$start >= as.numeric(x[2]) & bism_chr$end <= as.numeric(x[3]), ])
        
        toBind <- lmrs_chr
        colnames(toBind) <- gsub("n_CpGs", "n_CpGs_inRef", colnames(toBind))
        toBind$nCpGs_inSamp <- unlist(lapply(dfList, nrow))
        
        # percMethylation function from methylKit obtains the percentage of methylation
        # of the tiles by dividing the number of total number of Cs of the CpGs
        # within a tile by the total coverage of the CpGs in the tile. So let's
        # do the same with the LMRs.
        toBind$percMethAll <- lapply(dfList, function(x) sum(x$numCs)/(sum(x$coverage)) * 100) %>% unlist
        
        lmrDF[[i]] <- toBind
}

lmrDF <- do.call(rbind, lmrDF)

outName <- sprintf("%s%s_LMRs.csv", outDir, samp)
writeCsvFst(lmrDF, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))