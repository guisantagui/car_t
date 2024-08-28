################################################################################
# CAR-T: Given a dataframe of genomic regions and a directory containing       #
# methylation bed files per sample, returns a matrix of methylation percentage #
# per region for each sample. Intended to be used with Youngblood bed files.   #
################################################################################

if(!require("argparser", quietly = T)){
        install.packages("argparser",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}

library(argparser)
library(dplyr)

# Parser
################################################################################
parser <- arg_parser("Given a dataframe of genomic regions and a directory
containing methylation bed files per sample, returns a matrix of methylation
percentage per region for each sample. Intended to be used with Youngblood bed
files.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--regs",
                               "--minCov",
                               "--minNumCpGs",
                               "--ref",
                               "--outName"),
                       help = c("Directory where bed files to be integrated are located.",
                                "File with regions to be obtained are. Has to have at least a chr, start and end columns.",
                                "Minimum coverage to consider a CpG.",
                                "Minimum number of CpGs to consider a region.",
                                "Reference genome that was used to generate the bed files.",
                                "Output name of the resulting file."),
                       flag = c(F, F, F, F, F, F))

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
bedDir <- parsed$input
regionFile <- parsed$regs
minCov <- as.numeric(parsed$minCov)
minNumCpGs <- as.numeric(parsed$minNumCpGs)
ref <- parsed$ref
outName <- parsed$outName
outDir <- dirname(outName)


bedDir <- addSlashIfNot(bedDir)
outDir <- addSlashIfNot(outDir)

createIfNot(outDir)

# Load data
################################################################################

regs <- readCsvFst(regionFile)
bedFiles <- list.files(bedDir, full.names = T)
bedFiles <- bedFiles[!grepl(".csv", bedFiles)]

# Get methylation percentage of the desired samples
################################################################################
chrs <- unique(regs$chr)
methPercRegsDF <- list()
for(j in seq_along(bedFiles)){
        bF <- bedFiles[j]
        bed <- readCsvFst(bF, row.names = F, header = F)
        colnames(bed) <- c("seqnames",
                           "start",
                           "end",
                           "width",
                           "strand",
                           "coverage",
                           "numCs",
                           "numTs")
        sampName <- basename(bF)
        sampName <- gsub(".gz", "", sampName)
        sampName <- gsub(".txt", "", sampName)
        sampName <- gsub(sprintf("_%s", ref), "", sampName)
        print(sprintf("Parsing regions in in %s for sample %s...",
                      basename(regionFile),
                      sampName))
        sampMethPercDF <- list()
        #chrs <- gsub("chr", "", chrs)
        pb <- txtProgressBar(min = 0, max = length(chrs), initial = 0, style = 3)
        for(i in seq_along(chrs)){
                #i <- 1
                setTxtProgressBar(pb, i)
                chr <- chrs[i]
                bed_chr <- bed[bed$seqnames == chr, ]
                regs_chr <- regs[regs$chr == chr, ]
                dfList <- apply(regs_chr,
                                1,
                                function(x) bed_chr[bed_chr$start >= as.numeric(x[2]) & bed_chr$end <= as.numeric(x[3]), ])
                
                # Filter out CpGs with less than the indicated minimum coverage
                dfList <- lapply(dfList, function(x) x[x$coverage >= minCov, ])
                
                # Remove dataframes that have less than the minimum number of CpGs 
                # for considering a region
                dfList <- lapply(dfList, function(x){
                        if(nrow(x) < minNumCpGs){
                                x[rep(F, nrow(x)), ]
                        }else{
                                x
                        }
                })
                
                toBind <- regs_chr[, c("chr", "start", "end")]
                
                # percMethylation function from methylKit obtains the percentage of methylation
                # of the tiles by dividing the number of total number of Cs of the CpGs
                # within a tile by the total coverage of the CpGs in the tile. So let's
                # do the same with the LMRs.
                toBind$percMethAll <- lapply(dfList,
                                             function(x) sum(x$numCs)/(sum(x$coverage)) * 100) %>% unlist
                
                sampMethPercDF[[i]] <- toBind
        }
        
        sampMethPercDF <- do.call(rbind, sampMethPercDF)
        
        sampMethPercDF$region <- paste(sampMethPercDF$chr,
                                       paste(sampMethPercDF$start,
                                             sampMethPercDF$end,
                                             sep = "_"),
                                       sep = ".")
        toBindBig <- data.frame(matrix(sampMethPercDF$percMethAll,
                                       ncol = length(sampMethPercDF$percMethAll),
                                       nrow = 1,
                                       dimnames = list(sampName,
                                                       sampMethPercDF$region)))
        toBindBig[apply(toBindBig, 2, is.nan)] <- NA
        methPercRegsDF[[j]] <- toBindBig
}

methPercRegsDF <- do.call(rbind, methPercRegsDF)
methPercRegsDF

write.csv(methPercRegsDF, outName)
print(sprintf("%s saved in %s.", basename(outName), dirname(outName)))