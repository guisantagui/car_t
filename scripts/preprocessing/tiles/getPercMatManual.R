################################################################################
# CAR-T: given a the tiled RDS object obtained with the bismark files, obtains #
# methylation percentage without the need of uniting the regions that are      #
# found across all the samples. If a region is not present in one sample it is #
# assigned to a NA.                                                            #
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

library(methylKit)
library(argparser)

# Parser
################################################################################
parser <- arg_parser("Given a the tiled RDS object obtained with the bismark
# files, obtains methylation percentage without the need of uniting the regions
# that are found across all the samples. If a region is not present in one
# sample it is assigned to a NA.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--outName"),
                       help = c("Rds file with the tiled samples, but without being united.",
                                "Name for the resulting CSV file."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Add chr to chromosome section of regionIDs that are autosome (i.e. they are
# a number).
chrInAutosome <- function(regIDs){
        allChrs <- gsub("\\..*", "", regIDs)
        isAutosome <- !is.na(suppressWarnings(as.numeric(regIDs)))
        regIDs[isAutosome] <- paste0("chr", regIDs[isAutosome])
        return(regIDs)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

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

# Directory stuff
################################################################################
tileFile <- parsed$input
outName <- parsed$outName

outDir <- dirname(outName)
createIfNot(outDir)

# Load data
################################################################################
tiles <- readRDS(tileFile)


# Integrate samples
################################################################################

# Obtain all unique tiles across all the samples
allTiles <- lapply(tiles,
                   function(x) paste(getData(x)$chr,
                                     paste(getData(x)$start,
                                           getData(x)$end,
                                           sep = "_"),
                                     sep = "."))

allTiles <- unique(unlist(allTiles))

# Add chr to autosomes
allTiles<- chrInAutosome(allTiles)

# Obtain the methylation percentage of all tiles in each sample
print("Integrating tiled samples into a single methylation percentage matrix...")
methMat <- data.frame(matrix(nrow = 0,
                             ncol = length(allTiles),
                             dimnames = list(NULL,
                             allTiles)))
pb <- txtProgressBar(min = 0, max = length(tiles), initial = 0, style = 3)
for(i in seq_along(tiles)){
        setTxtProgressBar(pb, i)
        samp <- getSampleID(tiles[[i]])
        tileDF <- getData(tiles[[i]])
        tileDF$tileID <- paste(tileDF$chr,
                               paste(tileDF$start,
                                     tileDF$end,
                                     sep = "_"),
                               sep = ".")
        # Add chr to autosomes
        tileDF$tileID <- chrInAutosome(tileDF$tileID)
        # Compute methylation percentage of each region.
        tileDF$methPerc <- tileDF$numCs/tileDF$coverage * 100
        methPercVec <- tileDF$methPerc[match(allTiles,
                                             tileDF$tileID)]
        toBind <- data.frame(matrix(methPercVec,
                                    nrow = 1,
                                    ncol = length(allTiles),
                                    dimnames = list(samp,
                                                    allTiles)))
        methMat <- rbind.data.frame(methMat, toBind)
}
close(pb)

# Integrate the methylation percentages into a single dataset.
methMat <- do.call(rbind.data.frame, methMat)

writeCsvFst(methMat, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))
