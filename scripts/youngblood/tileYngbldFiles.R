################################################################################
# CAR-T: Tiles the Youngblood files with methylkit.                            #
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

bF <- "/home/users/gsantamaria/projects/car_t/data/youngblood_hg38Lift/GSM5677818_2-1245637_hg38.txt.gz"
bed <- readCsvFst(bF, header = F, row.names = F)
colnames(bed) <- c("chr", "start", "end", "width", "strand", "coverage", "numCs", "numTs")

mRaw <- methRead()