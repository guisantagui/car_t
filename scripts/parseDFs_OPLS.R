################################################################################
#                                                                              #
# Parse CAR-T cell metadata and percentage of methylation datato have a        #
# version to work with WGBS OPLS-DA RFE (removing samples that are not in the  #
# WGBS data and mocks at zero).                                                #
#                                                                              #
################################################################################

library(readxl)
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("RFE for omic data classification")

parser <- add_argument(parser = parser,
                       arg = c("--metDatFile",
                               "--exhFile",
                               "--percMethFile",
                               "--outDir"),
                       help = c("metadata XLSX file",
                                "Exhaustion score CSV file",
                                "Percentage of methylation CSV file",
                                "Output directory"),
                       flag = c(F, F, F, F))

# Directory stuff
################################################################################
#metDat <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/data/Tcells2023WGBSMeta_Links_collab.xlsx"
#exhFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/results/exh_measure/exh_score.csv"
#percMethFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/results/preprocessing/percMeth.csv"

#outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/results/preprocessing/"

metDat <- parsed$metDatFile
exhFile <- parsed$exhFile
percMethFile <- parsed$percMethFile

outDir <- parsed$outDir

# Load data 
################################################################################

metDat <- data.frame(read_xlsx(metDat))
#metDat <- metDat[ c("Sample", "Donor", "CAR", "Time")]

# This column is for indicating the samples where bisulphite conversion didn't
# work well.
metDat$isGood <- rep(T, nrow(metDat))
metDat$isGood[c(4, 11, 18)] <- F


exh <- read.csv(exhFile, row.names = 1)
percMeth <- read.csv(percMethFile, row.names = 1)

# Parse datasets
################################################################################

sampInfo <- metDat
sampInfo$exh_score <- exh$exh_score[match(sampInfo$Sample, exh$sample)]

sampInfo$WGBS_Cov <- sapply(basename(sampInfo$WGBS_Cov),
                            function(x) paste(strsplit(x,
                                                       ".",
                                                       fixed = T)[[1]][1:(length(strsplit(x,
                                                                                          ".",
                                                                                          fixed = T)[[1]]) - 1)],
                                              collapse = "."))

percMeth <- percMeth[rownames(percMeth) %in% sampInfo$Sample[!is.na(sampInfo$exh_score)], ]

colnames(sampInfo) <- gsub("Sample", "sample", colnames(sampInfo))

# Write files to the output directory
################################################################################

write.csv(sampInfo, file = paste0(outDir, "sample_info.csv"))

percMethOutName <- sprintf("%s%s_4OPLS.csv", outDir, gsub(".csv", "", basename(percMethFile)))
write.csv(percMeth, file = percMethOutName)