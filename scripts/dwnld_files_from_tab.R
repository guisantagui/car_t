################################################################################
# CAR-T project: given a table in txt or xlsx format, downloads the files and  #
# generates a sample_info CSV file.                                            #
################################################################################

if(!require("argparser", quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
library(readxl)

parser <- arg_parser("This script downloads files given a file with the URLs (xlsx or txt, separated with tabs)")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--expName",
                               "--outDir"),
                       help = c("Directory where toRem files are stored",
                                "Name for the experiment (for attaching it to sample info file).",
                                "Path to directory where the downloaded files will be kept."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)



# Directory stuff
################################################################################

inFile <- parsed$input
expName <- parsed$expName
outDir <- parsed$outDir

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
        print(sprintf("%s was created.", outDir))
}

# This prevents that in case that the output directory name is not ended by "/"
# the output name is not generated correctly.
if(substr(outDir, nchar(outDir), nchar(outDir)) != "/"){
        outDir <- paste0(outDir, "/")
}

# Parse input table for the sample_info dataframe and to obtain URLs
################################################################################


if(grepl(".xlsx", inFile)){
        metaDF <- read_xlsx(inFile)
        
        metaDF <- data.frame(metaDF)
        metaDF$T_cell <- rep("Mixed", nrow(metaDF))
        metaDF$exp_name <- rep(expName, nrow(metaDF))
        metaDF$isGood <- rep(T, nrow(metaDF))
        metaDF$isGood[c(4, 11, 18)] <- F
        metaDF <- metaDF[, c("Sample", "Donor", "CAR", "T_cell",
                             "Time", "exp_name", "isGood", "WGBS_Cov")]
        URLs <- metaDF$WGBS_Cov
}else if(grepl(".txt", inFile)){
        URLs <- read.table(inFile, sep = "\t")$V1
        samps <- gsub("\\..*", "", basename(URLs))
        
        metaDF <- data.frame(Sample = samps,
                             Donor = paste0("d",
                                            sapply(samps,
                                                   function(x) strsplit(x,
                                                                        split = "-")[[1]][4])),
                             CAR = sapply(samps,
                                          function(x) strsplit(x,
                                                               split = "-")[[1]][2]),
                             T_cell = sapply(samps,
                                             function(x) strsplit(x,
                                                                  split = "-")[[1]][1]),
                             Time = sapply(samps,
                                           function(x) strsplit(x,
                                                                split = "-")[[1]][3]),
                             exp_name = rep(expName, length(samps)),
                             isGood = rep(T, length(samps)),
                             WGBS_Cov = URLs)
}



for(link in URLs){
        cmd <- sprintf("wget %s -P %s", link, outDir)
        fileNam <- basename(link)
        if(!file.exists(sprintf("%s%s", outDir, fileNam))){
                system(cmd)
        }else{
                print(sprintf("%s is already in %s", fileNam, outDir))
        }
}
samp_info_name <- sprintf("%ssamp_info_%s.csv", outDir, expName)

write.csv(metaDF, file = samp_info_name)
print(sprintf("%s saved in %s", samp_info_name, outDir))