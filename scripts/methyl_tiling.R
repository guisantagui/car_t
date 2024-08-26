if(!require("BiocManager", quietly = T)) install.packages("BiocManager",
                                                          repos = "https://pbil.univ-lyon1.fr/CRAN/")
if(!require("readxl", quietly = T)) install.packages("readxl",
                                                     repos = "https://pbil.univ-lyon1.fr/CRAN/")
library(readxl)
if(!require("methylKit", quietly = T)) BiocManager::install("methylKit",
                                                            update = F)

library(parallel)
library(methylKit)

# Directory stuff
################################################################################

rootDir <- "/home/users/gsantamaria/projects/car_t/"
tileSize <- 100
minCov <- 5
minCpGs_inTile <- 5
uncomp <- F


dataDir <- paste0(rootDir, "data/")
wgbsDir <- paste0(dataDir, "wgbs/")
resuDir <- paste0(rootDir, "results/")
preprocDir <- paste0(resuDir, "preprocessing/")

wgbsFiles <- list.files(wgbsDir)[!R.utils::isDirectory(list.files(wgbsDir,
                                                                  full.names = T))]
print(wgbsFiles)

createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

createIfNot(preprocDir)

# Load and uncompress coverage data
################################################################################

# Load methylation metadata
metDat <- data.frame(read_xlsx(paste0(dataDir,
                                      "Tcells2023WGBSMeta_Links_collab.xlsx")))

# Add a vector indicating if the sample is good, based on the colors of the 
# excel file (as read_xlsx can't see the colors)
metDat$isGood <- rep(T, nrow(metDat))
metDat$isGood[c(4, 11, 18)] <- F

# Uncompress files 

if(uncomp){
        for(f in wgbsFiles){
                fSplit <- strsplit(f, split = ".", fixed = T)[[1]]
                fUncomp <- paste(fSplit[1:(length(fSplit) - 1)], collapse = ".")
                cmd <- sprintf("gzip -d %s%s", wgbsDir, f)
                if(!file.exists(sprintf("%s%s", wgbsDir, fUncomp))){
                        system(cmd)
                }
        }
        # Keep the names of the uncompressed COV files
        wgbsFilesUncomp <- list.files(wgbsDir)
        wgbsFilesUncomp <- wgbsFilesUncomp[!grepl(".cov.gz",
                                                  wgbsFilesUncomp,
                                                  fixed = T)]
}else{
        wgbsFilesUncomp <- list.files(wgbsDir)
        #wgbsFilesUncomp <- wgbsFilesUncomp[grepl(".cov.gz",
        #                                         wgbsFilesUncomp,
        #                                         fixed = T)]
}



# Keep only the files that are good quality (based on not being red in the 
# metadata file)
wgbsFilesUncomp <- wgbsFilesUncomp[sapply(wgbsFilesUncomp,
                                          function(x) any(grepl(x,
                                                                basename(metDat$WGBS_Cov)[metDat$isGood])))]

# Create a sample info file
sampInfo <- metDat[metDat$isGood, ]

sampInfo$WGBS_Cov <- sapply(basename(sampInfo$WGBS_Cov),
                            function(x) paste(strsplit(x,
                                                       ".",
                                                       fixed = T)[[1]][1:(length(strsplit(x,
                                                                                          ".",
                                                                                          fixed = T)[[1]]) - 1)],
                                              collapse = "."))

# Do the tiling with a given length of BPs
################################################################################

# Load Bismark coverage files
wgbsFilesUncomp_list <- as.list(sprintf("%s%s", wgbsDir, wgbsFilesUncomp))
sampList <- as.list(sampInfo$Sample[match(wgbsFilesUncomp, sampInfo$WGBS_Cov)])
print(sprintf("Loading WGBS files from %s...", wgbsDir))
bismCovListRaw <- methRead(wgbsFilesUncomp_list,
                           sample.id = sampList,
                           #assembly = as.list(wgbsFilesUncomp),
                           assembly = "hg18",
                           pipeline = "bismarkCoverage",
                           context = "CpG",
                           mincov = minCov,
                           treatment = as.numeric(!grepl("M", unlist(sampList))))

# Do a tile analysis
nCores <- detectCores() - 2

print(sprintf("Performing methylation dataset integration with tile size of %s bp...",
              as.character(tileSize)))

tileInteg <- tileMethylCounts(bismCovListRaw,
                              win.size = tileSize,
                              step.size = tileSize,
                              cov.bases = minCpGs_inTile,
                              mc.cores = nCores)

saveRDS(tileInteg,
        file = sprintf("%swgbs_%sbp_tileInt.rds",
                       preprocDir,
                       as.character(tileSize)))

print(sprintf("wgbs_%sbp_tileInt.rds saved at %s.",
              as.character(tileSize),
              preprocDir))

print(sprintf("Uniting the %s bp tiles across samples...",
              as.character(tileSize)))

tileInteg_unit <- unite(tileInteg)

tileInteg_unit_dat <- getData(tileInteg_unit)

tileNames <- paste(paste(paste0("chr", tileInteg_unit_dat[, 1]),
                         tileInteg_unit_dat[, 2], sep = ":"),
                   tileInteg_unit_dat[, 3], sep = "-")



percMat <- percMethylation(tileInteg_unit)

rownames(percMat) <- tileNames

percMat <- t(percMat)

write.csv(percMat,
          sprintf("%swgbs_%sbp_percMeth.csv",
                  preprocDir,
                  as.character(tileSize)))

saveRDS(percMat,
        sprintf("%swgbs_%sbp_percMeth.rds",
                  preprocDir,
                  as.character(tileSize)))

print(sprintf("wgbs_%sbp_percMeth.rds saved at %s.",
              as.character(tileSize),
              preprocDir))