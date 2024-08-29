################################################################################
# CAR-T: given a folder where WGBS Bismark COV files are, returns a tiled      #
# object with the specified tile size and a matrix of methylation percentage   #
# per tile                                                                     #
################################################################################


if(!require("BiocManager", quietly = T)){
        install.packages("BiocManager",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("readxl", quietly = T)){
        install.packages("readxl",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("argparser", quietly = T)){
        install.packages("argparser",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}
if(!require("methylKit", quietly = T)) BiocManager::install("methylKit",
                                                            update = F)

library(readxl)
library(parallel)
library(methylKit)
library(argparser)

# Parser
################################################################################
parser <- arg_parser("Given a directory of WGBS Bismark files, generates tiled methylKit object and percentage of methylation per tile.")

parser <- add_argument(parser = parser,
                       arg = c("--wgbsDir",
                               "--tileSize",
                               "--minCov",
                               "--minCpGs_inTile",
                               "--outName",
                               "--remPrevFile",
                               "--outDir"),
                       help = c("Directory where WGBS COV files are.",
                                "Size of the tiles.",
                                "Minimum coverage for considering a CpG.",
                                "Minimum number of CpGs within the tile for considering it.",
                                "Name to include in the output.",
                                "If previously generated united RDS file should be removed in order to uniting it again.",
                                "Output directory."),
                       flag = c(F, F, F, F, F, T, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Directory stuff
################################################################################

wgbsDir <- addSlashIfNot(parsed$wgbsDir)
tileSize <- as.numeric(parsed$tileSize)
minCov <- as.numeric(parsed$minCov)
minCpGs_inTile <- as.numeric(parsed$minCpGs_inTile)
outName <- parsed$outName
remPrevFile <- parsed$remPrevFile
outDir <- addSlashIfNot(parsed$outDir)


createIfNot(outDir)

# Get filenames for COV files
################################################################################

wgbsFiles <- list.files(wgbsDir, full.names = T)
wgbsFiles <- wgbsFiles[grepl(".cov.gz",
                             wgbsFiles,
                             fixed = T)]

#wgbsFiles <- wgbsFiles[!grepl("CD8-Mock-8-69.deduplicated.bismark.cov",
#                              wgbsFiles)]

# Do the tiling with a given length of BPs
################################################################################

MK_tileName <- sprintf("%s%s_MK_%sbp_tiled.rds",
                       outDir,
                       outName,
                       as.character(tileSize))

MK_tileUnitName <- sprintf("%s%s_MK_%sbp_tiledUnited.rds",
                           outDir,
                           outName,
                           as.character(tileSize))

if(remPrevFile & file.exists(MK_tileName) & file.exists(MK_tileUnitName)){
        print(sprintf("Removing %s and %s from %s, generated in a previous run.",
                      basename(MK_tileName),
                      basename(MK_tileUnitName),
                      dirname(MK_tileName)))
        rmComm <- sprintf("rm %s", MK_tileName)
        system(rmComm)
        rmComm <- sprintf("rm %s", MK_tileUnitName)
        system(rmComm)
}

if(!file.exists(MK_tileName)){
        # Load Bismark coverage files
        wgbsFiles_list <- as.list(wgbsFiles)

        samps <- gsub("_bismark_bt2_pe",
                      "",
                      gsub("\\..*",
                           "",
                           basename(wgbsFiles)))

        sampList <- as.list(samps)
        print(sprintf("Loading WGBS files from %s...", wgbsDir))
        
        

        bismCovListRaw <- methRead(wgbsFiles_list,
                                   sample.id = sampList,
                                   assembly = "hg38",
                                   pipeline = "bismarkCoverage",
                                   context = "CpG",
                                   mincov = minCov,
                                   treatment = as.numeric(!grepl("M",
                                                                 unlist(sampList))))

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
                file = MK_tileName)

        print(sprintf("%s saved at %s.",
                      basename(MK_tileName),
                      dirname(MK_tileName)))

        print(sprintf("Uniting the %s bp tiles across samples...",
                      as.character(tileSize)))

        tileInteg_unit <- unite(tileInteg)

        saveRDS(tileInteg_unit,
                file = MK_tileUnitName)

        print(sprintf("%s saved at %s.",
                      basename(MK_tileUnitName),
                      dirname(MK_tileUnitName)))
}else{  
        print(sprintf("Loading %s from %s...",
                      basename(MK_tileName),
                      dirname(MK_tileName)))
        tileInteg_unit <- readRDS(MK_tileName)
}

tileInteg_unit_dat <- getData(tileInteg_unit)

chrVec <- as.character(tileInteg_unit_dat[, 1])
strtVec <- as.character(tileInteg_unit_dat[, 2])
endVec <- as.character(tileInteg_unit_dat[, 3])

chrVec[nchar(chrVec) < 3] <- paste0("chr",
                                    chrVec[nchar(chrVec) < 3])

tileNames <- paste(paste(chrVec,
                         strtVec, sep = "."),
                   endVec, sep = "_")



percMat <- percMethylation(tileInteg_unit)

rownames(percMat) <- tileNames

percMat <- t(percMat)

percMatName <- sprintf("%s%s_%sbp_percMeth.csv",
                       outDir,
                       outName,
                       as.character(tileSize))

write.csv(percMat,
          percMatName)

print(sprintf("%s saved at %s.",
              basename(percMatName),
              dirname(percMatName)))