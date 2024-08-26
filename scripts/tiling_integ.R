if(!require("BiocManager", quietly = T)) install.packages("BiocManager",
                                                          repos = "https://pbil.univ-lyon1.fr/CRAN/")

if(!require("methylKit", quietly = T)) BiocManager::install("methylKit",
                                                            update = F)

library(methylKit)

# Directory stuff
################################################################################
rootDir <- "/home/users/gsantamaria/projects/car_t/"

dataDir <- paste0(rootDir, "data/")
resuDir <- paste0(rootDir, "results/")
prepDir <- paste0(resuDir, "preprocessing/")


# Load tile dataset
tileFile <- paste0(prepDir, "wgbs_100bp_tileInt.rds")

tileObj <- readRDS(tileFile)

# Look waht positions are found across all the samples (some were left out
# in some samples due to bad coverage) and generate a matrix of methylation
# percentage of the common positions.
posList <- lapply(tileObj,
                  function(x) paste(getData(x)$chr,
                                    getData(x)$start,
                                    getData(x)$end,
                                    sep = "_"))

allPositions <- sort(unique(unlist(posList)))

sampNames <- getSampleID(tileObj)
tileIntPercDF <- data.frame(matrix(nrow = length(tileObj),
                                   ncol = 0,
                                   dimnames = list(sampNames, NULL)))

print("Generating integrated dataframe...")
for(i in seq_along(allPositions)){
        p <- allPositions[i]
        isInAll <- all(unlist(lapply(posList,
                                     function(x) p %in% x)))
        if(isInAll){
                posPercVec <- c()
                for(j in seq_along(tileObj)){
                        samp <- getSampleID(tileObj[[j]])
                        sampTile <- getData(tileObj[[j]])
                        posVec <- paste(sampTile$chr,
                                        sampTile$start,
                                        sampTile$end,
                                        sep = "_")
                        sampSlice <- sampTile[posVec == p, ]
                        sampPerc <- sampSlice$numCs/sampSlice$coverage * 100
                        posPercVec <- c(posPercVec, sampPerc)
                }
                toBind <- data.frame(matrix(posPercVec,
                                            nrow = length(posPercVec),
                                            ncol = 1,
                                            dimnames = list(sampNames, p)))
                tileIntPercDF <- cbind.data.frame(tileIntPercDF, toBind)
                #print(tileIntPercDF)
        }
        percProc <- i/length(allPositions) * 100
        percProc <- as.character(round(percProc, digits = 3))
        print(sprintf("%s %%", percProc))
}

# Save the resulting dataframe
write.csv(tileIntPercDF, file = sprintf("%sintegTile100bp.csv", prepDir))
print(sprintf("The number fo common positions found across all the samples is %s.",
              as.character(ncol(tileIntPercDF))))
print(sprintf("integTile100bp.csv saved at %s", prepDir))
