################################################################################
# CAR-T: given a two tiled RDS object obtained with the bismark files, merges  #
# them into a single one.                                                      #
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
parser <- arg_parser("Given a two tiled RDS object obtained with the bismark files, merges  #
# them into a single one.")

parser <- add_argument(parser = parser,
                       arg = c("--file1",
                               "--file2",
                               "--outName"),
                       help = c("Rds file1 with the tiled samples.",
                                "Rds file2 with the tiled samples.",
                                "Name for the resulting RDS file."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Directory stuff
################################################################################

file1 <- "/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_2023/car_t_2023_MK_100bp_tiled.rds"
file2 <- "/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_2024/car_t_2024_new_MK_100bp_tiled.rds"
outName <- "/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_comb/car_t_comb_100bp_tiled.rds"

file1 <- parsed$file1
file2 <- parsed$file2
outName <- parsed$outName

outDir <- dirname(outName)
createIfNot(outDir)
# Load data
################################################################################
dat1 <- readRDS(file1)
dat2 <- readRDS(file2)

# Combine datasets and save result
################################################################################
dat_comb <- c(dat1, dat2)

saveRDS(dat_comb, outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))
