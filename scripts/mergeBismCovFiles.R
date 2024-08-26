################################################################################
# CAR-T: Add up bismark cov files in a new one. Intended to use to sum up      #    
# the counts of the new reseq data of 2024 analyses.                           #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Add up two bismark cov files.")

parser <- add_argument(parser = parser,
                       arg = c("--cov1",
                               "--cov2",
                               "--outName"),
                       help = c("COV file 1",
                                "COV file 2",
                                "path/name of the resulting file"),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

cov1_file <- parsed$cov1
cov2_file <- parsed$cov2
outName <- parsed$outName


if(!dir.exists(dirname(outName))){
        dir.create(dirname(outName), recursive = T)
}

# Functions
################################################################################

readBismCov <- function(f){
        dat <- data.frame(data.table::fread(f))
        colnames(dat) <- c("chr",
                           "start",
                           "end",
                           "perc_meth",
                           "count_meth",
                           "count_unmeth")
        return(dat)
}


# Load data
################################################################################
cov_1 <- readBismCov(cov1_file)
cov_2 <- readBismCov(cov2_file)

# Add up datasets
################################################################################

print("Finding unique positions in:")
chroms <- sort(unique(c(cov_1$chr, cov_2$chr)))
uniq_pos <- c()
strt_vec <- c()
chr_vec <- c()
for(i in seq_along(chroms)){
        chr <- chroms[i]
        print(sprintf("Chromosome %s", chr))
        positions <- sort(unique(c(cov_1$start[cov_1$chr == chr],
                                   cov_2$start[cov_2$chr == chr])))
        chr_vec <- c(chr_vec, rep(chr, length(positions)))
        strt_vec <- c(strt_vec, positions)
}

print("Merging COV files...")
comb_cov <- data.frame(chr = chr_vec,
                       start = strt_vec,
                       end = strt_vec,
                       perc_meth = rep(NA, length(chr_vec)),
                       count_meth = rep(0, length(chr_vec)),
                       count_unmeth = rep(0, length(chr_vec)))

match_1 <- match(paste(cov_1$chr, cov_1$start, sep = "_"),
                 paste(comb_cov$chr, comb_cov$start, sep = "_"))

match_2 <- match(paste(cov_2$chr, cov_2$start, sep = "_"),
                 paste(comb_cov$chr, comb_cov$start, sep = "_"))

comb_cov[match_1, c("count_meth", "count_unmeth")] <- comb_cov[match_1, c("count_meth", "count_unmeth")] + cov_1[, c("count_meth", "count_unmeth")]
comb_cov[match_2, c("count_meth", "count_unmeth")] <- comb_cov[match_2, c("count_meth", "count_unmeth")] + cov_2[, c("count_meth", "count_unmeth")]
comb_cov$perc_meth <- comb_cov$count_meth/(comb_cov$count_meth + comb_cov$count_unmeth) * 100

# Save and compress resulting dataset
################################################################################

outName <- gsub(".gz", "", outName)
write.table(comb_cov, col.names = F, row.names = F, quote = F, file = outName, sep = "\t")
print("Writing file...")
system(sprintf("gzip -f %s && rm %s", outName, outName))
print("Compressing file...")
print(sprintf("%s.gz saved at %s", basename(outName), dirname(outName)))