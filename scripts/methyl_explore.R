
if(!require("BiocManager", quietly = T)) install.packages("BiocManager")
if(!require("readxl", quietly = T)) install.packages("readxl")
library(readxl)
if(!require("ggplot2", quietly = T)) install.packages("ggplot2")
library(ggplot2)
if(!require("ggpubr", quietly = T)) install.packages("ggpubr")
library(ggpubr)
if(!require("ggrepel", quietly = T)) install.packages("ggrepel")
library(ggrepel)
if(!require("methylKit", quietly = T)) BiocManager::install("methylKit",
                                                            update = F)
library(parallel)
library(methylKit)

rootDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/car_t/"

dataDir <- paste0(rootDir, "data/")
wgbsDir <- paste0(dataDir, "wbgs/")

wgbsFiles <- list.files(wgbsDir)
wgbsFiles <- wgbsFiles[grep(".cov.gz", wgbsFiles)]

plotDir <- paste0(rootDir, "plots/")

createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

createIfNot(plotDir)

# Functions
################################################################################

# Define a function for obtaining the plot
pcBiplot <- function(PC, x="PC1", y="PC2", varPlotFilt = NULL, biPlot = T){
        #x <- x
        #y <- y
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        data <- data[, c("obsnames", x, y)]
        data$time <- sampInfo$Time[match(rownames(data),
                                         sampInfo$Sample)]
        data$time <- factor(as.character(data$time), levels = c("0", "8", "38"))
        data$car <- sampInfo$CAR[match(rownames(data),
                                       sampInfo$Sample)]
        data$car <- factor(data$car, levels = c("Mock", "CD19", "HA"))
        data$donor <- sampInfo$Donor[match(rownames(data), sampInfo$Sample)]
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        
        plot <- ggplot(data, aes(x = data[, x], 
                                 y = data[, y], 
                                 label = donor, 
                                 color = time,
                                 shape = car)) + 
                scale_discrete_manual("Time",
                                      aesthetics = "colour",
                                      values = c("cyan",
                                                 "dodgerblue",
                                                 "dodgerblue4")) +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) +
                geom_text_repel() +
                #theme_minimal()
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if(biPlot){
                datapc <- data.frame(varnames=rownames(PC$rotation), 
                                     PC$rotation)
                mult <- min(
                        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
                )
                datapc <- transform(datapc,
                                    v1 = .7 * mult * (get(x)),
                                    v2 = .7 * mult * (get(y))
                )
                datapc$x0 <- rep(0, nrow(datapc))
                datapc$y0 <- rep(0, nrow(datapc))
                if(!is.null(varPlotFilt)){
                        datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
                }
                plot <- plot +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3) + 
                        geom_segment(data = datapc, aes(x=x0, 
                                                        y=y0, 
                                                        xend=v1, 
                                                        yend=v2, 
                                                        label = varnames),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20), 
                                     alpha=0.75, 
                                     color="black", 
                                     size = 0.5)
        }
        return(plot)
}

# Plots Multi PC score/biplot
doPCAMultiPlot <- function(PC, nComps, varPlotFilt = NULL, biPlot = F){
        plotList <- list()
        for(j in 2:(nComps)){
                for(i in 1:(nComps - 1)){
                        if(j > i){
                                scPlot <- pcBiplot(PC,
                                                   biPlot = F,
                                                   x = sprintf("PC%s", i),
                                                   y = sprintf("PC%s", j),
                                                   varPlotFilt = varPlotFilt)
                                if(j < nComps){
                                        scPlot <- scPlot +
                                                theme(axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                                }
                                if(i > 1){
                                        scPlot <- scPlot +
                                                theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank())
                                }
                        }else{
                                scPlot <- NA
                        }
                        plotList[[sprintf("PC%s_PC%s", i, j)]] <- scPlot
                }
        }
        multPlot <- ggarrange(plotlist = plotList,
                              common.legend = T,
                              ncol = nComps - 1,
                              nrow = nComps - 1,
                              widths = c(1, rep(.8, nComps-2)),
                              heights = c(rep(.8, nComps-2), 1))
        return(multPlot)
}


# Load methylation metadata
metDat <- data.frame(read_xlsx(paste0(dataDir, "Tcells2023WGBSMeta_Links_collab.xlsx")))

# Add a vector indicating if the sample is good, based on the colors of the 
# excel file (as read_xlsx can't see the colors)
metDat$isGood <- rep(T, nrow(metDat))
metDat$isGood[c(4, 11, 18)] <- F


# Uncompress files 
for(f in wgbsFiles){
        #f <- wgbsFiles[1]
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
basename(metDat$WGBS_Cov)[metDat$isGood]

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

# Analysis of common CpGs accross all the samples
################################################################################

# Load Bismark coverage files
wgbsFilesUncomp_list <- as.list(sprintf("%s%s", wgbsDir, wgbsFilesUncomp))
sampList <- as.list(sampInfo$Sample[match(wgbsFilesUncomp, sampInfo$WGBS_Cov)])
bismCovListRaw <- methRead(wgbsFilesUncomp_list,
                           sample.id = sampList,
                           #assembly = as.list(wgbsFilesUncomp),
                           assembly = "hg18",
                           pipeline = "bismarkCoverage",
                           context = "CpG",
                           mincov = 10,
                           treatment = as.numeric(!grepl("M", unlist(sampList))))

bismCov_united <- unite(bismCovListRaw)

getData(bismCov_united)

getSampleID(bismCov_united)



# Do PCA
commCpGs_percMat <- t(percMethylation(bismCov_united))

colnames(commCpGs_percMat) <- paste(paste0("chr", as.character(getData(bismCov_united)$chr)),
                                    as.character(getData(bismCov_united)$start),
                                    sep = "_")

commCpGs_percMat_std <- commCpGs_percMat[, apply(commCpGs_percMat, 2, sd) != 0]

commCpGs_percMat_std <- apply(commCpGs_percMat_std,
                              2,
                              function(x) (x - mean(x))/sd(x))

commCpGs_percMat_pca <- prcomp(commCpGs_percMat_std, scale. = F, center = F)

# Plot multiPC scoreplot
doPCAMultiPlot(commCpGs_percMat_pca, nComps = 5)
ggsave(filename = sprintf("%sCpG_bpRes_PCA.pdf", plotDir))





DF <- commCpGs_percMat
filt <- "CD19"
filtAcross <- "CAR"
compAcross <- "Time"
comp <- c(8, 38)
filtPVals <- .05

if(!is.null(filtAcross)){
        DF <- DF[sampInfo$Sample[sampInfo[, filtAcross] == filt], ]
}


pb <- txtProgressBar(min = 0,
                     max = ncol(DF),
                     #max = 80,
                     style = 3)
pVals <- c()
for(j in seq_along(colnames(DF))){
        j <- 1
        setTxtProgressBar(pb, j)
        pos <- colnames(DF)[j]
        #DF[, pos]
        compSamps <- list()
        print(j)
        print(comp)
        for(i in 1:length(comp)){
                cmp <- comp[i]
                print(cmp)
                smpsInGrp <- rownames(DF)[rownames(DF) %in% sampInfo$Sample[sampInfo[, compAcross] == cmp]]
                compSamps[[i]] <- smpsInGrp
        }
        wTst <- wilcox.test(DF[compSamps[[1]], pos], DF[compSamps[[2]], pos])
        pVal <- wTst$p.value
        pVals <- c(pVals, pVal)
}
close(pb)
resuDF <- data.frame(position = colnames(DF),
                     p.value = pVals, 
                     p.adj = p.adjust(pVals, method = "BH"))

if(!is.null(filtPVals)){
        resuDF <- resuDF[!is.nan(resuDF$p.value), ]
        resuDF <- resuDF[resuDF$p.adj <= filtPVals, ]
}


# Does mann whitney test across the groups defined by the user (either two
# CAR types, or two time points). Allows to prefilter according to a timepoint
# or a CAR type
doUnivarComps <- function(DF,
                          filt,
                          filtAcross,
                          comp,
                          compAcross,
                          filtPVals = NULL){
        #DF <- commCpGs_percMat
        #filt <- "HA"
        #filtAcross <- "CAR"
        #compAcross <- "Time"
        #comp <- c(8, 38)
        #filtPVals <- .05
        
        if(!is.null(filtAcross)){
                DF <- DF[sampInfo$Sample[sampInfo[, filtAcross] == filt], ]
        }
        
        
        pb <- txtProgressBar(min = 0,
                             max = ncol(DF),
                             #max = 80,
                             style = 3)
        pVals <- c()
        for(j in seq_along(colnames(DF))){
                setTxtProgressBar(pb, j)
                pos <- colnames(DF)[j]
                compSamps <- list()
                for(i in 1:length(comp)){
                        cmp <- comp[i]
                        smpsInGrp <- rownames(DF)[rownames(DF) %in% sampInfo$Sample[sampInfo[, compAcross] == cmp]]
                        compSamps[[i]] <- smpsInGrp
                }
                emptyBool <- unlist(lapply(compSamps,
                                           function(x) length(x) == 0))
                if(any(emptyBool)){
                        emptyComp <- comp[emptyBool]
                        print("\n")
                        print(sprintf("There are not enough samples of %s = %s",
                              compAcross, emptyComp))
                        stop()
                }
                wTst <- wilcox.test(DF[compSamps[[1]], pos], DF[compSamps[[2]], pos])
                pVal <- wTst$p.value
                pVals <- c(pVals, pVal)
        }
        close(pb)
        resuDF <- data.frame(position = colnames(DF),
                             p.value = pVals, 
                             p.adj = p.adjust(pVals, method = "BH"))
        
        if(!is.null(filtPVals)){
                resuDF <- resuDF[!is.nan(resuDF$p.value), ]
                resuDF <- resuDF[resuDF$p.adj <= filtPVals, ]
        }
        return(resuDF)
}


HA_8_vs_38_bpRes_univSign <- doUnivarComps(commCpGs_percMat,
                                           filt = "HA",
                                           filtAcross = "CAR",
                                           comp = c(8, 38),
                                           compAcross = "Time",
                                           filtPVals = .05)
# Nothing significant

CD19_8_vs_38_bpRes_univSign <- doUnivarComps(commCpGs_percMat,
                                             filt = "CD19",
                                             filtAcross = "CAR",
                                             comp = c(8, 38),
                                             compAcross = "Time",
                                             filtPVals = .05)

Mock_8_vs_38_bpRes_univSign <- doUnivarComps(commCpGs_percMat,
                                             filt = "Mock",
                                             filtAcross = "CAR",
                                             comp = c(8, 38),
                                             compAcross = "Time",
                                             filtPVals = .05)


getMethylationStats(bismCovListRaw[[1]], plot=TRUE,both.strands=FALSE)

# Do a tile analysis
nCores <- detectCores() - 2
tiles100bp <- tileMethylCounts(bismCovListRaw, win.size = 100, step.size = 100, cov.bases = 10,
                               mc.cores = nCores)


bismList <- list()
for(f in wgbsFilesUncomp){
        samp <- sampInfo$Sample[match(f, sampInfo$WGBS_Cov)]
        bismCov <- methRead(sprintf("%s%s", wgbsDir, f2),
                            sample.id = samp,
                            assembly = f,
                            pipeline = "bismarkCoverage")
        bismList[[samp]] <- bismCov
}
getCoverageStats(bismList)
