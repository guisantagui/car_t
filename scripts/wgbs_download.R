################################################################################
# CAR-T project: download the Bismark cov files                                #
################################################################################

rootDir <- "/home/users/gsantamaria/projects/car_t/"

dataDir <- paste0(rootDir, "data/")
wgbsDir <- paste0(dataDir, "wbgs/")

metaDF <- readxl::read_xlsx(paste0(dataDir,
                                   "Tcells2023WGBSMeta_Links_collab.xlsx"))

metaDF <- data.frame(metaDF)


for(link in metaDF$WGBS_Cov){
        #link <- metaDF$WGBS_Cov[2]
        cmd <- sprintf("wget %s -P %s", link, wgbsDir)
        fileNam <- basename(link)
        if(!file.exists(sprintf("%s%s", wgbsDir, fileNam))){
                system(cmd)
        }else{
                print(sprintf("%s is already in %s", fileNam, wgbsDir))
        }
}