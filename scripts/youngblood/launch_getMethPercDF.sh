#!/bin/bash

#conda activate r-4.3.1

bedDir="/home/users/gsantamaria/projects/car_t/data/youngblood_hg38Lift/"
regs="/home/users/gsantamaria/projects/car_t/results/model_build/100bpTiles/bestMod_regInfo.csv"
minCov=5
minNumCpGs=0
ref="hg38"
outName="/home/users/gsantamaria/projects/car_t/data/youngblood_hg38Lift/best100bpRegs_youngblood.csv"

Rscript getMethPercDF.R "$bedDir" --regs "$regs" --minCov "$minCov" --minNumCpGs "$minNumCpGs" --ref "$ref" --outName "$outName"

#conda deactivate