#! /bin/bash

# command line input
dir=${1:-"testcmsswval"} # Main output dir name
user=${2:-"legianni"} # lxplus user name

mkdir -p ${dir}

# Make SimTrack Validation directories
simdir=("CMSSWVAL" "CMSSWVALSIM" "CMSSWVALSIMSIG" "CMSSWVALNOSIM")

for((i=0;i<${#simdir[@]};++i));do

mkdir -p ${dir}/${simdir[i]}
mkdir -p ${dir}/${simdir[i]}/eff
mkdir -p ${dir}/${simdir[i]}/ineff
mkdir -p ${dir}/${simdir[i]}/dr
mkdir -p ${dir}/${simdir[i]}/mkfiteff
mkdir -p ${dir}/${simdir[i]}/diffs
mkdir -p ${dir}/${simdir[i]}/nHits
mkdir -p ${dir}/${simdir[i]}/score

mv *mkfiteff*${simdir[i]}.png ${dir}/${simdir[i]}/mkfiteff
mv *ineff*${simdir[i]}.png ${dir}/${simdir[i]}/ineff
mv *eff*${simdir[i]}.png ${dir}/${simdir[i]}/eff
mv *dr*${simdir[i]}.png ${dir}/${simdir[i]}/dr
mv *_d*${simdir[i]}.png ${dir}/${simdir[i]}/diffs
mv *_nHits_*${simdir[i]}.png ${dir}/${simdir[i]}/nHits
mv *_score_*${simdir[i]}.png ${dir}/${simdir[i]}/score


mv validation_SKL-SP_CMSSW_TTbar_PU50_STD_${simdir[i]}/*txt ${dir}/${simdir[i]}
mv validation_SKL-SP_CMSSW_TTbar_PU50_STD_${simdir[i]}/*.txt ${dir}/${simdir[i]}

done

#cp index.php into all subdirectories
find ${dir} -mindepth 0 -type d -exec cp web/index.php {} \; 

# Then copy to lxplus
echo "Moving plots and text files remotely to lxplus"
./web/tarAndSendToLXPLUS.sh ${dir} forPR afs ${user}

# Final cleanup of directory
echo "Removing local files"
./xeon_scripts/trashSKL-SP.sh 0 
rm -rf ${dir}

# Final message
echo "Finished moving benchmark plots to LXPLUS!"
