source xeon_scripts/common-variables.sh forPR
source xeon_scripts/init-env.sh
make distclean -j 32 WITH_ROOT:=1 AVX_512:=1
make clean
make -j 32 WITH_ROOT:=1 AVX_512:=1
mkdir -p tmp

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-std --cmssw-val-mkfiteff
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVAL.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-ce --cmssw-val-mkfiteff
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVAL.root

############### OPTS ########################

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-std --cmssw-val-simmatch
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALSIM.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-ce --cmssw-val-simmatch
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALSIM.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-std --cmssw-val-simsignalmatch --cmssw-val-mkfiteffsimsignal
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALSIMSIG.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 1000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-ce --cmssw-val-simsignalmatch --cmssw-val-mkfiteffsimsignal
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALSIMSIG.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 10000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-std --cmssw-val-nomatch
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALNOSIM.root

./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 10000 --remove-dup --cmssw-val-fhit-bhit --read-cmssw-tracks --try-to-save-sim-info --backward-fit-pca --build-ce --cmssw-val-nomatch
hadd -O valtree.root valtree_*.root
rm valtree_*.root
mv valtree.root tmp/valtree_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALNOSIM.root

make clean -j 32 WITH_ROOT:=1 AVX_512:=1
mv tmp/valtree_*.root .
rm -rf tmp

root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVAL",1,1,0,"pdf",1)'
root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVAL",1,1,0,"pdf",1)'
root -b -q -l 'plotting/makeValidation.C("SKL-SP_CMSSW_TTbar_PU50","_CMSSWVAL",1,"forPR",1)'

root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALSIM",1)'
root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALSIM",1)'
root -b -q -l 'plotting/makeValidation.C("SKL-SP_CMSSW_TTbar_PU50","_CMSSWVALSIM",1,"forPR")'

root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALSIMSIG",1,1,0,"pdf",1)'
root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALSIMSIG",1,1,0,"pdf",1)'
root -b -q -l 'plotting/makeValidation.C("SKL-SP_CMSSW_TTbar_PU50","_CMSSWVALSIMSIG",1,"forPR",1)'

root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_STD_CMSSWVALNOSIM",1)'
root -b -q -l 'plotting/runValidation.C("_SKL-SP_CMSSW_TTbar_PU50_CE_CMSSWVALNOSIM",1)'
root -b -q -l 'plotting/makeValidation.C("SKL-SP_CMSSW_TTbar_PU50","_CMSSWVALNOSIM",1,"forPR")'

make distclean -j 32 WITH_ROOT:=1 AVX_512:=1


# # # 	"  --cmssw-val-fprm-bprm    use CMSSW validation with track parameter based matching for forward built tracks\n"

