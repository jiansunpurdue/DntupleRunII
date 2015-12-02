#DIRECTORYOUTPUT="/afs/cern.ch/work/g/ginnocen/DmesonNtuple/pp/nt_HiForest_MinimumBias1_Run2015E_PromptReco_pp_5020GeV"
DIRECTORYOUTPUT="./"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

count=0 ; for i in `cat /afs/cern.ch/user/v/velicanu/public/forgm/lxplusdoga.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q 1nh $i ; done

#bjobs
#cat doga.txt | sed 's@/mnt/hadoop/cms@root://xrootd.cmsaf.mit.edu/@g' > lxplusdoga.txt

