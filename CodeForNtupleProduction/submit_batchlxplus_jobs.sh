####REMEMBER TO MOVE TO PP OR PbPb MODE
DIRECTORYOUTPUT="/afs/cern.ch/work/g/ginnocen/HeavyFlavour/DataHFAnalysis/CMSSW_7_5_5_patch4/src/DntupleRunII/CodeForNtupleProduction"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

rm  mylistfinal.txt

eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/694/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/694/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/695/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/695/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/703/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/703/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/726/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/726/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/735/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/AOD/000/262/735/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done

#cmscaf1nd
