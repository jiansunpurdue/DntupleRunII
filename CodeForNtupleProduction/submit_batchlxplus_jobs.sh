####REMEMBER TO MOVE TO PP OR PbPb MODE
DIRECTORYOUTPUT="./"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

rm mylist.txt mylistfinal.txt
eos ls /store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/Merged >>mylist.txt
awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIHardProbes/Merged/" $0}' mylist.txt >>mylistfinal.txt
rm mylist.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q 8nh $i ; done

