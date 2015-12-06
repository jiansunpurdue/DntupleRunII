DIRECTORYINPUT="/afs/cern.ch/work/g/ginnocen/HeavyFlavour/DataHFAnalysis/CMSSW_7_5_5_patch4/src/DntupleRunII/CodeForNtupleProduction"
DIRECTORYOUTPUT="/data/dmeson2015/PbPbNtuple/ntuple_HIRun2015HIHardProbesAOD12062015_test"
cd $DIRECTORYINPUT
cmsenv
mkdir $DIRECTORYOUTPUT
cd  $DIRECTORYOUTPUT
cp $DIRECTORYINPUT/ntuple_*.root .
hadd ntuple_merged.root ntuple_*.root