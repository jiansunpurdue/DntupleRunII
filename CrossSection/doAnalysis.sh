#!/bin/bash
source clean.sh
DOFITS=1
FONLLDATINPUT="pp_d0meson_5TeV_y1"
FONLLDATINPUTBtoD="pp_Btod0meson_5TeV_y1"
FONLLDATINPUTB="pp_Bmeson_5TeV_y1"

FONLLOUTPUTFILE="output_pp_d0meson_5TeV_y1.root"
FONLLOUTPUTFILEBtoD="output_pp_Btod0meson_5TeV_y1.root"
FONLLOUTPUTFILEInclusiveD="output_inclusiveDd0meson_5TeV_y1.root"
FONLLOUTPUTFILEB="output_pp_Bmeson_5TeV_y1.root"


INPUTMCPP="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi/ntD_EvtBase_20160107_Dfinder_20151229_pp_Pythia8_prompt_D0pt30p0_Pthat30_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20151212_dPt1tkPt1_D0Ds.root"
INPUTDATAPP="/data/HeavyFlavourRun2/DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/merged_ntuple.root"
LUMIPP=13.9
ISMCPP=0
ISDOWEIGHTPP=1
CUTPP="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5"
SELGENPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPP="((HLT_DmesonPPTrackingGlobal_Dpt15_v1&&Dpt>25&&Dpt<40)||(HLT_DmesonPPTrackingGlobal_Dpt30_v1&&Dpt>40&&Dpt<60)||(HLT_DmesonPPTrackingGlobal_Dpt50_v1&&Dpt>60))"
LABELPP="PP"
OUTPUTFILEPP="hPtSpectrumDzeroPP.root"


INPUTMCPbPb="/data/dmeson2015/MCDntuple/ntD_20151115_DfinderMC_20151110_EvtMatching_Pythia_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_7415_v20_1116_Pthat5_15_35merged.root"
INPUTDATAPbPb="/afs/cern.ch/work/g/ginnocen/HeavyFlavourRun2/crab_DfinderData_PbPb_20151227_dPt10tkPt2p5_D0Dstar3p5p_CameliaJSON/merged_ntuple.root"
LUMIPbPb=0.000014930
ISMCPbPb=0
ISDOWEIGHTPbPb=1
CUTPbPb="Dy>-1.&&Dy<1.&&(Dtrk1Pt>8.5&&Dtrk2Pt>8.5)&&((Dpt<20&&(DsvpvDistance/DsvpvDisErr)>4.50&&Dchi2cl>0.05&&Dalpha<0.12)||(Dpt>20&&Dpt<40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.2&&Dalpha<0.12)||(Dpt>40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.05&&Dalpha<0.12))"
SELGENPbPb="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPbPb="((HLT_HIDmesonHITrackingGlobal_Dpt20_v1&&Dpt>20&&Dpt<40)||(HLT_HIDmesonHITrackingGlobal_Dpt40_v1&&Dpt>40&&Dpt<60)||(HLT_HIDmesonHITrackingGlobal_Dpt60_v1&&Dpt>60))"
LABELPbPb="PbPb"
OUTPUTFILEPbPb="hPtSpectrumDzeroPbPb.root"

OUTPUTPrescalePP="prescalePP.root"
OUTPUTPrescalePbPb="prescalePbPb.root"

g++ Dzerodsigmadpt.cc $(root-config --cflags --libs) -g -o Dzerodsigmadpt.exe 
./Dzerodsigmadpt.exe "$FONLLDATINPUT"  "$FONLLOUTPUTFILE" 
./Dzerodsigmadpt.exe "$FONLLDATINPUTBtoD"  "$FONLLOUTPUTFILEBtoD" 

g++ BplusAlldsigmadpt.cc $(root-config --cflags --libs) -g -o BplusAlldsigmadpt.exe 
./BplusAlldsigmadpt.exe "$FONLLDATINPUTB"  "$FONLLOUTPUTFILEB" 

if [ $DOFITS -eq 1 ]; then      

#g++ RatioFeedDown.cc $(root-config --cflags --libs) -g -o RatioFeedDown.exe 
#./RatioFeedDown.exe "$FONLLOUTPUTFILE"  "$FONLLOUTPUTFILEBtoD" "$FONLLOUTPUTFILEInclusiveD"

g++ triggercombination.cc $(root-config --cflags --libs) -g -o triggercombination.exe 
./triggercombination.exe "$LABELPP"  "$INPUTDATAPP" "$OUTPUTPrescalePP"

g++ triggercombination.cc $(root-config --cflags --libs) -g -o triggercombination.exe 
./triggercombination.exe "$LABELPbPb"  "$INPUTDATAPbPb" "$OUTPUTPrescalePbPb"

g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPP"

g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   "$LUMIPbPb"   "$ISDOWEIGHTPbPb"  "$LABELPbPb"  "$OUTPUTFILEPbPb"

g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILE"  "$OUTPUTFILEPP"  "$OUTPUTFILEPbPb"  "$OUTPUTPrescalePP"  "$OUTPUTPrescalePbPb"
 
fi