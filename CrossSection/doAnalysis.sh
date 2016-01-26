#!/bin/bash
#source clean.sh
DOFONLL=1
DOTRGCOMBINATION=1
DOFEEDDOWN=1
DOFITSPP=0
DOFITSPPMB=1
DOFITSPbPb=0
DOCrossSectionPP=0
DOCrossSectionPPMB=1
DORAA=0

FONLLDATINPUT="pp_d0meson_5TeV_y1"
FONLLDATINPUTBtoD="pp_Btod0meson_5TeV_y1"
FONLLDATINPUTB="pp_Bmeson_5TeV_y1"
FONLLOUTPUTFILE="output_pp_d0meson_5TeV_y1.root"
FONLLOUTPUTFILEBtoD="output_pp_Btod0meson_5TeV_y1.root"
FONLLOUTPUTFILEInclusiveD="output_inclusiveDd0meson_5TeV_y1.root"
FONLLOUTPUTFILEB="output_pp_Bmeson_5TeV_y1.root"
NTUPLAPYTHIA="/data/HeavyFlavourRun2/BtoDPythia/treefile_ptall_11january2016.root"

LUMIPP=26.31
INPUTMCPP="/afs/cern.ch/work/w/wangj/public/Dmeson/ntD_20151110_DfinderMC_20151110_EvtMatching_Pythia_D0pt15p0_Pthat15_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_MBseed_twang-Pythia_1107.root"
INPUTDATAPP="/data/dmeson2015/DataDntuple/nt_20160112_DfinderData_pp_20160111_dPt0tkPt1_D0Dstar3p5p_DCSJSON_v2.root"
ISMCPP=0
ISDOWEIGHTPP=1
CUTPP="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5"
SELGENPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPP="((HLT_DmesonPPTrackingGlobal_Dpt15_v1&&Dpt>15&&Dpt<40)||(HLT_DmesonPPTrackingGlobal_Dpt30_v1&&Dpt>40&&Dpt<60)||(HLT_DmesonPPTrackingGlobal_Dpt50_v1&&Dpt>60))"
LABELPP="PP"
OUTPUTFILEPP="hPtSpectrumDzeroPP.root"
USEPRESCALEPP=1
                                         
LUMIPPMB=0.0047         #assuming sigma=70mb, Nevents=367 millions, estimated MB efficiency=0.90  Lumi=0.367*10^9/(70*10^9pb)=0.367/70*0.9=0.0047
INPUTMCPPMB="/afs/cern.ch/work/w/wangj/public/Dmeson/ntD_20150924_MC_merge_withoutweight_skim.root"
INPUTDATAPPMB="/data/yjlee/dmeson/2015/trigger/mb.root"
ISMCPPMB=0
ISDOWEIGHTPPMB=1
CUTPPMB="Dy>-1.&&Dy<1.&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5"
SELGENPPMB="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPPMB="(HLT_L1MinimumBiasHF1OR_part1_v1||HLT_L1MinimumBiasHF1OR_part2_v1||HLT_L1MinimumBiasHF1OR_part3_v1||HLT_L1MinimumBiasHF1OR_part4_v1||HLT_L1MinimumBiasHF1OR_part5_v1||HLT_L1MinimumBiasHF1OR_part6_v1||HLT_L1MinimumBiasHF1OR_part7_v1||HLT_L1MinimumBiasHF1OR_part8_v1||HLT_L1MinimumBiasHF1OR_part9_v1||HLT_L1MinimumBiasHF1OR_part10_v1||HLT_L1MinimumBiasHF1OR_part11_v1||HLT_L1MinimumBiasHF1OR_part12_v1||HLT_L1MinimumBiasHF1OR_part13_v1||HLT_L1MinimumBiasHF1OR_part14_v1||HLT_L1MinimumBiasHF1OR_part15_v1||HLT_L1MinimumBiasHF1OR_part16_v1||HLT_L1MinimumBiasHF1OR_part17_v1||HLT_L1MinimumBiasHF1OR_part18_v1||HLT_L1MinimumBiasHF1OR_part19_v1)"
LABELPPMB="PPMB"
OUTPUTFILEPPMB="hPtSpectrumDzeroPPMB.root"
USEPRESCALEPPMB=0


INPUTMCPbPb="/afs/cern.ch/work/w/wangj/public/Dmeson/ntD_20151110_DfinderMC_20151110_EvtMatching_Pythia_D0pt15p0_Pthat15_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_MBseed_twang-Pythia_1107.root"
INPUTDATAPbPb="/data/HeavyFlavourRun2/crab_DfinderData_PbPb_20151227_dPt10tkPt2p5_D0Dstar3p5p_CameliaJSON/merged_ntuple.root"
LUMIPbPb=0.000014930
ISMCPbPb=0
ISDOWEIGHTPbPb=1
CUTPbPb="Dy>-1.&&Dy<1.&&(Dtrk1Pt>8.5&&Dtrk2Pt>8.5)&&((Dpt<20&&(DsvpvDistance/DsvpvDisErr)>4.50&&Dchi2cl>0.05&&Dalpha<0.12)||(Dpt>20&&Dpt<40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.2&&Dalpha<0.12)||(Dpt>40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.05&&Dalpha<0.12))"
SELGENPbPb="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPbPb="((HLT_HIDmesonHITrackingGlobal_Dpt20_v1&&Dpt>20&&Dpt<50)||(HLT_HIDmesonHITrackingGlobal_Dpt40_v1&&Dpt>50&&Dpt<70)||(HLT_HIDmesonHITrackingGlobal_Dpt60_v1&&Dpt>70))"
LABELPbPb="PbPb"
OUTPUTFILEPbPb="hPtSpectrumDzeroPbPb.root"

OUTPUTPrescalePP="prescalePP.root"
OUTPUTPrescalePbPb="prescalePbPb.root"

if [ $DOFONLL -eq 1 ]; then      

g++ Dzerodsigmadpt.cc $(root-config --cflags --libs) -g -o Dzerodsigmadpt.exe 
./Dzerodsigmadpt.exe "$FONLLDATINPUT"  "$FONLLOUTPUTFILE" 
./Dzerodsigmadpt.exe "$FONLLDATINPUTBtoD"  "$FONLLOUTPUTFILEBtoD" 

g++ BplusAlldsigmadpt.cc $(root-config --cflags --libs) -g -o BplusAlldsigmadpt.exe 
./BplusAlldsigmadpt.exe "$FONLLDATINPUTB"  "$FONLLOUTPUTFILEB" 
fi 

if [ $DOFEEDDOWN -eq 1 ]; then      

g++ RatioFeedDown.cc $(root-config --cflags --libs) -g -o RatioFeedDown.exe 
./RatioFeedDown.exe "$FONLLOUTPUTFILE"  "$FONLLOUTPUTFILEBtoD" "$FONLLOUTPUTFILEInclusiveD"

g++ plotFeedDown.C $(root-config --cflags --libs) -g -o plotFeedDown.exe 
./plotFeedDown.exe "$FONLLOUTPUTFILE"  "$FONLLOUTPUTFILEB" "$NTUPLAPYTHIA" 1 0 100 

fi

if [ $DOTRGCOMBINATION -eq 1 ]; then      

g++ triggercombination.cc $(root-config --cflags --libs) -g -o triggercombination.exe 
./triggercombination.exe "$LABELPP"  "$INPUTDATAPP" "$OUTPUTPrescalePP"

g++ triggercombination.cc $(root-config --cflags --libs) -g -o triggercombination.exe 
./triggercombination.exe "$LABELPbPb"  "$INPUTDATAPbPb" "$OUTPUTPrescalePbPb"
fi

if [ $DOFITSPP -eq 1 ]; then      
g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPP"
fi

if [ $DOFITSPPMB -eq 1 ]; then      
g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPPMB"  "$INPUTMCPPMB"  "$TRGPPMB" "$CUTPPMB"   "$SELGENPPMB"   "$ISMCPPMB"   "$LUMIPPMB"   "$ISDOWEIGHTPPMB"   "$LABELPPMB"  "$OUTPUTFILEPPMB"
fi


if [ $DOFITSPbPb -eq 1 ]; then      
g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   "$LUMIPbPb"   "$ISDOWEIGHTPbPb"  "$LABELPbPb"  "$OUTPUTFILEPbPb"
fi

if [ $DOCrossSectionPP -eq 1 ]; then      
g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILEInclusiveD"  "$OUTPUTFILEPP" "$OUTPUTPrescalePP" "$USEPRESCALEPP"
fi

if [ $DOCrossSectionPPMB -eq 1 ]; then      
g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILEInclusiveD"  "$OUTPUTFILEPPMB" "$OUTPUTPrescalePP" "$USEPRESCALEPPMB"
fi


if [ $DORAA -eq 1 ]; then      
g++ NuclearModification.C $(root-config --cflags --libs) -g -o NuclearModification.exe 
./NuclearModification.exe "$OUTPUTFILEPP" "$OUTPUTFILEPbPb" "$OUTPUTPrescalePP" "$OUTPUTPrescalePbPb"
fi


