#!/bin/bash
source clean.sh
FONLLDATINPUT="pp_d0meson5_5TeV_y1"
FONLLOUTPUTFILE="output_pp_d0meson5_5TeV_y1.root"

INPUTMCPP="/data/dmeson2015/MCDntuple/ntD_20151115_DfinderMC_20151110_EvtMatching_Pythia_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_7415_v20_1116_Pthat5_15_35merged.root"
INPUTDATAPP="/data/HeavyFlavourRun2/DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/merged_ntuple.root"
LUMIPP=13.9
ISMCPP=0
ISDOWEIGHTPP=1
CUTPP="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5"
SELGENPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPP="((HLT_DmesonPPTrackingGlobal_Dpt15_v1&&Dpt>15&&Dpt<50)||(HLT_DmesonPPTrackingGlobal_Dpt50_v1&&Dpt>50))"
LABELPP="PP"
OUTPUTFILEPP="hPtSpectrumDzeroPP.root"

INPUTMCPbPb="/data/dmeson2015/MCDntuple/ntD_20151115_DfinderMC_20151110_EvtMatching_Pythia_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_7415_v20_1116_Pthat5_15_35merged.root"
INPUTDATAPbPb="/afs/cern.ch/work/g/ginnocen/HeavyFlavourRun2/crab_DfinderData_PbPb_20151227_dPt10tkPt2p5_D0Dstar3p5p_CameliaJSON/merged_ntuple.root"
LUMIPbPb=0.000014930
ISMCPbPb=0
ISDOWEIGHTPbPb=1
CUTPbPb="Dy>-1.&&Dy<1.&&(Dtrk1Pt>8.5&&Dtrk2Pt>8.5)&&((Dpt<20&&(DsvpvDistance/DsvpvDisErr)>4.50&&Dchi2cl>0.05&&Dalpha<0.12)||(Dpt>20&&Dpt<40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.2&&Dalpha<0.12)||(Dpt>40&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dchi2cl>0.05&&Dalpha<0.12))"
SELGENPbPb="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
TRGPbPb="(HLT_HIDmesonHITrackingGlobal_Dpt60_v1)"
LABELPbPb="PbPb"
OUTPUTFILEPbPb="hPtSpectrumDzeroPbPb.root"

g++ Dzerodsigmadpt5TeV.cc $(root-config --cflags --libs) -g -o Dzerodsigmadpt5TeV.exe 
./Dzerodsigmadpt5TeV.exe "$FONLLDATINPUT"  "$FONLLOUTPUTFILE" 

g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPP"

g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   "$LUMIPbPb"   "$ISDOWEIGHTPbPb"  "$LABELPbPb"  "$OUTPUTFILEPbPb"

g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILE"  "$OUTPUTFILEPP"  "$OUTPUTFILEPbPb"
