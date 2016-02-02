#!/bin/bash
#source clean.sh

cp config/parametersHighpt.h parameters.h 

DOFONLL=0
DOTRGCOMBINATION=0
DOFEEDDOWN=0
DOFITSPP=0
DOFITSPbPb=0
DOMCstudyPP=0
DOMCstudyPbPb=0
DOCrossSectionPP=0
DOCrossSectionPbPb=0
DORAA=1

FONLLDATINPUT="pp_d0meson_5TeV_y1"
FONLLDATINPUTBtoD="pp_Btod0meson_5TeV_y1"
FONLLDATINPUTB="pp_Bmeson_5TeV_y1"

FONLLOUTPUTFILE="output_pp_d0meson_5TeV_y1.root"
FONLLOUTPUTFILEBtoD="output_pp_Btod0meson_5TeV_y1.root"
FONLLOUTPUTFILEInclusiveD="output_inclusiveDd0meson_5TeV_y1.root"
FONLLOUTPUTFILEB="output_pp_Bmeson_5TeV_y1.root"
OUTPUTFILEPP="hPtSpectrumDzeroPP.root"
OUTPUTFILEPlotPP="CrossSectionFONLLPP.root"
OUTPUTFILEPbPb="hPtSpectrumDzeroPbPb.root"
OUTPUTFILEPlotPbPb="CrossSectionFONLLPbPb.root"

NTUPLAPYTHIA="/data/HeavyFlavourRun2/BtoDPythia/treefile_ptall_11january2016.root"
INPUTMCPP="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi/ntD_EvtBase_20160125_Dfinder_20151229_pp_Pythia8_prompt_D0_withweight.root"
INPUTDATAPP="/data/dmeson2015/DataDntuple/nt_20160112_DfinderData_pp_20160111_dPt0tkPt1_D0Dstar3p5p_DCSJSON_v2.root"
INPUTMCPbPb="/data/wangj/MC2015/Dntuple/PbPb/ntD_PbPb_Dzero_kpi/ntD_EvtBase_20160127_Dfinder_20151229_Pythia8_prompt_D0_noweight.root"
INPUTDATAPbPb="/data/dmeson2015/DataDntuple/nt_skimmed_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6_partialstats_v2.root"

LUMIPP=26.31
ISMCPP=0
ISDOWEIGHTPP=1
SELGENPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
SELGENPPACCPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))&&abs(Gtk1eta)<2.0&&abs(Gtk2eta)<2.0&&Gtk1pt>2.0&&Gtk2pt>2.0"
RECOONLYPP="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11"
CUTPP="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&(DsvpvDistance/DsvpvDisErr)>3.5&&(DlxyBS/DlxyBSErr)>1.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11"
TRGPP="((HLT_DmesonPPTrackingGlobal_Dpt15_v1&&Dpt>20&&Dpt<40)||(HLT_DmesonPPTrackingGlobal_Dpt30_v1&&Dpt>40&&Dpt<60)||(HLT_DmesonPPTrackingGlobal_Dpt50_v1&&Dpt>60))"
LABELPP="PP"
USEPRESCALEPP=1
                                         
LUMIPbPb=0.000404 #404*0.001*0.001
ISMCPbPb=0
ISDOWEIGHTPbPb=1
SELGENPbPb="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
SELGENPPACCPbPb="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))&&abs(Gtk1eta)<2.0&&abs(Gtk2eta)<2.0&&Gtk1pt>2.0&&Gtk2pt>2.0"
RECOONLYPbPb="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>8.5&&Dtrk2Pt>8.5&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&(Dtrk1PixelHit+Dtrk1StripHit)>=11"
CUTPbPb="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>8.5&&Dtrk2Pt>8.5&&(DsvpvDistance/DsvpvDisErr)>3.5&&(DlxyBS/DlxyBSErr)>2.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&(Dtrk1PixelHit+Dtrk1StripHit)>=11"
TRGPbPb="((HLT_HIDmesonHITrackingGlobal_Dpt20_v1&&Dpt>20&&Dpt<45)||(HLT_HIDmesonHITrackingGlobal_Dpt40_v1&&Dpt>45&&Dpt<70)||(HLT_HIDmesonHITrackingGlobal_Dpt60_v1&&Dpt>70))"
LABELPbPb="PbPb"
USEPRESCALEPbPb=1

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

if [ $DOFITSPbPb -eq 1 ]; then      
g++ fitD.C $(root-config --cflags --libs) -g -o fitD.exe 
./fitD.exe "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   "$LUMIPbPb"   "$ISDOWEIGHTPbPb"  "$LABELPbPb"  "$OUTPUTFILEPbPb"
fi

if [ $DOMCstudyPP -eq 1 ]; then      
g++ MCefficiency.C $(root-config --cflags --libs) -g -o MCefficiency.exe 
./MCefficiency.exe "$INPUTMCPP"  "$SELGENPP" "$SELGENPPACCPP"  "$RECOONLYPP" "$CUTPP"  "$LABELPP"
fi

if [ $DOMCstudyPbPb -eq 1 ]; then      
g++ MCefficiency.C $(root-config --cflags --libs) -g -o MCefficiency.exe 
./MCefficiency.exe "$INPUTMCPbPb"  "$SELGENPbPb" "$SELGENPPACCPbPb"  "$RECOONLYPbPb" "$CUTPbPb"  "$LABELPbPb"
fi

if [ $DOCrossSectionPP -eq 1 ]; then      
g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILEInclusiveD"  "$OUTPUTFILEPP" "$OUTPUTPrescalePP" "$USEPRESCALEPP" "$OUTPUTFILEPlotPP" "$LABELPP"
fi

if [ $DOCrossSectionPbPb -eq 1 ]; then      
g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILEInclusiveD"  "$OUTPUTFILEPbPb" "$OUTPUTPrescalePbPb" "$USEPRESCALEPbPb" "$OUTPUTFILEPlotPbPb" "$LABELPbPb"
fi

if [ $DORAA -eq 1 ]; then      
g++ NuclearModificationFactor.C $(root-config --cflags --libs) -g -o NuclearModificationFactor.exe 
./NuclearModificationFactor.exe "$OUTPUTFILEPlotPP" "$OUTPUTFILEPlotPbPb" 
fi


