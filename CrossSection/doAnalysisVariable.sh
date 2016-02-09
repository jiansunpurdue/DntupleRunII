#!/bin/bash
#source clean.sh

DOFITSPPVariable=1

cp config/parametersVariables.h parameters.h 

INPUTMCPP="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi/ntD_EvtBase_20160203_Dfinder_20160201_pp_Pythia8_prompt_D0_dPt0tkPt0p5_pthatweight.root"
INPUTMCNPPP="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi_nonprompt/ntD_EvtBase_20160203_Dfinder_20160201_pp_Pythia8_nonprompt_D0_dPt0tkPt0p5_pthatweight.root"
#INPUTMCPP="/data/wangj/MC2015/Dntuple/backup/ntD_EvtBase_20160125_Dfinder_20151229_pp_Pythia8_prompt_D0_pthatweight.root"
#INPUTMCNPPP="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi_nonprompt/ntD_EvtBase_20160203_Dfinder_20160201_pp_Pythia8_nonprompt_D0_dPt0tkPt0p5_pthatweight.root"
INPUTDATAPP="/data/dmeson2015/DataDntuple/nt_20160112_DfinderData_pp_20160111_dPt0tkPt1_D0Dstar3p5p_DCSJSON_v2.root"

VARIABLE="DsvpvDistance/DsvpvDisErr";
VARIABLEPLOT="d_{xy}/#sigma(d_{xy})"
#VARIABLE="Dpt";
#VARIABLEPLOT="Dpt"
LABELPP="MCPrompt"
LUMIPP=1
ISMCPP=1
ISDOWEIGHTPP=1
SELGENPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
SELGENPPACCPP="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))&&abs(Gtk1eta)<2.0&&abs(Gtk2eta)<2.0&&Gtk1pt>2.0&&Gtk2pt>2.0"
CUTPP="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&(DsvpvDistance/DsvpvDisErr)>3&&(DlxyBS/DlxyBSErr)>1.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11"
TRGPP="(1)"
OUTPUTFILEPPVariable="MCCutVariablePrompt.root"

g++ fitDVariable.C $(root-config --cflags --libs) -g -o fitDVariable.exe 
./fitDVariable.exe "$VARIABLE" "$VARIABLEPLOT" "$INPUTMCPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPPVariable"

#OUTPUTFILEPPVariable="MCCutVariableNonPrompt.root"
#LABELPP="MCNonPrompt"
#./fitDVariable.exe "$VARIABLE" "$VARIABLEPLOT" "$INPUTMCNPPP"  "$INPUTMCNPPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPPVariable"

#ISMCPP=0
#TRGPP="HLT_DmesonPPTrackingGlobal_Dpt15_v1"
#OUTPUTFILEPPVariable="DataCutVariable.root"
#LABELPP="Data"
#./fitDVariable.exe "$VARIABLE" "$VARIABLEPLOT" "$INPUTDATAPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   "$LUMIPP"   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPPVariable"
