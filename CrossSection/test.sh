FONLLOUTPUTFILE="output_pp_d0meson5_5TeV_y1.root"
OUTPUTFILEPP="hPtSpectrumDzeroPP.root"
OUTPUTFILEPbPb="hPtSpectrumDzeroPbPb.root"

g++ CrossSectionRatio.C $(root-config --cflags --libs) -g -o CrossSectionRatio.exe 
./CrossSectionRatio.exe "$FONLLOUTPUTFILE"  "$OUTPUTFILEPP"  "$OUTPUTFILEPbPb"
