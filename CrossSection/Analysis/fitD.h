TString cmcut = "Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)";
TString ptcut0 = "(DsvpvDistance/DsvpvDisErr)>3.50&&Dchi2cl>0.05&&Dalpha<0.12";
TString cut0 = Form("%s&&%s",cmcut.Data(),ptcut0.Data());
TString selmcgen = Form("(GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1)");

//TString ptcut1 = "Dalpha<0.12&&((Dpt>1.0&&Dpt<3.5&&(DsvpvDistance/DsvpvDisErr)>5.90&&Dchi2cl>0.248) || (Dpt>3.5&&Dpt<4.5&&(DsvpvDistance/DsvpvDisErr)>5.81&&Dchi2cl>0.200) || (Dpt>4.5&&Dpt<5.5&&(DsvpvDistance/DsvpvDisErr)>5.10&&Dchi2cl>0.191) || (Dpt>5.5&&Dpt<7.0&&(DsvpvDistance/DsvpvDisErr)>4.62&&Dchi2cl>0.148) || (Dpt>7.0&&Dpt<9.0&&(DsvpvDistance/DsvpvDisErr)>4.46&&Dchi2cl>0.102) || (Dpt>9.0&&Dpt<11.&&(DsvpvDistance/DsvpvDisErr)>4.39&&Dchi2cl>0.080) || (Dpt>11.&&Dpt<13.&&(DsvpvDistance/DsvpvDisErr)>4.07&&Dchi2cl>0.073) || (Dpt>13.&&Dpt<16.&&(DsvpvDistance/DsvpvDisErr)>3.88&&Dchi2cl>0.060) || (Dpt>16.&&Dpt<20.&&(DsvpvDistance/DsvpvDisErr)>3.67&&Dchi2cl>0.055) || (Dpt>20.&&Dpt<28.&&(DsvpvDistance/DsvpvDisErr)>3.25&&Dchi2cl>0.054) || (Dpt>28.&&Dpt<40.&&(DsvpvDistance/DsvpvDisErr)>2.55&&Dchi2cl>0.043))";
