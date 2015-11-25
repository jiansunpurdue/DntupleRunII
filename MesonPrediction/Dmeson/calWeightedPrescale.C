#include "../Bmeson/Bplusdsigmadpt_1ptbins.h"
#include "iostream"
void calWeightedPrescale(){
    string fonlldata="../fonllData/D_pp_pt_rap24_5p5_pt0to200.dat";
    double diffXsec1 = 0;
    double diffXsec2 = 0;
    double diffXsec3 = 0;
    Bplusdsigmadpt_1ptbins(false, diffXsec1, 7, 8, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec2, 8, 9, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec3, 7, 9, fonlldata);
	diffXsec1=diffXsec1*1;
	diffXsec2=diffXsec2*1;
	diffXsec3=diffXsec3*2;
	cout<<diffXsec3/((diffXsec1/750)+(diffXsec2/140))<<endl;

    Bplusdsigmadpt_1ptbins(false, diffXsec1, 13, 15, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec2, 15, 16, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec3, 13, 16, fonlldata);
	diffXsec1=diffXsec1*2;
	diffXsec2=diffXsec2*1;
	diffXsec3=diffXsec3*3;
	cout<<diffXsec3/((diffXsec1/140)+(diffXsec2/33))<<endl;

    Bplusdsigmadpt_1ptbins(false, diffXsec1, 28, 30, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec2, 30, 40, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec3, 28, 40, fonlldata);
	diffXsec1=diffXsec1*2;
	diffXsec2=diffXsec2*10;
	diffXsec3=diffXsec3*12;
	cout<<diffXsec3/((diffXsec1/10)+(diffXsec2/5))<<endl;

    Bplusdsigmadpt_1ptbins(false, diffXsec1, 40, 50, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec2, 50, 60, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec3, 48, 60, fonlldata);
	diffXsec1=diffXsec1*10;
	diffXsec2=diffXsec2*10;
	diffXsec3=diffXsec3*20;
	cout<<diffXsec3/((diffXsec1/5)+(diffXsec2/1))<<endl;

    Bplusdsigmadpt_1ptbins(false, diffXsec1, 28, 30, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec2, 30, 40, fonlldata);
    Bplusdsigmadpt_1ptbins(false, diffXsec3, 28, 40, fonlldata);
	diffXsec1=diffXsec1*2;
	diffXsec2=diffXsec2*10;
	diffXsec3=diffXsec3*12;
	cout<<diffXsec3/((diffXsec1/50.)+(diffXsec2/34.46))<<endl;
}
