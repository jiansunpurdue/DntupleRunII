# Instructions for pp references
In this repository you can find the code to build the pp reference of the B meson measurements in pPb collisions at 5.02 TeV starting from the FONLL calculations. The input files .dat were produced using http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html using the following parameters for the calculations:

Collider: LHC(pp, 5.5), Heavy quark: botton ,PDFs: CTEQ6.6

Perturbative order: FONLL, Final state: D0 hadron, Further decay: -

Cross section type: dsigma/dpt (or dsigma/dy) uncertainty range from scales and masses

Include PDFs uncertainties: (missing for the moment)

The file fo_pp_d0meson5_5TeV.dat was produced with the kinematic ranges:

pt min = 2 GeV
pt max = 60 GeV
y min= -1
y max= 1
Use y:check
npoints=393
