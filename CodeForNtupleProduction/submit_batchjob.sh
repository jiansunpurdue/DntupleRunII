#!/bin/sh

#need to know total number of events in the file first
python submit_Ntuple_batch_cfg.py -i root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestExpress_run262656.root -n 10000 -j 12 -o /store/group/phys_heavyions/jisun/Dntuple2015run/PbPbExpress
