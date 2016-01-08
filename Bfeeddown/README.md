# README FILE TO GENERATE B->D PYTHIA

1) Make sure you have the latest PythiaMomDauFilter decayer.
    First git cms-addpkg GeneratorInterface/GenFilters. Then:
    -  go to GeneratorInterface/GenFilters/src 
        curl -O https://raw.githubusercontent.com/taweiXcms/cmssw/_PR_AddPythiaMomDauFilter_20151218/GeneratorInterface/GenFilters/src/PythiaMomDauFilter.cc
    - GeneratorInterface/GenFilters/interface
        curl -O https://raw.githubusercontent.com/taweiXcms/cmssw/_PR_AddPythiaMomDauFilter_20151218/GeneratorInterface/GenFilters/interface/PythiaMomDauFilter.h

2) Copy the non prompt python config files in your Configuration/Generator/python. 
    If you dont have it yet do:  git cms-addpkg Configuration/Generator
    cd Configuration/Generator/python
    and then curl all these files
    curl -O https://raw.githubusercontent.com/taweiXcms/HFAnaGenFrags/master/Run2Ana/DAnaPbPb/python/Pythia8_nonprompt_D0pt0p0_Pthat0_TuneCUETP8M1_5020GeV_cfi_evtgen130.py
    curl -O https://raw.githubusercontent.com/taweiXcms/HFAnaGenFrags/master/Run2Ana/DAnaPbPb/python/Pythia8_nonprompt_D0pt15p0_Pthat15_TuneCUETP8M1_5020GeV_cfi_evtgen130.py
    curl -O https://raw.githubusercontent.com/taweiXcms/HFAnaGenFrags/master/Run2Ana/DAnaPbPb/python/Pythia8_nonprompt_D0pt30p0_Pthat30_TuneCUETP8M1_5020GeV_cfi_evtgen130.py
    curl -O https://raw.githubusercontent.com/taweiXcms/HFAnaGenFrags/master/Run2Ana/DAnaPbPb/python/Pythia8_nonprompt_D0pt50p0_Pthat50_TuneCUETP8M1_5020GeV_cfi_evtgen130.py

3) Get the modified version of the HydjetAnalyzer in this folder and move it to GeneratorInterface/HydjetInterface/test.
    If you dont have it yet do git cms-addpkg GeneratorInterface/HydjetInterface
    then cp the modified version

4) Create the files with cmsDriver as below: 
    
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt0p0_Pthat0_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat0.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat0.py --no_exec -n 30000
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt15p0_Pthat15_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat15.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat15.py --no_exec -n 30000
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt30p0_Pthat30_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat30.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat30.py --no_exec -n 30000
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt50p0_Pthat50_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat50.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat50.py --no_exec -n 30000
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt70p0_Pthat70_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat70.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat70.py --no_exec -n 30000
    cmsDriver.py Configuration/Generator/python/Pythia8_nonprompt_D0pt90p0_Pthat90_TuneCUETP8M1_5020GeV_cfi_evtgen130.py --fileout file:pp_GENSIM_pthat90.root --eventcontent=RAWSIM --datatier GEN-SIM --step GEN,SIM --conditions auto:mc --processName GENERATOR python_filename pthat90.py --no_exec -n 30000


4b) Change the input and output numbers 

     input = cms.untracked.int32(-1),
     output  = cms.untracked.int32(1000)

It will run until it gets to 1000 events.

5) Then replace/attach these lines at the end of the file to be able to run the Hydjetanalyzer at the same time:


process.ProductionFilterSequence = cms.Sequence(process.generator+process.partonfilter+process.D0Daufilter)

process.ana = cms.EDAnalyzer('HydjetAnalyzer'
                             )
                             

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('treefile_pthat0.root')
                                                                      )
                                                                      

process.genAna = cms.Path(process.ana)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.genAna,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 









