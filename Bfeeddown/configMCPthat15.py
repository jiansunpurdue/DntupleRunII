from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MCproductionsForBtoDpthat15'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'Pythia8_nonprompt_D0pt15p0_Pthat15_TuneCUETP8M1_5020GeV_cfi_evtgen130_py_GEN_SIM.py'

config.Data.outputPrimaryDataset = 'Pythia'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 2500  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2'
config.Data.publication = True
config.Data.outputDatasetTag = 'MCproductionsForBtoDpthat15'
config.Site.storageSite = 'T2_CH_CERN'
