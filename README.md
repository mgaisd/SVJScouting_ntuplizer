# SVJScouting_ntuplizer
Produce ntuples in a PFNanoAOD-like format from AOD or MiniAOD


# Introduction
The code is based on the [repository](https://github.com/mgaisd/SVJScouting) and runs on MiniAOD to produce NanoAOD-like format files.
For the production of Run2 scouting MiniAOD for SVJ signal look at the [repository](https://github.com/cms-svj/SVJProduction).


# To check out
```
cmsrel CMSSW_10_6_26 
cd CMSSW_10_6_26/src
git cms-init
git clone -b main https://github.com/CMS-SVJ-scouting/SVJScouting_ntuplizer.git PhysicsTools/SVJScouting
git clone -b master https://github.com/tresreid/PatUtils.git PhysicsTools/PatUtils
git cms-addpkg CommonTools/PileupAlgos
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg JetMETCorrections
git cms-addpkg CommonTools/RecoAlgos
git cms-addpkg Configuration/ProcessModifiers

# for importing theory weights for backgrounds (not really necessary and not tested)
git cms-merge-topic -u https TreeMaker:GenWeightRefactor_106X_squash
```

# To setup and compile
```
cd $CMSSW_BASE/src
cmsenv
scram b
```

# Running the ntuplizer from MiniAOD or AODs to NanoAOD-like format
To run from MiniAODs:

```
cmsRun SVJScouting/test/ScoutingNanoAOD_fromMiniAOD_cfg.py inputFiles=file:miniaod_file.root outputFile=flatscouting_signal.root maxEvents=-1 isMC=True era=<year> signal=True MatrixElementInfo=True
```

To run from AODs:

```
cmsRun SVJScouting/test/ScoutingNanoAOD_fromAOD_cfg.py inputFiles=file:aod_file.root outputFile=flatscouting_background.root maxEvents=-1 isMC=True era=<year> 
```


To run from Scouting data:

```
cmsRun SVJScouting/test/ScoutingNanoAOD_fromData_cfg.py inputFiles=file:data_file.root outputFile=flatscouting_data.root maxEvents=-1 era=<year>
```

To produce a reduced ScoutingAOD format containing only scouting collections:
```
cmsRun SVJScouting/test/ScoutingAOD_cfg.py inputFiles=file:aod_file.root outputFile=scouting_aod_file.root maxEvents=-1
```

To run the ntuplizer on the ScoutingAOD format:
```
cmsRun SVJScouting/test/ScoutingNanoAOD_fromAOD_cfg.py inputFiles=file:scouting_aod_file.root outputFile=flatscouting_background.root maxEvents=-1 isMC=True era=<year> onlyScouting=True
```
