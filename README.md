
# Introduction
The code is based on the [repository](https://github.com/mgaisd/SVJScouting) and runs on MiniAOD to produce NanoAOD-like format files.
For the production of Run2 scouting MiniAOD for SVJ signal look at the [repository](https://github.com/cms-svj/SVJProduction).


# To check out
```
cmsrel CMSSW_10_6_26 #you can use CMSSW_11_1_0
cd CMSSW_10_6_26/src
git cms-init
git clone -b mods https://github.com/cesarecazzaniga/SVJScouting.git PhysicsTools/SVJScouting
git clone -b master https://github.com/tresreid/PatUtils.git PhysicsTools/PatUtils
git cms-addpkg CommonTools/PileupAlgos
git cms-addpkg PhysicsTools/PatAlogs
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
cmsRun SVJScouting/test/ScoutingNanoAOD_fromMiniAOD_cfg.py inputFiles=file:miniaod_file.root outputFile=flatscouting_signal.root maxEvents=-1 isMC=true era=<year> signal=True
```

To run from AODs:

```
cmsRun SVJScouting/test/ScoutingNanoAOD_fromAOD_cfg.py inputFiles=file:miniaod_file.root outputFile=flatscouting_signal.root maxEvents=-1 isMC=true era=<year> signal=True
```

