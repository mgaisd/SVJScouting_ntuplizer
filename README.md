# To check out test
```
cmsrel CMSSW_10_6_26 #you can use CMSSW_11_1_0
cd CMSSW_10_6_26/src
git clone -b mods https://github.com/SUEPPhysics/SUEPScouting.git PhysicsTools/SUEPScouting
git clone -b master https://github.com/tresreid/PatUtils.git PhysicsTools/PatUtils
```

# To setup and compile
do this each time you login
```
cd $CMSSW_BASE/src
cmsenv
scram b
```

# To run QCD for year={2018,2017 or 2016} 2016 turns off scouting for qcd as this is not availible in the UL samples
```
cmsRun SUEPScouting/test/ScoutingNanoAOD_cfg.py inputFiles=file:qcd.root outputFile=flatscouting_qcd.root maxEvents=-1 isMC=true era=<year>
```
# To run Signal (note that "SUEP" must be in the input file name)
```
cmsRun SUEPScouting/test/ScoutingNanoAOD_cfg.py inputFiles=file:SUEP.root outputFile=flatscouting_signal.root maxEvents=-1 isMC=true era=<year>
```
# To run Data
```
cmsRun SUEPScouting/test/ScoutingNanoAOD_cfg.py inputFiles=file:data.root outputFile=flatscouting_data.root maxEvents=-1 isMC=false era=<year>
```
