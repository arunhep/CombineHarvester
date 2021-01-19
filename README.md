# CombineHarvester

Full documentation: http://cms-analysis.github.io/CombineHarvester

## Quick start

A new full release area can be set up and compiled in the following steps:

    export SCRAM_ARCH=slc7_amd64_gcc700
    scram project CMSSW CMSSW_10_2_13
    cd CMSSW_10_2_13/src
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit/
    git fetch origin
    git checkout v8.1.0
    git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
    scramv1 b -j4

