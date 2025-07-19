# skimming_3_stations_MUonE

On lxplus8:

cd /afs/cern.ch/work/p/pasenov/ML_for_MUonE/skimming

source /afs/cern.ch/work/p/pasenov/fairsoft/install/bin/FairRootConfig_forCondor.sh
OR:
source /afs/cern.ch/work/p/pasenov/fairsoft/install/bin/FairRootConfig.sh

root -l -b -q 'runConvertMC.C("minbias6969.root")'

Take the produced file MC ROOT file and run the macro that converts it to datalike format. Change the path in that macro. It is called runConvertMC.C

compile.sh : change paths from ehess to pasenov and the folder of my FairMUonE installation

mkdir run_6969

OR

mkdir run_6976

mv /eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/minbias6969_convertedTFormat.root /eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/run_6969/muedaq01-6969.root

OR:
mv /eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/minbias6976_convertedTFormat.root /eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/run_6976/muedaq01-6976.root

In the runs folder make a text file called, e.g., run_6969_chunks_dgm.txt

cp run_6969_chunks_dgm.txt run_6969_files.txt

In this file, add the list of the ROOT files you want to use in the form:
muedaq01-6969.root

skim.cc : change 	string filepath = "/eos/experiment/mu-e/staging/daq/2025/decoded/8/"; to /eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/ the folder containing the folder run_6976

source ./compile.sh

source ./Run.sh run_6969 0 1001
       
OR:

source ./Run.sh run_6976 0 1001

