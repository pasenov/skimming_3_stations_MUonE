#!/bin/bash

input_run=$1 # run_x
nsplit=$2    # split in groups of nsplit files; 0: do not split 
iskim=$3     # 1: skim clean single-mu int.; 40: PU (2+3+>=4) int.; 41=1+40; 60: PU (2+3+4) int.; 61=1+60

DIR=`pwd`
ORDER_INPUT=${DIR}/order_input
SPLIT=${DIR}/split_input
SKIM=${DIR}/skim
SUMUP=${DIR}/add_output
runlist=${DIR}/runs/${input_run}_chunks_dgm.txt
run=${DIR}/runs/${input_run}_files.txt

echo "INPUT RUN: " ${input_run}

if [[ (-f ${runlist}) ]]; then
    echo "Input run file list: " $runlist
else
    echo "Do not find input for:"$runlist
    return 1
fi

WORKDIR=job_`date +"%y-%m-%d_%T"`
echo "Working directory is: " ${WORKDIR}
mkdir ${WORKDIR}
cd ${WORKDIR}
ln -s ../runs .

mkdir inputs; mkdir outputs

indir=inputs/${input_run}
mkdir $indir
outdir=outputs/${input_run}
mkdir $outdir

#########################################
echo "ordering input files along SuperID"
${ORDER_INPUT}  ${input_run}
#########################################

if [ $nsplit -ne 0 ]; then
###########################################################
    echo "splitting input in groups of " ${nsplit} " files"
    ${SPLIT}  ${input_run}  ${nsplit}
###########################################################
else
    cp ${run} inputs/${input_run}/inputs
fi

fname="histos"
skimname=""

case ${iskim} in
    0)
	echo "analysis without saving skims..."
	;;
    1)
	echo "skimming clean single-mu interactions..."
	skimname="skim_interactions_singlemu"
        ;;
	1001)
	echo "skimming clean single-mu interactions for MC..."
	skimname="skim_interactions_singlemu_MC"
        ;;
    40)
	echo "skimming clean pile-mu interactions..."
	skimname="skim_interactions_pileupmu"
	;;
    41)
	echo "skimming any clean single or pileup mu interactions..."
	skimname="skim_interactions_anymu"
	;;
    60)
	echo "skimming clean 2,3,4 pile-mu interactions..."
	skimname="skim_interactions_pileup234mu"
	;;
    61)
	echo "skimming any clean single or pileup mu interactions..."
	skimname="skim_interactions_any1234mu"
	;;
    100)
	echo "skimming single golden mu in S0..."
	skimname="skim_goldenS0"
	;;
    101)
	echo "skimming single golden mu in S0 (preselected with one hit in S1)..."
	skimname="skim_goldenS0_presel"
	;;
    *)
	echo "Select a skim type among: 1, 40, 41, 60, 61, 100, 101 or 0 for no skim !";;
esac

#-----------------------------------------------------------------
if [ $nsplit -eq 0 ]; then
    
    file=inputs
    ##############################################################
    time ${SKIM}  ${input_run} ${file} ${iskim}   >&  ${file}.log
    ##############################################################
    
else

    i=0
    while [[ -f  inputs/${input_run}/inputs_$i ]] 
    do
        file=inputs_$i	
	##############################################################
	time ${SKIM}  ${input_run} ${file} ${iskim}   >&  ${file}.log
	##############################################################
	
	hfile=outputs/${input_run}/${fname}_${input_run}_inputs_$i.root
	hfile_old=outputs/${input_run}/${fname}_${input_run}_sumold.root
	hfile_base=outputs/${input_run}/${fname}_${input_run}_sum
	hfile_new=${hfile_base}.root	
	
	if [ $i -eq 0 ]; then
	    cp $hfile $hfile_new
	    cp relative_efficiency_S1_goodS0.txt relative_efficiency_S1_goodS0_sum.txt
	    cp relative_efficiency_S0_goodS1.txt relative_efficiency_S0_goodS1_sum.txt
	else
	    mv $hfile_new $hfile_old
	    #+++++++++++++++++++++++++++++++++++++++++++++++++++
	    $ROOTSYS/bin/hadd -k  $hfile_new  $hfile_old  $hfile
	    #+++++++++++++++++++++++++++++++++++++++++++++++++++
	    
	    mv relative_efficiency_S1_goodS0_sum.txt relative_efficiency_S1_goodS0_sumold.txt
	    mv relative_efficiency_S0_goodS1_sum.txt relative_efficiency_S0_goodS1_sumold.txt
	    
	    histpath=outputs/${input_run}
	    histbase=${fname}_${input_run}_sum
	    total=skim_${input_run}_counts
	    
            ######################################
	    $SUMUP  $histpath/$histbase  >& $total
            ######################################
	fi
	
	((i++))
	
#	if [ $i -gt 1 ]; then
#	    echo "*** i = " $i
#	    break
#	fi
	
    done
    
fi
#-----------------------------------------------------------------

if [ $nsplit -ne 0 ]; then
    rm ${outdir}/*sumold.root
    rm relative_efficiency*sumold.txt
elif [ $iskim -gt 0 ]; then
    # add the BranchList to automatically splitted files (by ROOT)
    basename=${outdir}/${skimname}_${input_run}_inputs
    nfiles=`ls ${basename}* | wc -l`

    if [ $nfiles -gt 1 ]; then
	mv ${basename}.root ${basename}_0.root
	root -l -b -q ''${DIR}'/addBranchList.C("'$basename'", '$nfiles')' # horrible syntax
    fi
fi

echo "all Done!"
