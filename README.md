# Fluctuating-Selection

The computer program FlucSel22.cpp estimates the long-term steady-state distribution of allele frequencies for a biallelic locus under the joint forces of drift, mutation, and fluctuating selection (with a designated mean and standard deviation, assumed to be Gaussian distributed with no temporal autocorrelation, and applied randomly to each allele each generation).

The underlying details of the population-genetic model, which is Wright-Fisher in form, can be found at the top of the code.

On lines 58-76 of the code, the user enters the mutation bias, the numbers of major and minor sites, the mean and temporal SD of the selection coefficient at each site, the number of generations elapsed between censuses, the burnin period before statistics are taken, and the number of inervals between printouts of the cumulative statistics to the slurm files.

Although the program can, in principle, be run with large numbers of both types of loci, a number of issues will have to be resolved as to how selection should vary among different loci. So, for now, it is best to set the numbers of both loci = 1, and treat one as a neutral site. 

Lines 111, 118, and 284 can be edited to name the output files.

The code is written so that runs are made in parallel for 21 different combinations of population sizes and mutation rates, which can be modified on lines 289 to 340. Lines 347 to 369 allow the user to designate how long the runs should proceed; this can be checked to the user’s satisfaction that the final results have stabilized. 

If it is desired to alter the number of runs, edit line 116; and if >21 are to be run, the above-noted arrays will need to be expanded.

There will be two sets of output files: “slurm” files containing the cumulative statistics allow the user to determine whether the runs have proceeded for sufficient time to equilibrate; “dataout” files give the final results for each of the 21 conditions. 


To run the program, enter on the unix command line:

module load intel/2019.4
icc -o FlucSel22 FlucSel22.cpp -lm -lgsl
sbatch –array=1-21 FlucSel22.sh

(The first line may need to be modified, depending on the system involved. The 21 in the final line needs to be edited if a different number of runs is being made).


A batch shell file (FlucSel22.sh) must be provided in the local file space to set up the series of runs. For example:

#! /bin/bash

#SBATCH -A mlynch11
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH -n 1
#SBATCH -t 11-4:00

echo "Running the script in parallel"
./FlucSel22 $SLURM_ARRAY_TASK_ID


Here, the #SBATCH lines will need to be modified to the user’s specifications, FlucSel22 is the folder within which the .cpp and .sh files sit.

If all is operating properly, upon submission of the job, the slurm and dataout files should immediately appear, and the slurm files will begin to be periodically updated with the cumulative statistics following a burn-in period, with the dataout files becoming populated after each run completes. The completion times will vary depending on the user-defined run lengths.  


Summary file of results:

Upon completion of all runs, a sbatch file called concatenate.sh (for example; lines below) can be submitted, which will create an output called summary.txt that contains the comma-delimited stacked set of final results, which can be imported to a spreadsheet. “dataoutxxx” needs to be edited to give the prefix of the output file names, which will come out in parallel as dataoutxxx_1.txt to dataoutxxx_9.txt.

To run this file, type on the command line: sbatch concatenate.sh, making the appropriate change to the username. Summary.txt will then appear in your file space. 


#! /bin/bash

#SBATCH -A mlynch11

#SBATCH -n 1
#SBATCH -t 0-4:00

myFiles1=`ls dataoutxxx_?.txt`
myFiles2=`ls dataoutxxx_??.txt`
cat $myFiles1 $myFiles2 > summary.txt

 
