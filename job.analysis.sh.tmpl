#!/bin/sh

#job
#PBS -N runAnalysis-<TMPL_VAR NAME=conditionSetTimeStamp>-<TMPL_VAR NAME=iJob>

#user
#PBS -P <TMPL_VAR NAME=linuxUser>:<TMPL_VAR NAME=linuxUser>

#notification
#PBS -M <TMPL_VAR NAME=emailAddress>
#PBS -m a

#resources
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=22gb
#PBS -l vmem=22gb

#log
#PBS -o <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/out.analysis-<TMPL_VAR NAME=iJob>.log
#PBS -e <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/err.analysis-<TMPL_VAR NAME=iJob>.log
#PBS -W umask=002

#schedule
#PBS -W depend=afterany<TMPL_VAR NAME=afterany>

#set environment
export MATLAB_PREFDIR=/tmp/emptydir
export MCR_CACHE_ROOT=/tmp/mcr_cache_$PBS_JOBID
mkdir -p $MCR_CACHE_ROOT

#setup
cd <TMPL_VAR NAME=baseDir>

#set permissions
sudo chown -R <TMPL_VAR NAME=linuxUser>:<TMPL_VAR NAME=linuxUser> <TMPL_VAR NAME=outDir>/singleGeneDeletions/
chmod -R 775 <TMPL_VAR NAME=outDir>/singleGeneDeletions/

#run analysis
<TMPL_IF NAME=firstStageAnalysis>
./run.sh <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/bin/runAnalysis \
  <TMPL_VAR NAME=pathToRunTime> \
  <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/matlab.analysis.log \
  <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp> <TMPL_VAR NAME=iJob> <TMPL_VAR NAME=nFirstStageJobs>
</TMPL_IF>

#compile summary
./compileSummary.pl <TMPL_VAR NAME=simulationIdx> <TMPL_VAR NAME=conditionSetTimeStamp>

#set permissions
sudo chown -R <TMPL_VAR NAME=linuxUser>:<TMPL_VAR NAME=linuxUser> <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/
chmod -R 775 <TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/

sudo chown -R <TMPL_VAR NAME=linuxUser>:<TMPL_VAR NAME=linuxUser> <TMPL_VAR NAME=outDir>/singleGeneDeletions/
chmod -R 775 <TMPL_VAR NAME=outDir>/singleGeneDeletions/

#mail results
<TMPL_IF NAME=firstStageAnalysis>
<TMPL_ELSE>
./mailSimulationsSummary.pl <TMPL_VAR NAME=simulationIdx> <TMPL_VAR NAME=conditionSetTimeStamp>
</TMPL_IF>

#cleanup
rm -rf $MCR_CACHE_ROOT/*

#resources
echo ""
echo "=============="
echo "=== status ==="
echo "=============="
qstat -f $PBS_JOBID

#status
if [[ -f "<TMPL_VAR NAME=outDir>/<TMPL_VAR NAME=conditionSetTimeStamp>/err.analysis.log" ]]
then
  exit 1
fi
exit 0
