#!/bin/sh
#
# Compiles and runs tests. Arguments:
#   1. target name (e.g. "runMediumTests" for runMediumTests.m)
#   2. path to MATLAB Compiler Runtime (e.g. "/share/apps/MATLAB/MATLAB_Compiler_Runtime/v714")
#   3. File to write MATLAB stdout 

pathToBuildAndBuild=$1
if [[ $pathToBuildAndBuild =~ "^(.*?)/([^/]+)$" ]]
then
  build=${BASH_REMATCH[2]}
  pathToBuild=${BASH_REMATCH[1]}
else
  build=$pathToBuildAndBuild
  pathToBuild='bin'
fi

if [ ! -f "$build.m" ]
then
  echo No such file: $build.m
  exit 1
fi

echo Running $build...
export MCR_CACHE_ROOT=/tmp/mcr_cache_$$
export XAPPLRESDIR=/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/X11/app-defaults
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/runtime/glnxa64:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/os/glnxa64:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64
mkdir -p $MCR_CACHE_ROOT
rm -rf $MCR_CACHE_ROOT/*

pathToRunTime=$2
outFile=$3
shift; shift; shift;
if [ "$outFile" == "" ]; then
  $pathToBuild/$build/run_$build.sh $pathToRunTime
else
  $pathToBuild/$build/run_$build.sh $pathToRunTime $* | tee $outFile
fi
status=$?

rm -rf $MCR_CACHE_ROOT

if [[ $status -ne 0 ]]
then
  exit 1
fi
exit 0
