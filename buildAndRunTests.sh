#!/bin/sh
#
# Compiles and runs tests. Arguments:
#   1. target name (e.g. "runMediumTests" for runMediumTests.m)
#   2. path to MATLAB Compiler Runtime (e.g. "/share/apps/MATLAB/MATLAB_Compiler_Runtime/v714")

if [ ! -f "$1.m" ]
then
  echo No such file: $1.m
  exit 1
fi

#build
./build.sh $1
if [[ $? -ne 0 ]]
then
  exit 1
fi

#run
build=$1
pathToRunTime=$2
shift; shift;
fname=$(mktemp $(dirname $(readlink -f $0))/tmp/test.out.XXXXXX)
./run.sh $build $pathToRunTime $fname $*

#parse output
result="$(tail -2 $fname)"
if [[ "$result" =~ "NO DISPLAY$" ]]
then
  result="$(tail -4 $fname)"
  if [[ "$result" =~ "[Pp][Aa][Ss][Ss][Ee][Dd] in [0-9. a-z]*[.\s]*" ]]
  then
    exit 0
  fi
else
  if [[ "$result" =~ "[Pp][Aa][Ss][Ss][Ee][Dd] in [0-9. a-z]*[.\s]*$" ]]
  then
    exit 0
  fi
fi
exit 1
