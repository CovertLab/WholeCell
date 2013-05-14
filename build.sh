#!/bin/sh
#
# Compiles and runs tests. Arguments:
#   1. target name (e.g. "runMediumTests" for runMediumTests.m)

if [ ! -f "$1.m" ]
then
  echo No such file: $1.m
  exit 1
fi

echo Compiling $1...

rm -rf bin/$1
mkdir -p bin/$1  
./runMcc.pl runSimulation.prj $1
chmod -R 775 bin/$1

if [ $? -ne 0 ]
then
  echo Compile failed.
  exit 1
fi

echo Compile succeeded
exit 0
