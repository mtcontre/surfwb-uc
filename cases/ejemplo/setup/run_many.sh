#!/bin/bash
for i in 4096 8192 16384
do
  octave --eval 'setrun_dbseco('$i')'
  cd ..
  ./xsurf
  cd setup
done