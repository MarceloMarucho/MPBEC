#!/bin/bash
var1="$(find /home -type d -name mpbecini 2> /dev/null)/mpbecini.m"
echo "$var1"
matlab -nosplash -nodesktop -r "run $var1"
