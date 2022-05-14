#!/bin/bash
var1="$(find ~ -name mpbecini 2> /dev/null)/mpbecini.m"
"$(find /Applications -name matlab 2> /dev/null | grep "bin")" -nosplash -nodesktop -r "run $var1"
