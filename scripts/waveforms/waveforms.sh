#!/usr/bin/env bash

set -e
set -x

ANALYSIS_HOME=$(cd $(dirname $0); pwd)

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/simulation [tmax]"
    exit 1
fi
target=$1
shift

cd $target

WAVEFORMS=$(ls -1 output-0000/waveforms/)

mkdir -p waveforms
for path in $WAVEFORMS; do
    name=$(basename $path)
    python $ANALYSIS_HOME/../collate.py -i -o waveforms/$name \
        output-????/waveforms/$name $@
done

python $ANALYSIS_HOME/waveforms.py --verbose waveforms

touch waveforms.done
