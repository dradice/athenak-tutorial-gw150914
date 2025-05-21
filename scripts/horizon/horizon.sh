#!/usr/bin/env bash

ANALYSIS_HOME=$(cd $(dirname $0)/../..; pwd)

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/simulation"
    exit 0
fi
TARGET_SIMULATION=$1

set -e
set -x

cd $TARGET_SIMULATION
for segment in output-????; do
    if [ ! -f ${segment}.hor.done ]; then
        for hor in $segment/horizon_?; do
            for out in $hor/output_*; do
                (cd $out && 
                ${ANALYSIS_HOME}/ahfind/exe/cactus_einsteintoolkitanalysis \
                    ET_analyze_BHaH_data_horizon.par &> ET_analysis.log
                )
            done
        done
        touch ${segment}.hor.done
    fi
done

mkdir -p horizons
${ANALYSIS_HOME}/scripts/workflow/collate.py -i \
    -o horizons/horizon_BH_0_ahf_ihf_diags.txt \
    output-????/horizon_0/output_*/horizon_BH_0_ahf_ihf_diags.txt
${ANALYSIS_HOME}/scripts/workflow/collate.py -i \
    -o horizons/horizon_BH_1_ahf_ihf_diags.txt \
    output-????/horizon_1/output_*/horizon_BH_1_ahf_ihf_diags.txt
