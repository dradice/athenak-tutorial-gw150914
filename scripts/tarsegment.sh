#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "Usage: $0 /segment/to/archive"
    exit 0
fi
segment=$1

set -e
set -x

tar cf ${segment}.min.tar \
    --exclude ${segment}/bin \
    --exclude ${segment}/horizon_0 \
    --exclude ${segment}/horizon_1 \
    --exclude ${segment}/cce \
    --exclude ${segment}/rst \
    ${segment}
tar cf ${segment}.bin.tar ${segment}/bin
if [ -d ${segment}/cce ]; then
    tar cf ${segment}.cce.tar ${segment}/cce
fi
if [ -d ${segment}/horizon_0 ]; then
    tar cf ${segment}.hor.tar ${segment}/horizon_?
fi
tar cf ${segment}.rst.tar ${segment}/rst/$(ls -r ${segment}/rst | head -n 1)
touch ${segment}.tar.done
