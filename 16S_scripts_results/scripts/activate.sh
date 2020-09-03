#!/bin/bash

export PATH=/biodata/dep_psl/grp_psl/tools/bin/:/biodata/dep_psl/grp_psl/tools/usearch/:$PATH

export LD_LIBRARY_PATH=/biodata/dep_psl/grp_psl/tools/lib/:$LD_LIBRARY_PATH

export PYTHONPATH=/biodata/dep_psl/grp_psl/tools/lib/python2.7/site-packages/:$PYTHONPATH

echo -e "$(whoami)"\\t"$(date -u)".>> /biodata/dep_psl/common/garridoo/logs/activate.log

