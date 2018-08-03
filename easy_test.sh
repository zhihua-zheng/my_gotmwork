#!/bin/bash

cd /Users/Danny/Documents/GitLab/GOTM_dev/gotmwork
./build_src.sh -clean
./build_src.sh -build
cd /Users/Danny/Documents/GitLab/GOTM_dev/gotmwork/cases/OSMOSIS
./case_run
