#!/bin/bash

set -e

# Created 2025-10-10 15:41:21

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_118.FHISTC_LTso.ne120_L58.test001"

/glade/work/hannay/cesm_tags/cam6_4_118/cime/scripts/create_newcase --compset FHISTC_LTso --res ne120pg3_ne120pg3_mt13 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=6144

./case.setup

