#!/bin/bash

set -e

# Created 2025-12-09 14:49:36

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_134.FHISTC_LTso.ne30_L58.c6aero.001"

/glade/work/hannay/cesm_tags/cam6_4_134/cime/scripts/create_newcase --compset FHISTC_LTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=2176

./case.setup

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./preview_namelists

./case.build

