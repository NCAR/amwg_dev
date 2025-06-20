#!/bin/bash

set -e

# Created 2025-06-20 15:36:04

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_097.FHISTMT_LTso.ne30_setting166.001"

/glade/work/hannay/cesm_tags/cam6_4_097/cime/scripts/create_newcase --compset FHISTC_MTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./xmlchange RUN_STARTDATE=1990-01-01

