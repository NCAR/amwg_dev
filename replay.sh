#!/bin/bash

set -e

# Created 2025-06-23 14:28:26

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_097.FHISTMT_LTso.ne30_setting166_smooth-phis.001"

/glade/work/hannay/cesm_tags/cam6_4_097/cime/scripts/create_newcase --compset FHISTC_MTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./preview_namelists

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=1,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

