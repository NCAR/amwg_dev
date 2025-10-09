#!/bin/bash

set -e

# Created 2025-10-09 09:50:44

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_097.FHISTC_LTso.ne30.baseline_DGLC_zref500.001"

/glade/work/hannay/cesm_tags/cam6_4_097/cime/scripts/create_newcase --compset HIST_CAM70%LT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_DGLC%NOEVOLVE_SWAV_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2176

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=11,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=special --force

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=11,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

