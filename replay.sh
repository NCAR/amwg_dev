#!/bin/bash

set -e

# Created 2024-07-11 16:55:26

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e23_alpha17f.FLTHIST_ne30.moving_mtn.099"

/glade/work/hannay/cesm_tags/cesm2_3_alpha17f/cime/scripts/create_newcase --compset HIST_CAM%DEV%LT%GHGMAM4_CLM51%SP_CICE%PRES_DOCN%DOM_MOSART_CISM2%GRIS-NOEVOLVE_SWAV --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./preview_namelists

./case.build

./preview_namelists

./case.build

./case.build

./case.build

./case.build

./xmlchange PROJECT=P93300642,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange CHARGE_ACCOUNT=P93300642

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./case.build --clean-all

./case.build

./case.submit

