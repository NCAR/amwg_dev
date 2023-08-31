#!/bin/bash

set -e

# Created 2023-08-31 09:42:09

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.e23_alpha16b.FLT2000.ne30.033"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16b/cime/scripts/create_newcase --compset 2000_CAM%DEV%LT%GHGMAM4_CLM51%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001

./xmlchange RUN_REFDATE=1994-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./case.build

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=1,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.submit

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=0001-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001

./xmlchange RUN_REFDATE=1994-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=1,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.build

./case.build

./case.build --clean-all

./case.build

./case.submit

