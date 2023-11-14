#!/bin/bash

set -e

# Created 2023-10-10 14:27:09

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16b.BLT1850.ne30_t232.053"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16b/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange ROF2OCN_LIQ_RMAPNAME=/glade/work/gmarques/cesm/tx2_3/runoff_mapping/map_r05_to_tx2_3_nnsm_e250r250_230914.nc

./xmlchange GLC2OCN_LIQ_RMAPNAME=/glade/work/gmarques/cesm/tx2_3/runoff_mapping/map_gland4km_to_tx2_3_nnsm_e250r250_230914.nc

./case.build

./case.build

./case.build

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=12,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=12,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange RESUBMIT=10

