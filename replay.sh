#!/bin/bash

set -e

# Created 2023-10-05 16:16:40

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16b.BLT1850.ne30_t232.050"

/glade/work/gmarques/cesm.sandboxes/cesm2_3_alpha16b_t025/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case b.e23_alpha16b.BLT1850.ne30_t232.050 --run-unsupported --project CESM0023

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange ROF2OCN_LIQ_RMAPNAME=/glade/work/gmarques/cesm/tx2_3/runoff_mapping/map_r05_to_tx2_3_nnsm_e250r250_230914.nc

./xmlchange GLC2OCN_LIQ_RMAPNAME=/glade/work/gmarques/cesm/tx2_3/runoff_mapping/map_gland4km_to_tx2_3_nnsm_e250r250_230914.nc

./case.setup

./case.build

./check_case

./pelayout

./case.submit

./case.submit

./case.setup --reset

./case.build --clean-all

./case.build

./case.submit

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.build

./case.submit

./case.submit

./case.build

./case.submit

./case.submit

./case.submit

