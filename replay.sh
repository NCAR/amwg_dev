#!/bin/bash

set -e

# Created 2023-06-09 17:30:47

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_112.FLT1850.ne30.landspinup.001"

/glade/work/hannay/cesm_tags/cam6_3_112/cime/scripts/create_newcase --compset 1850_CAM%DEV%LT%GHGMAM4_CLM51%BGC-CROP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./case.setup

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.build

./case.submit

./case.build

./case.build

./case.build --clean-all

./case.build

./case.submit

./case.build

./case.submit

./case.submit

./case.submit

./case.submit

./xmlchange RESUBMIT=15

