#!/bin/bash

set -e

# Created 2022-06-06 13:41:14

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom.B1850WscMOM.ne30_L58_t061.007"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./case.setup

./case.build

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=10,STOP_OPTION=nyears

./xmlchange JOB_QUEUE=premium

./case.submit

./xmlchange JOB_QUEUE=regular

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=10,STOP_OPTION=nyears

./case.submit

./case.build

./xmlchange STOP_N=1

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.submit

./case.submit

