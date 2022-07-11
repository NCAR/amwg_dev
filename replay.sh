#!/bin/bash

set -e

# Created 2022-07-09 21:52:42

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_c.B1850WscMOM.f09_L58_t061.011"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_c/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res f09_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./case.setup

./case.build

./case.build

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange RESUBMIT=4

