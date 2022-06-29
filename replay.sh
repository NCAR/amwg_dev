#!/bin/bash

set -e

# Created 2022-05-04 10:23:57

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam041_mom.B1850WcMOM.ne30_L58_t061.005"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_041_MOM3/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./xmlchange ATM_GRID=ne30np4.pg3

./case.setup

./case.build

./preview_namelists

./case.build

./case.submit

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./case.submit

./preview_namelists

./case.submit

./xmlchange RESUBMIT=3

./case.submit

./xmlchange RESUBMIT=20

./case.submit

./xmlchange RESUBMIT=0

