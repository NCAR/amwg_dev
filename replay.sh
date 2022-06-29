#!/bin/bash

set -e

# Created 2022-05-23 14:16:34

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam041_mom.B1850WcMOM.ne30_L58_t061.cice5.cam6.006"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_041_MOM3/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam6 -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./xmlchange ATM_GRID=ne30np4.pg3

./case.setup

./case.build

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./case.submit

./preview_namelists

./xmlchange RESUBMIT=4

./xmlchange RESUBMIT=3

./xmlchange RESUBMIT=10,STOP_N=1

./case.build

./case.submit

./case.submit

./case.setup --reset

./case.build --clean

./case.build

./case.submit

./case.build

./case.submit

./case.setup --reset

./case.build --clean

./case.build

./case.submit

./xmlchange JOB_QUEUE=premium

./case.submit

./xmlchange CONTINUE_RUN=TRUE

./xmlchange CONTINUE_RUN=TRUE,STOP_N=1,STOP_OPTION=nyears

./xmlchange CONTINUE_RUN=TRUE,STOP_N=1,STOP_OPTION=nyears,JOB_QUEUE=regular,RESUBMIT=10

./case.submit

./xmlchange JOB_QUEUE=premium

./case.submit

./case.submit

./xmlchange JOB_QUEUE=regular

./xmlchange RESUBMIT=0

