#!/bin/bash

set -e

# Created 2022-05-23 09:24:19

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam041_mom.B1850WcMOM.ne30_L58_t061.cice5.cam6.005"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_041_MOM3/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam6 -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./xmlchange ATM_GRID=ne30np4.pg3

./case.setup

./case.build

./case.build

./case.build

./case.build

./case.build

./case.submit

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=12,STOP_OPTION=nmonths

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=12,STOP_OPTION=nmonths

./case.submit

./case.build

./case.build

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./case.submit

./preview_namelists

./preview_namelists

./case.submit

./xmlchange RESUBMIT=4

./xmlchange STOP_N=1,RESUBMIT=13

./xmlchange CREATE_ESMF_PET_FILES=TRUE

./xmlchange JOB_QUEUE=premium

./xmlchange JOB_QUEUE=premium,STOP_N=1

./case.submit

./xmlchange JOB_QUEUE=regular

./case.submit

./case.setup --reset

./case.build --clean

./case.build

./case.submit

