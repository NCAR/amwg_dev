#!/bin/bash

set -e

# Created 2022-08-09 13:08:32

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L32_t061.camdev_cice5.022"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_e/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.build --clean-all

./case.build

./case.submit

./case.setup --reset

./case.build --clean-all

./case.build

./case.submit

./case.submit

./case.submit

./xmlchange RESUBMIT=6

