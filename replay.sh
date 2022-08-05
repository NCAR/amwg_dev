#!/bin/bash

set -e

# Created 2022-08-05 13:35:25

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_e.B1850MOM.ne30_L32_t061.camdev_cice5.021"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_e/cime/scripts/create_newcase --compset 1850_CAM60_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2

./case.build

./preview_namelists

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

