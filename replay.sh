#!/bin/bash

set -e

# Created 2022-07-25 16:50:35

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_e.B1850MOM.f09_L32_t061.cam6_cice5.016"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_e/cime/scripts/create_newcase --compset 1850_CAM60_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res f09_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange CAM_CONFIG_OPTS=-phys cam6

./case.build

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange STOP_N=5

./case.submit

