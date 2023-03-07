#!/bin/bash

set -e

# Created 2022-10-25 19:49:02

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026g"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_e/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE5_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026g --run-unsupported --project CESM0002

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_QUEUE=premium

./case.build

./case.setup

./case.build

./case.setup --reset

./case.build

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./case.build

./case.setup --reset

./case.build

./case.build --clean-all

./case.build

./case.submit

./case.submit

./xmlchange JOB_QUEUE=premium,RESUBMIT=10,STOP_N=3,STOP_OPTION=nyears

./case.submit

./case.submit

./xmlchange RESUBMIT=11

./case.submit

./xmlchange JOB_QUEUE=regular

./xmlchange RESUBMIT=9

./case.submit

./case.submit

./xmlchange PROJECT=cesm0023

./xmlchange RESUBMIT=12

./case.submit

./xmlchange RESUBMIT=10

./case.submit

./xmlchange RESUBMIT=12

./case.submit

./xmlchange RESUBMIT=9

./case.submit

./case.submit

