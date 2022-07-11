#!/bin/bash

set -e

# Created 2022-07-08 13:45:22

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.010"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM_c/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE_MOM6_MOSART_CISM2%GRIS-NOEVOLVE_SWAV_SESP_BGC%BDRD --res ne30pg3_t061 --case b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.010 --run-unsupported

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./case.setup

./xmlchange JOB_QUEUE=regular,STOP_N=1,STOP_OPTION=ndays

./xmlchange JOB_QUEUE=premium

./xmlchange RUN_TYPE=hybrid

./case.build

./xmlchange RUN_REFCASE=gmom.e23.GJRAv3.TL319_t061_zstar_N75.nuopc.cice6salt.001

./xmlchange RUN_REFDATE=0060-01-01

./xmlchange DOUT_S=False

./check_case

./case.submit

./xmlchange JOB_WALLCLOCK_TIME=01:00:00

./case.submit

./xmlchange RUN_TYPE=initial

./xmlchange RUN_TYPE=startup

./case.submit

./case.submit

./case.submit

./xmlchange STOP_OPTION=nyears

./xmlchange DOUT_S=True

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./xmlchange JOB_QUEUE=regular

./xmlchange STOP_N=3

./xmlchange RESUBMIT=9

./case.submit

