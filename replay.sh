#!/bin/bash

set -e

# Created 2021-12-07 10:51:36

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_035.tphysac_zm2_reorder_clubbdiffusion.001.hf2"

/glade/work/hannay/cesm_tags/cam6_3_035.tphysac/cime/scripts/create_newcase --compset HIST_CAM60%WCSC_CLM50%BGC-CROP_CICE%PRES_DOCN%DOM_MOSART_CISM2%NOEVOLVE_SWAV_SIAC_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS= -pcols 9 --append

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -age_of_air_trcs -chem waccm_sc_mam4 -cppdefs -Dwaccm_debug -nlev 58

./xmlchange CAM_CONFIG_OPTS

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev --append

./xmlchange CAM_CONFIG_OPTS

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -age_of_air_trcs -chem waccm_sc_mam4 -cppdefs -Dwaccm_debug -nlev 58

./case.setup

./xmlchange RUN_STARTDATE=1979-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.e21.FWscHIST_BGC.ne30_ne30_mg17_L48_revert-J.001

./xmlchange RUN_REFDATE=1989-01-01

./xmlchange GET_REFCASE=TRUE

./case.build

./case.build

./case.build

