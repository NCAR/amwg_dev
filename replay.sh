#!/bin/bash

set -e

# Created 2021-10-19 16:19:05

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e21.BWsc1850.ne30_L48_BL10_cesm2_3_alpha05c_cam6_3_028_cam6_parcel_zm.004_zm2.hf"

/glade/work/hannay/cesm_tags/cesm2_3_alpha05c_cam6_3_028_cam6_parcel_zm/cime/scripts/create_newcase --compset 1850_CAM60%WCSC_CLM50%BGC-CROP_CICE_POP2%ECO_MOSART_CISM2%NOEVOLVE_WW3_SIAC_SESP_BGC%BDRD --res ne30pg3_g17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS_ATM=1152

./xmlchange NTASKS_CPL=1152

./xmlchange NTASKS_OCN=24

./xmlchange NTASKS_WAV=36

./xmlchange NTASKS_GLC=1152

./xmlchange NTASKS_ICE=288

./xmlchange NTASKS_ROF=828

./xmlchange NTASKS_LND=828

./xmlchange NTASKS_ESP=1

./xmlchange ROOTPE_ATM=0

./xmlchange ROOTPE_CPL=0

./xmlchange ROOTPE_OCN=1152

./xmlchange ROOTPE_WAV=1116

./xmlchange ROOTPE_GLC=0

./xmlchange ROOTPE_ICE=828

./xmlchange ROOTPE_ROF=0

./xmlchange ROOTPE_LND=0

./xmlchange ROOTPE_ESP=0

./xmlchange NTHRDS=1

./xmlchange NTASKS_OCN=144

./xmlchange NTASKS_ATM=1800,NTASKS_CPL=1800,NTASKS_LND=900,NTASKS_ROF=900,NTASKS_ICE=900

./xmlchange ROOTPE_OCN=1800,ROOTPE_ICE=900

./xmlchange NTASKS_LND=1080,NTASKS_ROF=1080,NTASKS_ICE=720,ROOTPE_ICE=1080

./xmlchange CAM_CONFIG_OPTS= -pcols 9 --append

./xmlchange CAM_CONFIG_OPTS=-phys cam6 -age_of_air_trcs -chem waccm_sc_mam4 -cppdefs -Dwaccm_debug -nlev 58

./case.setup

./xmlchange RUN_STARTDATE=0001-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=b.e21.B1850.f09_g17.CMIP6-piControl.001

./xmlchange RUN_REFDATE=0501-01-01

./xmlchange GET_REFCASE=TRUE

./preview_namelists

./xmlchange --noecho CPL_I2O_PER_CAT=TRUE

./preview_namelists

./case.build

./preview_namelists

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./case.build

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./case.build

./case.build --clean-all

./case.build

./case.submit

