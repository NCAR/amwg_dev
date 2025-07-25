#!/bin/bash

set -e

# Created 2025-07-25 16:29:57

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_beta06.FHISTC_LTso.f09.b183.001"

/glade/derecho/scratch/cacraig/cam6_4_089/cime/scripts/create_newcase --compset HIST_CAM70%LT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_DGLC%NOEVOLVE_SWAV_SESP --res f09_f09_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

