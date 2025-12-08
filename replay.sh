#!/bin/bash

set -e

# Created 2025-12-08 11:22:22

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_134.FHISTC_LTso.ne30.baseline_DGLC.001"

/glade/work/hannay/cesm_tags/cam6_4_134/cime/scripts/create_newcase --compset HIST_CAM70%LT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_DGLC%NOEVOLVE_SWAV_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2176

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

