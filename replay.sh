#!/bin/bash

set -e

# Created 2025-07-25 18:17:19

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_078.FHISTC_LTso.ne30.baseline156_icefrac0.5.001"

/glade/work/hannay/cesm_tags/cam6_4_078/cime/scripts/create_newcase --compset HIST_CAM70%LT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_DGLC%NOEVOLVE_SWAV_SESP --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

