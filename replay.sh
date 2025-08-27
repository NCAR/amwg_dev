#!/bin/bash

set -e

# Created 2025-08-27 17:03:56

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_107.FHISTC_LTso.ne30.lwoff.001"

/glade/work/hannay/cesm_tags/cam6_4_107/cime/scripts/create_newcase --compset FHISTC_LTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

