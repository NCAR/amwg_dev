#!/bin/bash

set -e

# Created 2023-04-13 11:51:38

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_106.FLTHIST_v0a.ne30.dcs_non-ogw_lcp-moistF.001"

/glade/work/hannay/cesm_tags/cam6_3_106/cime/scripts/create_newcase --compset FLTHIST_v0a --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_100.FWscHIST.ne30_L58.001

./xmlchange RUN_REFDATE=1991-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=1,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

