#!/bin/bash

set -e

# Created 2024-09-27 19:55:59

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_036.FMTHIST.ne120_L93.001"

/glade/work/hannay/cesm_tags/cam6_4_036/cime/scripts/create_newcase --compset FMTHIST --res ne120pg3_ne120pg3_mt13 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange NTASKS=6144

./case.setup

./preview_namelists

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=p03010039,RESUBMIT=10,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=p03010039

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.build

./case.submit

./case.build

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=p03010039,RESUBMIT=10,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=p03010039

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.submit

./case.submit

./xmlchange CONTINUE_RUN=TRUE

./case.submit

./xmlchange CONTINUE_RUN=FALSE

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=p03010039,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=p03010039

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./xmlchange CONTINUE_RUN=FALSE

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=p03010039,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=p03010039

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./xmlchange PROJECT=CESM0024

./case.submit

./case.build

./case.submit

./case.submit

./case.submit

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.submit

./xmlchange JOB_PRIORITY=regular

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.submit

./xmlchange JOB_PRIORITY=premium

./case.submit

./preview_namelists

./case.submit

./case.submit

./xmlchange NTASKS=6144

./xmlchange NTASKS=7680

./case.setup --reset

./case.build --clean

./preview_namelists

./case.build

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.submit

./case.submit

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./preview_namelists

./xmlchange CONTINUE_RUN=FALSE

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=29,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.build

