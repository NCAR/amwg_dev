#! /bin/csh -f

# -----------------------------------------------------------------------------
# NOTICE: This script was custom-generated for /glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_019_plus_CESM2.2.002_zm2_zmke_4.hf.
#
#         DO NOT COPY to another case!
#
#         Use this script ONLY in /glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_019_plus_CESM2.2.002_zm2_zmke_4.hf.
# -----------------------------------------------------------------------------

set sourcemod_dir  = /glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_019_plus_CESM2.2.002_zm2_zmke_4.hf/SourceMods/src.cism/source_cism
set cism_confIOdir = /glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_019_plus_CESM2.2.002_zm2_zmke_4.hf/Buildconf/cismIOconf

# -----------------------------------------------------------------------------
# NOTE: If you are viewing this script within the bld subdirectory of the cism
# code directory, please note that this is not a complete script. Instead, it
# is embedded in a script that is created by cism.cpl7.template (in the parent
# directory). That is where some variables are defined (cism_confIOdir,
# sourcemod_dir). This is done because cism.cpl7.template has access to the
# CASEROOT and CASEBUILD environment variables, whereas this script (which is
# meant to be run as a standalone script -- NOT part of the cesm build) does
# not necessarily know the values of these variables.
#
# If you are viewing this script from within your CASE directory, then the
# above note does not apply.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# generate new CISM _io.F90 files
# -----------------------------------------------------------------------------
cd $cism_confIOdir

# NOTE(wjs, 2015-04-03): glint is no longer used by CESM. However, I'm keeping
# the glint stuff here for now so that we can keep the glint default i/o files
# up-to-date (since I use this mechanism to regenerate the default i/o files)
foreach file (glide glad glad_mbal glint glint_mbal)
  set file_varsdef = ${file}_vars.def
  set file_ioF90 = ${file}_io.F90
  if (-f ${file_varsdef}) then
    # ---------------------------------------------------------------------------
    #  create new _io.F90 file using CISM's python script
    # ---------------------------------------------------------------------------
    $PYTHON generate_ncvars.py $file_varsdef ncdf_template.F90.in

    if (-f ${file_ioF90}) then
      # ---------------------------------------------------------------------------
      #  compare new _io.F90 file with current version in the objdir (if it exists)
      #  if different, copy the new one to the objdir
      # ---------------------------------------------------------------------------
      cp ${file_ioF90} ${sourcemod_dir}/${file_ioF90}
    else
      # ---------------------------------------------------------------------------
      #  if new _io.F90 file not created for some reason, exit
      # ---------------------------------------------------------------------------
      echo ERROR: CISM python script failed to produce new file: ${file_ioF90}
      exit 2
    endif

  else
    echo ERROR: missing CISM variable definition file: ${file_varsdef}
    exit 2
  endif
end


