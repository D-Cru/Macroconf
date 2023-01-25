#!/bin/bash -e

# Adapted from https://github.com/derrickstolee/sparse-checkout-example

SCRIPTNAME=$0
die() {
	echo "$SCRIPTNAME: $1"
	exit 1
}

TEAM=$1

case $TEAM in
"dataset")
	FOLDERS="data docs"
	;;
"workflow")
	FOLDERS="data docs workflow/docs workflow/envs workflow/hpc workflow/libs workflow/misc workflow/notebooks workflow/reports workflow/rules workflow/scripts workflow/src workflow/tests"
	;;
"workflow-examples")
	FOLDERS="data docs workflow/data/external workflow/data/interim/example workflow/data/processed/example workflow/docs workflow/envs workflow/hpc workflow/libs workflow/misc workflow/notebooks workflow/reports workflow/rules workflow/scripts workflow/src"
	;;
"workflow-tests")
	FOLDERS="data docs workflow/data/external workflow/data/interim/tests workflow/data/processed/tests workflow/docs workflow/envs workflow/hpc workflow/libs workflow/misc workflow/notebooks workflow/reports workflow/rules workflow/scripts workflow/src workflow/tests"
	;;
"workflow-full")
	FOLDERS="data docs workflow"
	;;
*)
	die "please specify a valid input: 'dataset', 'workflow', 'workflow-examples', 'workflow-full'"
	;;
esac

echo "Running 'git sparse-checkout init --cone'"
git sparse-checkout init --cone

echo "Running 'git sparse-checkout set $FOLDERS'"
git sparse-checkout set $FOLDERS

# if [ $TEAM = "workflow-examples" ]; then
# 	echo "Touching missing output files that are too big to supply via Github (full trajectories and md log files)"
# 	touch "workflow/data/interim/example/22/H2O/11_GaMD_full/2000/0/210a1ea8aa678b16_md.out"
# 	touch "workflow/data/interim/example/22/H2O/11_GaMD_full/2000/0/210a1ea8aa678b16_gamd.log"
# 	touch "workflow/data/interim/example/22/H2O/11_GaMD_full/2000/0/210a1ea8aa678b16_traj.netcdf"
# 	touch "workflow/data/interim/example/22/H2O/11_GaMD_full/2000/0/210a1ea8aa678b16_traj.ncdf"
# 	touch "workflow/data/interim/example/22/H2O/11_GaMD_full/2000/0/210a1ea8aa678b16_weights.dat"

# 	touch "workflow/data/interim/example/22/H2O/8_cMD/2000/0/586db4c575bef492_traj.netcdf"
# 	touch "workflow/data/interim/example/22/H2O/8_cMD/2000/0/586db4c575bef492_traj.ncdf"
# 	touch "workflow/data/interim/example/22/H2O/8_cMD/2000/0/586db4c575bef492_md.out"

# 	touch "workflow/data/interim/example/22/H2O/9_aMD_cMD/2000/0/3595ce0609206d95_md.out"
# 	touch "workflow/data/interim/example/22/H2O/9_aMD_cMD/2000/0/3595ce0609206d95_traj.netcdf"
# 	touch "workflow/data/interim/example/22/H2O/9_aMD_cMD/2000/0/3595ce0609206d95_traj.ncdf"
# 	touch "workflow/data/interim/example/22/H2O/9_aMD_cMD/2000/0/3595ce0609206d95_Epot.dat"
# 	touch "workflow/data/interim/example/22/H2O/9_aMD_cMD/2000/0/3595ce0609206d95_Edih.dat"

# 	touch "workflow/data/interim/example/22/H2O/10_aMD_prod/2000/0/3595ce0609206d95_md.out"
# 	touch "workflow/data/interim/example/22/H2O/10_aMD_prod/2000/0/3595ce0609206d95_traj.netcdf"
# 	touch "workflow/data/interim/example/22/H2O/10_aMD_prod/2000/0/3595ce0609206d95_traj.ncdf"
# 	touch "workflow/data/interim/example/22/H2O/10_aMD_prod/2000/0/3595ce0609206d95_aMD.log"
# 	touch "workflow/data/interim/example/22/H2O/10_aMD_prod/2000/0/3595ce0609206d95_weights.dat"

# fi