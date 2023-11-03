#!/bin/bash
# ------------------------------------------------------------
# This script runs all Serious Game DRIFTERS test cases
#
# (D. Bruciaferri)
# ------------------------------------------------------------

echo "SERIOUS GAME DRIFTERS TEST CASE STARTED"

OUT_DIR_EUL="../OUTPUT/EULER_SG_drifter_test_case"
OUT_DIR_RUN="../OUTPUT/RUNGE_SG_drifter_test_case"

./euler_RUN_drifter_TESTCASE.sh

./runge_RUN_drifter_TESTCASE.sh

python ../PLOT/plot_traj_BEACON.py
python ../PLOT/plot_traj_CODE.py
python ../PLOT/plot_traj_ISPHERE.py
python ../PLOT/plot_traj_SLDMB.py

mv eul*_trajectory.png ${OUT_DIR_EUL}
mv run*_trajectory.png ${OUT_DIR_RUN}

echo ""
echo "SERIOUS GAME DRIFTERS TEST CASE COMPLETED"
echo ""
