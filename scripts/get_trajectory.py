
import sys
sys.path.append("../common/")
sys.path.append("../measuring_volume/")
import measuring_volume.run_output_bonds as rob
import measuring_volume.convexhull_volume2 as cv
# import measuring_volume.get_top_data as gtd
import common.get_top_data as gtd
import sys
sys.path.append('../')
sys.path.append('.')
sys.path.append('measuring_volume/')
# from common import get_target_file as gtf
# from common import check_dir as cd
import common.get_target_file as gtf
import common.check_dir as cd
import subprocess
import os.path as osp
import sys


if __name__ == "__main__": 
    if len(sys.argv) != 3:
        #     subprocess.call(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/measuring_volume/get_trajectory.py", type_of_l, path])
        print("usage : python3 get_trajectory.py L3 target")
    target_dir=sys.argv[2]
    print("target_dir", target_dir, sys.argv)

    cd.new_input(target_dir)
    top = gtf.get_top(target_dir)
    conf = gtf.get_conf(target_dir)
    input = gtf.get_new_input(target_dir)
    traj = gtf.get_traj(target_dir)

    # "input/results/L1-GA100000-0.80-ERT-0_277_11/trajectory_L1-GA100000-0.80-ERT-0_277_11.dat"
    # export trajectory
    # path="/home/user/SA-EDS/scripts/measuring_volume/read_trajectory_simple.py"
    # result_path=target_dir
    # subprocess.call(["/home/user/venv/bin/python3", path, traj, result_path])
    # print(target_dir)
    basename = osp.splitext(osp.basename(traj))[0]
    path="/home/user/SA-EDS/scripts/measuring_volume/output_bonds_traj.py"
    for i in range(10):
        new_traj = target_dir + basename + "_" + str(i)
        index = str(i)
        subprocess.call(["/home/user/venv/bin/python3", path, input, new_traj, top, index])
