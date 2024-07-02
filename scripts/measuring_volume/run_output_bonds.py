import sys
sys.path.append('../')
sys.path.append('.')
sys.path.append('measuring_volume/')
import common.get_target_file as gtf
import common.check_dir as cd
import subprocess

def make_bonds(target_dir):
    cd.new_input(target_dir)
    top = gtf.get_top(target_dir)
    conf = gtf.get_conf(target_dir)
    input = gtf.get_new_input(target_dir)

    path="/home/user/SA-EDS/scripts/measuring_volume/output_bonds.py"
    
    # subprocess.call(["ls", "-l", "/home/user/SA-EDS/oxDNA/build/bin/" ,"|", "grep", "DNAanalysis"])
    subprocess.call(["/home/user/venv/bin/python3", path, input, conf, top])

def main():
    make_bonds("../input/results/oxdna_ked/seqA/A5/test_a5_200000_1")

if __name__ == '__main__':
    main()
