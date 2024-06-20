
import subprocess
import sys

if __name__ == '__main__':

    type_of_l = sys.argv[1]

    subprocess.call(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/get_x2.py", type_of_l, "initial"])
    subprocess.call(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/make_data.py", type_of_l, "initial"])
    subprocess.call(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/connectivity.py", type_of_l])
