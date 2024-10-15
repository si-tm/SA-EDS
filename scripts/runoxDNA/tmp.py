
if __name__ == '__main__':
    for i in range(8):
        print(f"cd /clusterhome/mhyakuzuka/oxdna/L3_1/knight{i}")
        print(f"nohup ./runDocker.sh &")