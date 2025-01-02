import sys
sys.path.append("analyze")
sys.path.append("measuring_volume")
sys.path.append("common")
import pickle
import glob
import analyze.count_strands as cs
# import convexhull_volume as cv
import measuring_volume.convexhull_volume2 as cv
import measuring_volume.run_output_bonds as rob
import common.check_dir as cd
import common.get_target_file as gtf
import datetime
from classify_seq import make_input_seq as mis
import get_x2
import load_results as lr
import time

def make_data_dir_initial(path, type_of_l, target):
    # input/results/fromQD/r20230613134109/r1686027494794-20/*
    dirs = glob.glob(path + "/*")
    dic = {}
    # y作る
    for index, d in enumerate(dirs):
        print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
        if cd.included_full_files(d) == False:
            continue
        # print(d)
        # make bonds
        rob.make_bonds(d)
        before, after = cs.count_strands(d)
        meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
        key = (d.split("/")[-1], d.split("_")[-2])
        dic[key] = {}
        dic[key]["ratio_of_volume"] = after/before
        dic[key]["mean_volume"] = meanv
        dic[key]["deviation_of_volume"] = devv
        print(d.split("/")[-1])        
    return dic

def extract_string_between_last_two_underscores(path):
    # 最後から1つ目の_と最後から2つ目の_に挟まれた文字列を抽出
    last_underscore_index = path.rfind('_')  # 最後から1つ目の_のインデックスを取得
    second_last_underscore_index = path.rfind('_', 0, last_underscore_index)  # 最後から2つ目の_のインデックスを取得
    if last_underscore_index != -1 and second_last_underscore_index != -1:
        return path[second_last_underscore_index + 1:last_underscore_index]
    else:
        return None
    
def make_initial():
    if len(sys.argv) != 3:
        print("usage : python3 make_data.py [type_of_path] [target]")
        print("ex)")
        print("type of l : L1")
        print("version : initial/")

    print(sys.argv)
    
    type_of_l=sys.argv[1]
    target=sys.argv[2]

    dirs = glob.glob("home/user/SA-EDS/" + target + "/" + type_of_l + "*")
    # dirs = glob.glob(target + "/" + type_of_l + "*")
    print(dirs)
    initial_dic = {}
    for dir in dirs:
        print(dir)
        path=dir
        dic = make_data_dir_initial(path, type_of_l, target)

        print(dic)
        initial_dic.update(dic)

    result_path =  "home/user/SA-EDS/dataset/" + type_of_l + "_data_" + target + "_each.pkl"
    print(result_path)

    with open(result_path, "wb") as tf:
        pickle.dump(initial_dic,tf)


if __name__ == '__main__':
    make_initial()
    # make_QD()
