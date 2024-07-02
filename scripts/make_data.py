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


# def make_data(path, type_of_l, version=4):
#     # input/results/oxdna_random_6/L1/d-0-1/L1_d-0-1_2023-01-27-083608/L1_d-0-1_2023-01-27-083608/bonds
#     dirs = glob.glob(path + "/" + type_of_l + "/*/*/*/")
#     dic = {}
#     for index, d in enumerate(dirs):
#         print(d)
#         print(datetime.datetime.now(), " : ", index, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         print(dic[d])
    
#     result_path =  "../data/dic/" + type_of_l + "_data_" + str(version) + ".pkl"
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic


# def make_data_fromQD(path, type_of_l, version=4):
#     # input/results/fromQD/r20230613134109/r1686027494794-20/*
#     dirs = glob.glob(path + "*/")
#     dic = {}
#     # y作る
#     for index, d in enumerate(dirs):
#         print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         print(d)
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         print(dic[d])
    
#     result_path =  "../data/dic/" + type_of_l + "_data_fromKED_" + version + ".pkl"
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic

# def make_data_one(version=4):
#     # '../input/results/fromQD/big/r1692694568274-908
#     type_of_l="L3"
#     dirs = ["../input/results/fromQD/big/r1692694568274-908"]
#     # dirs = glob.glob(path + "*")
#     dic = {}
#     # y作る
#     for index, d in enumerate(dirs):
#         print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         print(d)
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         print(dic[d])
    
#     result_path =  "../data/dic/" + type_of_l + "_data_r1692694568274-908.pkl"
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic

# def make_data_fromKED(path, type_of_l, target, temperature):
#     # input/results/fromQD/r20230613134109/r1686027494794-20/*
#     dirs = glob.glob(path + "/*")
#     dic = {}
#     print("dirs ", dirs)
#     # y作る
#     for index, d in enumerate(dirs):
#         print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         # print(d)
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         # print(dic[d])
    
#     result_path =  "../data/dic/" + type_of_l + "_data_fromKED_" + target + "_t" + temperature + ".pkl"
#     print(result_path)
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic
# def make_data_fromKED(path, type_of_l, target, temperature):
#     # input/results/fromQD/r20230613134109/r1686027494794-20/*
#     dirs = glob.glob(path + "/*")
#     dic = {}
#     print("dirs ", dirs)
#     # y作る
#     for index, d in enumerate(dirs):
#         print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         # print(d)
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         # print(dic[d])
    
#     result_path =  "../data/dic/" + type_of_l + "_data_fromKED_" + target + "_t" + temperature + ".pkl"
#     print(result_path)
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic

# def make_data_dir(path, type_of_l, target):
#     # input/results/fromQD/r20230613134109/r1686027494794-20/*
#     dirs = glob.glob(path + "/*")
#     dic = {}
#     print("dirs ", dirs)
#     # y作る
#     for index, d in enumerate(dirs):
#         print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
#         if cd.included_full_files(d) == False:
#             continue
#         # print(d)
#         # make bonds
#         rob.make_bonds(d)
#         before, after = cs.count_strands(d)
#         meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
#         dic[d] = {}
#         dic[d]["ratio_of_volume"] = after/before
#         dic[d]["mean_volume"] = meanv
#         dic[d]["deviation_of_volume"] = devv
        
#         # print(before, after)
#         # print(dic[d])
    
#     result_path =  "../data/dic/data_" + target.split("/")[-1] + ".pkl"
#     print(result_path)
    
#     with open(result_path, "wb") as tf:
#         pickle.dump(dic,tf)

#     return dic

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
        dic[d] = {}
        dic[d]["ratio_of_volume"] = after/before
        dic[d]["mean_volume"] = meanv
        dic[d]["deviation_of_volume"] = devv
        
        # print(before, after)
        # print(dic[d])
        
    
    return dic


# def test():
#     path="../input/results/oxdna_random_6"
#     type_of_l="L1"
#     make_data(path, type_of_l)

# def make_fromQD():
#     # path="../input/results/fromQD/r20230613134109/"
#     if len(sys.argv) != 3:
#         print("usage : python3 make_data.py [type_of_path] [version]")
#         print("ex)")
#         print("type of l : L2")
#         print("version : r20230613134109")
    
#     path="../input/results/fromQD/" + sys.argv[2] + "/"
#     type_of_l=sys.argv[1]
#     version=sys.argv[2]
#     print(path, type_of_l, version)
#     make_data_fromKED(path, type_of_l, version)

# def make_fromKED():
#     # path="../input/results/fromKED/L3/I_t308"
#     if len(sys.argv) != 4:
#         print("usage : python3 make_data.py [type_of_path] [target] [temperature]")
#         print("ex)")
#         print("type of l : L3")
#         print("version : I")
#         print("version : 308")
    
#     type_of_l=sys.argv[1]
#     target=sys.argv[2]
#     temperature=sys.argv[3]

#     path="../input/results/fromKED/" + type_of_l + "/" + target + "_t" + temperature + ""
#     print(path, type_of_l, target, temperature)
#     make_data_fromKED(path, type_of_l, target, temperature)

# def make_dir():
#     # path="../input/results/fromNPQD/L1_t277"
#     if len(sys.argv) != 4:
#         print("usage : python3 make_data.py [type_of_path] [target]")
#         print("ex)")
#         print("type of l : L1")
#         print("version : ../input/results/fromNPQD/L1_t277")
    
#     type_of_l=sys.argv[1]
#     target=sys.argv[2]

#     path=target
#     make_data_dir(path, type_of_l, target)

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
        tmp_dic = {}
        for key in dic:
            # 温度ごとに分ける
            temp = extract_string_between_last_two_underscores(key)
            if temp in tmp_dic:
                tmp_dic[temp].append(dic[key])
            else:
                tmp_dic[temp] = []
                tmp_dic[temp].append(dic[key])

        # 平均を作る
        for temp in tmp_dic:
            data = tmp_dic[temp]
            mean_ratio_of_volume = sum(item['ratio_of_volume'] for item in data) / len(data)
            mean_mean_volume = sum(item['mean_volume'] for item in data) / len(data)
            mean_deviation_of_volume = sum(item['deviation_of_volume'] for item in data) / len(data)
            initial_dic[(temp, dir)] = {'ratio_of_volume': mean_ratio_of_volume , 'mean_volume': mean_mean_volume , 'deviation_of_volume': mean_deviation_of_volume }
        # {(temp, dir) : {'ratio_of_volume': ... , 'mean_volume': ..., 'deviation_of_volume': ...}}
        # initial_dic[( , dir)] = {}

        result_path =  "home/user/SA-EDS/dataset/" + type_of_l + "_data_" + target + ".pkl"
        print(result_path)

        with open(result_path, "wb") as tf:
            pickle.dump(initial_dic,tf)

def make_QD():
    if len(sys.argv) != 4:
        print("usage : python3 make_data.py [type_of_path] [target]")
        print("ex)")
        print("type of l : L1")
        print("version : ../input/results/QD_1/L1_1/")
    
    print(sys.argv)
    
    type_of_l=sys.argv[1]
    target=sys.argv[2]

    dirs = glob.glob(target + "/" + "*")
    initial_dic = {}
    for dir in dirs:
        path=dir
        dic = make_data_dir_initial(path, type_of_l, target)
        tmp_dic = {}
        for key in dic:
            # 温度ごとに分ける
            temp = extract_string_between_last_two_underscores(key)
            if temp in tmp_dic:
                tmp_dic[temp].append(dic[key])
            else:
                tmp_dic[temp] = []
                tmp_dic[temp].append(dic[key])

        # 平均を作る
        for temp in tmp_dic:
            data = tmp_dic[temp]
            mean_ratio_of_volume = sum(item['ratio_of_volume'] for item in data) / len(data)
            mean_mean_volume = sum(item['mean_volume'] for item in data) / len(data)
            mean_deviation_of_volume = sum(item['deviation_of_volume'] for item in data) / len(data)
            initial_dic[(temp, dir)] = {'ratio_of_volume': mean_ratio_of_volume , 'mean_volume': mean_mean_volume , 'deviation_of_volume': mean_deviation_of_volume }
        # {(temp, dir) : {'ratio_of_volume': ... , 'mean_volume': ..., 'deviation_of_volume': ...}}
        # initial_dic[( , dir)] = {}
        print(initial_dic)

        target=sys.argv[2].split("/")[3]
        print(target)
        result_path =  "../data/dic/" + type_of_l + "_data_" + target + ".pkl"

        with open(result_path, "wb") as tf:
            pickle.dump(initial_dic,tf)
    

if __name__ == '__main__':
    make_initial()
    # make_QD()
