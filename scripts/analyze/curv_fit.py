#! /bin/python

import os.path as osp
import common.get_target_file as gtf
import measuring_volume.get_top_data as gtd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import pickle
import subprocess
import glob

def breakdown_trajectory_file(filename, path="./"):
    # print(filename, path)
    with open(filename, "r") as file:
        lines = file.readlines()
    i = 0
    expi = 0
    basename = osp.splitext(osp.basename(filename))[0]
    to_export = lines[0]
    size = len(lines)
    while i < size:
        i+= 1
        if i >= size or "t =" in lines[i]:
            with open(osp.join(path,basename+"_"+str(expi)),"w") as expfile:
                expfile.write(to_export)
                expi +=1
                if i <size:
                    to_export = lines[i]
        else:
            to_export += lines[i]

def readbonds_str2dic(line):
    # ['696', '697', '0.0128922', '0', '-0', '0', '0', '0', '0', '0.0128922\n']
    key_lst = ["id1", "id2", "FENE", "BEXC", "STCK", "NEXC", "HB", "CRSTCK", "CXSTCK", "total"]
    new_dic = {}
    for i, str_v in enumerate(line.split(" ")):
        if key_lst[i] == "id1" or key_lst[i] == "id2":
            new_dic[key_lst[i]] = int(str_v.replace("\n",""))
        else:
            new_dic[key_lst[i]] = float(str_v.replace("\n",""))
    return new_dic

def HB_connectivity(path, type_of_l):
    # bonds_0 - bonds_9について、HBで結合している個数を出力する
    id_pair_dic = {
        9 : set(),
        8 : set(),
        7 : set(),
        6 : set(),
        5 : set(),
        4 : set(),
        3 : set(),
        2 : set(),
        1 : set(),
        0 : set()
       
    }


    for i in range(10):
        target = 9 - i
        filename = path + "bonds_" + str(target)
        with open(filename, "r") as file:
            lines = file.readlines()
        # id1 id2 FENE BEXC STCK NEXC HB CRSTCK CXSTCK total
        for l in lines:
            if "#" not in l:
                bond_dic = readbonds_str2dic(l)
                if bond_dic["HB"] < 0:
                    id_pair_dic[target].add((bond_dic["id1"], bond_dic["id2"]))

    return id_pair_dic

def count_domain(target_dir):
    req_dir = gtf.get_req(target_dir)
    f = open(req_dir, 'r')
    nb_domain = 0
    for l in f:
        if l[0] == 's':
            nb_domain = nb_domain + (l.split(" ").index("@") - l.split(" ").index("=") - 1)*5

    return nb_domain

def stability(path, type_of_l):
    id_pair_dic = HB_connectivity(path, type_of_l)

    # trans_dic[step] = {"disociation" : 4, "binding" : 3, "stability" : 2}
    trans_dic = {}
    int_step_nb = 100000

    nb_domain = count_domain(path)


    for i in range(9):
        # 解離: Dissociation
        # 結合: Binding
        # 安定: Stability
        # domainの数で割る
        step_nb = (i + 1)*int_step_nb
        trans_dic[step_nb] = {}
        trans_dic[step_nb]["Dissociation"] = len(id_pair_dic[i] - id_pair_dic[i + 1]) / nb_domain
        trans_dic[step_nb]["Binding"] = len(id_pair_dic[i + 1] - id_pair_dic[i]) / nb_domain
        trans_dic[step_nb]["Stability"] = len(id_pair_dic[i] & id_pair_dic[i + 1]) / nb_domain

    xdata = np.linspace(int_step_nb, 9*int_step_nb, 9)
    int_xdata = [int(x) for x in xdata]
    # print(int_xdata, trans_dic.keys())
    y = [trans_dic[i]["Stability"] for i in int_xdata]
    return y


def fit_sigmoid(path, type_of_l):

    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path, type_of_l)

    # trans_dic[step] = {"disociation" : 4, "binding" : 3, "stability" : 2}
    trans_dic = {}
    int_step_nb = 100000

    nb_domain = count_domain(path)


    for i in range(9):
        # 解離: Dissociation
        # 結合: Binding
        # 安定: Stability
        # domainの数で割る
        step_nb = (i + 1)*int_step_nb
        trans_dic[step_nb] = {}
        trans_dic[step_nb]["Dissociation"] = len(id_pair_dic[i] - id_pair_dic[i + 1]) / nb_domain
        trans_dic[step_nb]["Binding"] = len(id_pair_dic[i + 1] - id_pair_dic[i]) / nb_domain
        trans_dic[step_nb]["Stability"] = len(id_pair_dic[i] & id_pair_dic[i + 1]) / nb_domain

    xdata = np.linspace(int_step_nb, 9*int_step_nb, 9)
    int_xdata = [int(x) for x in xdata]
    # print(int_xdata, trans_dic.keys())
    y = [trans_dic[i]["Stability"] for i in int_xdata]
    xdata = np.float64([float(x/int_step_nb) for x in xdata])
    try :
        popt, pcov = curve_fit(func, xdata, y, maxfev = 100)
    except RuntimeError:
        popt = [0, 0, 1]
    plt.figure()
    plt.plot(xdata, y, 'b-', label='Stably Bonded Base Pair Count / Domain Count')
    plt.plot(xdata, func(xdata, *popt), 'r-', label='Sigmoid Curve Approximation')
    # Add labels and a legend
    plt.xlabel('Steps (scaled)')
    plt.ylabel('Stability')
    plt.legend()
    plt.savefig(f"{path}curv_fit.png")  # グラフを保存
    plt.show()
    plt.close()

    return popt
    


def strand_connectivity(path, type_of_l):
    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path, type_of_l)

    # trans_dic[step] = {"disociation" : 4, "binding" : 3, "stability" : 2}
    trans_dic = {}
    nb_domain = count_domain(path)

    # for i in range(9):
        # 解離: Dissociation
        # 結合: Binding
        # 安定: Stability
        # domainの数で割る
        # print((i + 1)*20000, "\tsteps")
        # print("解離", len(id_pair_dic[i] - id_pair_dic[i + 1]), "/", len(particle2strand))
        # print("結合", len(id_pair_dic[i + 1] - id_pair_dic[i]), "/", len(particle2strand))
        # print("安定", len(id_pair_dic[i] & id_pair_dic[i + 1]), "/", len(particle2strand))
        # print("解離", len(id_pair_dic[i] - id_pair_dic[i + 1]), "/", nb_domain)
        # print("結合", len(id_pair_dic[i + 1] - id_pair_dic[i]), "/", nb_domain)
        # print("安定", len(id_pair_dic[i] & id_pair_dic[i + 1]), "/", nb_domain)

    # for i in range(9):
    #     G = nx.Graph()
    #     for strand in strands2particle:
    #         G.add_node(str(strand))
    #     for pair in id_pair_dic[i]:
    #         G.add_edge(str(particle2strand[pair[0]]), str(particle2strand[pair[1]]))
    #     # グラフの描画
    #     # pos = nx.spring_layout(G)  # レイアウトの計算
    #     pos = nx.circular_layout(G)
    #     nx.draw(G, pos, with_labels=True, node_size=500, font_size=16, font_color='white')

        # # 新しい図を作成
        # plt.figure()
        
        # # グラフの描画
        # pos = nx.circular_layout(G)
        # nx.draw(G, pos, with_labels=True, node_size=500, font_size=16, font_color='white')

        # # グラフのタイトルと保存
        # plt.title(f"{i + 1}*100000 step")
        # plt.savefig(f"{path}stable_connection_{(i + 1) * 100000}.png")  # ステップ数を含んだファイル名で保存
        # plt.close()  # 図を閉じる

def func(x, a, b, c):
    return a /(np.exp(-b * x) + c)

def calc_curv_fit(filename, extra_path, type_of_l):
    file_dic = gtf.file_dic(extra_path)
    if "trajectory1" not in file_dic.keys():
        breakdown_trajectory_file(filename, path=extra_path)
        subprocess.run(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/get_trajectory.py", type_of_l, extra_path])
    strand_connectivity(path=extra_path, type_of_l=type_of_l)
    print(extra_path)
    try:
        [a, b, c] = fit_sigmoid(path=extra_path, type_of_l=type_of_l)
    except RuntimeError:
        [a, b, c] = [0, 0, 0]

    return a, b, c

# def main():
#     if len(sys.argv) != 3:
#         print("usage : python3 connectivity L{1-3} target")
#     type_of_l = sys.argv[1]
#     target = sys.argv[2]
#     file_path=f"home/user/SA-EDS/dataset/x_{target}_{type_of_l}.pkl"
#     f = open(file_path, "rb")
#     data = pickle.load(f)
#     for (temp, dir) in data:
#         print(temp, dir)

#     with open(file_path, 'wb') as file:
#         pickle.dump(data, file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("usage : python3 curv_fig.py type_of_l target")

    # target = "int_initial"
    # type_of_l = "L3"
    type_of_l = sys.argv[1]
    target = sys.argv[2]
    base_dir = "/home/user/SA-EDS/"

    dirs = glob.glob(base_dir + target + "/" + type_of_l + "*")

    file_path=f"home/user/SA-EDS/dataset/{type_of_l}_data_{target}.pkl"
    f = open(file_path, "rb")
    data = pickle.load(f)

    temp_lst = ["277", "298", "308", "318", "328", "338", "348", "358"]

    for dir10 in dirs:
        dir = glob.glob(dir10 + "/*")
        print(dir10, dir)
        temp_dic = {}

        for temp in temp_lst:
            temp_dic[temp] = {"a" : 0, "b" : 0, "c" : 0, "nb" : 0}

        for dir1 in dir:
            try:
                filename = gtf.get_traj(dir1)
                extra_path = dir1 + "/"
                tmp_temp = extra_path.split("_")[-2]
                if tmp_temp == "348" or tmp_temp == "358":
                    a = 0
                    b = 0
                    c = 1
                else:
                    calc_curv_fit(filename, extra_path, type_of_l)
                    a, b, c = calc_curv_fit(filename, extra_path, type_of_l)
                temp_dic[tmp_temp]["a"] += a
                temp_dic[tmp_temp]["b"] += b
                temp_dic[tmp_temp]["c"] += c
                temp_dic[tmp_temp]["nb"] += 1
                print(a, b, c)
            except Exception as e:
                print(f"error: {e}")

            

        for temp in temp_lst:
            if temp_dic[temp]["nb"] != 0:
                data[(temp, dir10[1:])]["sigmoid"] = {"a" : temp_dic[temp]["a"] / temp_dic[temp]["nb"], "b" : temp_dic[temp]["b"] / temp_dic[temp]["nb"], "c" : temp_dic[temp]["c"] / temp_dic[temp]["nb"]}

    with open(file_path, 'wb') as file:
        pickle.dump(data, file)




    # filename = "int_initial/L3_initial_0/L3-GA100000-0.50-ERT-0_277_2/trajectory_L3-GA100000-0.50-ERT-0_277_2.dat"
    # extra_path = "int_initial/L3_initial_0/L3-GA100000-0.50-ERT-0_277_2/"
    # filename = base_dir + filename
    # extra_path= base_dir + extra_path
    # type_of_l = "L3"
    # print(sys.argv)
    # calc_curv_fit(filename, extra_path, type_of_l)

