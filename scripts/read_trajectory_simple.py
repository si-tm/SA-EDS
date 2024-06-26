#! /bin/python

import os.path as osp
import common.get_target_file as gtf
import measuring_volume.get_top_data as gtd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def breakdown_trajectory_file(filename, path="./"):
    print(filename, path)
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

def HB_connectivity(path):
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


def fit_sigmoid(path):
    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path)

    # trans_dic[step] = {"disociation" : 4, "binding" : 3, "stability" : 2}
    trans_dic = {}
    int_step_nb = 20000

    for i in range(9):
        # 解離: Dissociation
        # 結合: Binding
        # 安定: Stability
        # domainの数で割る
        step_nb = (i + 1)*int_step_nb
        trans_dic[step_nb] = {}
        trans_dic[step_nb]["Dissociation"] = len(id_pair_dic[i] - id_pair_dic[i + 1])
        trans_dic[step_nb]["Binding"] = len(id_pair_dic[i + 1] - id_pair_dic[i])
        trans_dic[step_nb]["Stability"] = len(id_pair_dic[i] & id_pair_dic[i + 1])

    xdata = np.linspace(int_step_nb, 10*int_step_nb, 10)
    int_xdata = [int(x) for x in xdata]
    print(int_xdata)
    y = [trans_dic[i]["Stability"] for i in int_xdata]
    plt.plot(xdata, y, 'b-', label='data')
    plt.show()
    
    # popt, pcov = curve_fit(func, xdata, ydata)


def strand_connectivity(path):
    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path)

    # trans_dic[step] = {"disociation" : 4, "binding" : 3, "stability" : 2}
    trans_dic = {}

    for i in range(9):
        # 解離: Dissociation
        # 結合: Binding
        # 安定: Stability
        # domainの数で割る
        print((i + 1)*20000, "\tsteps")
        print("解離", len(id_pair_dic[i] - id_pair_dic[i + 1]), "/", len(particle2strand))
        print("結合", len(id_pair_dic[i + 1] - id_pair_dic[i]), "/", len(particle2strand))
        print("安定", len(id_pair_dic[i] & id_pair_dic[i + 1]), "/", len(particle2strand))

    for i in range(9):
        G = nx.Graph()
        for strand in strands2particle:
            G.add_node(str(strand))
        for pair in id_pair_dic[i]:
            G.add_edge(str(particle2strand[pair[0]]), str(particle2strand[pair[1]]))
        # グラフの描画
        # pos = nx.spring_layout(G)  # レイアウトの計算
        pos = nx.circular_layout(G)
        nx.draw(G, pos, with_labels=True, node_size=500, font_size=16, font_color='white')

        # グラフの表示
        plt.title(f"{i + 1}*20000 step")
        plt.show()

def func(x, a, b, c):
    return a /(np.exp(-b * x) + c)


if __name__ == "__main__":
    import sys
    extra_path = "./"
    if len(sys.argv) > 2:
        extra_path = sys.argv[2]
    # breakdown_trajectory_file(sys.argv[1], path=extra_path)
    strand_connectivity(path=extra_path)
    # HB_connectivity(path=extra_path)
    # fit_sigmoid(path=extra_path)