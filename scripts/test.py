#! /bin/python

import os.path as osp
import common.get_target_file as gtf
import measuring_volume.get_top_data as gtd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.optimize import curve_fit
import sys
import pickle
import subprocess

def breakdown_trajectory_file(filename, path="./"):
    with open(filename, "r") as file:
        lines = file.readlines()
    i = 0
    expi = 0
    basename = osp.splitext(osp.basename(filename))[0]
    to_export = lines[0]
    size = len(lines)
    while i < size:
        i += 1
        if i >= size or "t =" in lines[i]:
            with open(osp.join(path, basename+"_"+str(expi)), "w") as expfile:
                expfile.write(to_export)
                expi += 1
                if i < size:
                    to_export = lines[i]
        else:
            to_export += lines[i]

def readbonds_str2dic(line):
    key_lst = ["id1", "id2", "FENE", "BEXC", "STCK", "NEXC", "HB", "CRSTCK", "CXSTCK", "total"]
    new_dic = {}
    for i, str_v in enumerate(line.split(" ")):
        if key_lst[i] == "id1" or key_lst[i] == "id2":
            new_dic[key_lst[i]] = int(str_v.replace("\n", ""))
        else:
            new_dic[key_lst[i]] = float(str_v.replace("\n", ""))
    return new_dic

def HB_connectivity(path, type_of_l):
    id_pair_dic = {i: set() for i in range(10)}
    for i in range(10):
        target = 9 - i
        filename = path + "bonds_" + str(target)
        with open(filename, "r") as file:
            lines = file.readlines()
        for l in lines:
            if "#" not in l:
                bond_dic = readbonds_str2dic(l)
                if bond_dic["HB"] < 0:
                    id_pair_dic[target].add((bond_dic["id1"], bond_dic["id2"]))
    return id_pair_dic

def count_domain(target_dir):
    req_dir = gtf.get_req(target_dir)
    with open(req_dir, 'r') as f:
        nb_domain = 0
        for l in f:
            if l[0] == 's':
                nb_domain += (l.split(" ").index("@") - l.split(" ").index("=") - 1) * 5
    return nb_domain

def fit_sigmoid(path, type_of_l):
    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path, type_of_l)
    trans_dic = {}
    int_step_nb = 100000
    nb_domain = count_domain(path)
    for i in range(9):
        step_nb = (i + 1) * int_step_nb
        trans_dic[step_nb] = {}
        trans_dic[step_nb]["Dissociation"] = len(id_pair_dic[i] - id_pair_dic[i + 1]) / nb_domain
        trans_dic[step_nb]["Binding"] = len(id_pair_dic[i + 1] - id_pair_dic[i]) / nb_domain
        trans_dic[step_nb]["Stability"] = len(id_pair_dic[i] & id_pair_dic[i + 1]) / nb_domain

    xdata = np.linspace(int_step_nb, 9 * int_step_nb, 9)
    int_xdata = [int(x) for x in xdata]
    y = [trans_dic[i]["Stability"] for i in int_xdata]
    popt, pcov = curve_fit(sigmoid, xdata, y)
    print(popt)
    plt.plot(xdata, y, 'b-', label='data')
    plt.plot(xdata, sigmoid(xdata, *popt), 'r-', label='sigmoid')
    plt.xlabel('Steps')
    plt.ylabel('Stability')
    plt.savefig(f"{path}curv_fit.png")
    plt.show()

def strand_connectivity(path, type_of_l):
    strands2particle, particle2strand = gtd.make_initial_strands_data(path)
    id_pair_dic = HB_connectivity(path, type_of_l)
    nb_domain = count_domain(path)
    for i in range(9):
        print((i + 1) * 100000, "\tsteps")
        print("解離", len(id_pair_dic[i] - id_pair_dic[i + 1]), "/", nb_domain)
        print("結合", len(id_pair_dic[i + 1] - id_pair_dic[i]), "/", nb_domain)
        print("安定", len(id_pair_dic[i] & id_pair_dic[i + 1]), "/", nb_domain)
    for i in range(9):
        G = nx.Graph()
        for strand in strands2particle:
            G.add_node(str(strand))
        for pair in id_pair_dic[i]:
            G.add_edge(str(particle2strand[pair[0]]), str(particle2strand[pair[1]]))
        pos = nx.circular_layout(G)
        nx.draw(G, pos, with_labels=True, node_size=500, font_size=16, font_color='white')
        plt.title(f"{i + 1}*100000 step")
        plt.savefig(f"{path}stable_connection.png")
        plt.show()

def sigmoid(x, L, x0, k, b):
    return L / (1 + np.exp(-k * (x - x0))) + b

def calc_curv_fit(filename, extra_path, type_of_l):
    breakdown_trajectory_file(filename, path=extra_path)
    subprocess.run(["/home/user/venv/bin/python3", "/home/user/SA-EDS/scripts/get_trajectory.py", type_of_l, extra_path])
    strand_connectivity(path=extra_path, type_of_l=type_of_l)
    fit_sigmoid(path=extra_path, type_of_l=type_of_l)

if __name__ == '__main__':
    base_dir = "/home/user/SA-EDS/"
    filename = "int_initial/L3_initial_0/L3-GA100000-0.50-ERT-0_277_2/trajectory_L3-GA100000-0.50-ERT-0_277_2.dat"
    extra_path = "int_initial/L3_initial_0/L3-GA100000-0.50-ERT-0_277_2/"
    filename = base_dir + filename
    extra_path= base_dir + extra_path
    type_of_l = "L3"
    calc_curv_fit(filename, extra_path, type_of_l)
