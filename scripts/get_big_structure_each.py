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
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
sys.path.append('common')
# from common import get_target_file as gtf
# import measuring_volume.get_target_file as gtf
import common.get_target_file as gtf
import csv 
import glob
import networkx as nx

def get_nb_node(graph):
    """
    グラフの連結成分を調べ、最大の連結成分のノード数とそのノードリストを返す関数。

    :param graph: dict, グラフの隣接リスト
    :return: tuple, 最大の連結成分のノード数とそのノードリスト
    """
    # NetworkXのグラフを作成
    G = nx.Graph()

    # エッジを追加
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            G.add_edge(node, neighbor)

    # 連結成分を取得
    connected_components = list(nx.connected_components(G))

    # 各連結成分のサイズとノードリストを取得
    component_sizes = [len(component) for component in connected_components]
    max_size = max(component_sizes)
    max_component = max(connected_components, key=len)

    # 結果を表示
    print("連結成分ごとのノード数:", component_sizes)
    print("最大の連結成分のノード:", max_component)
    print("最大の連結成分のノード数:", max_size)

    return max_size, max_component

def make_initial_strands_data(target_dir):
    # particle2strands[particle_id] = strand_id
    particle2strand = {} 
    # strands2particle[strand_id] = {particle_ids}
    strands2particle = {} 
    top_f = open(gtf.get_top(target_dir))

    col = 0
    particle_id = 0
    initial_strand_num = -1

    for l in top_f:
        if col == 0:
            initial_strand_num = int(l.split(" ")[-1].rstrip('\n'))
        else:
            strand_id = int(l.split(" ")[0])
            if strand_id in strands2particle:
                strands2particle[strand_id].add(particle_id)
            else:
                strands2particle[strand_id] = {particle_id}
            particle2strand[particle_id] = strand_id
            particle_id += 1
        col += 1
    
    # print(initial_strand_num)
    
    top_f.close()
    return strands2particle, particle2strand

def get_bonds_data(bonds_name):
    bonds_f = open(bonds_name)
    line = 0

    header = ["id1", "id2", "FENE", "BEXC", "STCK", "NEXC", "HB", "CRSTCK", "CXSTCK", "total"]
    # ex. bonds_dic[((id1, id2), "FENE")] = value
    bonds_dic = {}

    for l in bonds_f:
        if line > 0 and l[0] != '#':
            lst = l.split(" ")
            lst[-1] = lst[-1][:-2]
            id1 = int(lst[0])
            id2 = int(lst[1])
            for index, h in enumerate(header[2:]):
                # print(id1, id2, h, float(lst[index + 2]))
                if lst[index + 2] == '-':
                    continue
                bonds_dic[((id1, id2), h)] = float(lst[index + 2])
        line += 1
    return bonds_dic

def get_connection_strands(bonds_name, strands2particle, particle2strand):


    new_strands2particle = dict(strands2particle)  # 新しい辞書を作成してコピー
    new_particle2strand = dict(particle2strand)  # 新しい辞書を作成してコピー
    
    bonds_dic = get_bonds_data(bonds_name)
    # HB < 0.0同士のparticle id1, id2である時、strand idが小さい方に合体させる
    
    for b in bonds_dic:
        
        if b[1] == "HB" and bonds_dic[b] < 0.0:
            particle_id1 = b[0][0]
            particle_id2 = b[0][1]
            # print("bonds : ", particle_id1, particle_id2)
            strand_id1 = new_particle2strand[particle_id1]
            strand_id2 = new_particle2strand[particle_id2]
            # strand idが異なる場合、小さい方に合わせる
            # if strand_id2 == strand_id1:
            if strand_id2 != strand_id1:
                print(strand_id2, strand_id1)
                s_id1 = min(strand_id1, strand_id2)
                s_id2 = max(strand_id1, strand_id2)
                if s_id2 in new_strands2particle:
                    # s_id2 → s_id1へ
                    for particle in list(new_strands2particle[s_id2]):
                        new_particle2strand[particle] = s_id1

                    new_strands2particle[s_id1] |= new_strands2particle[s_id2]
                    new_strands2particle.pop(s_id2)

    return new_strands2particle, new_particle2strand

def get_connection_strands2(bonds_name, strands2particle, particle2strand):
    
    bonds_dic = get_bonds_data(bonds_name)
    # HB < 0.0同士のparticle id1, id2である時、接続している
    # {1:[2], 2: [1], 10:[14, 20], ...} strand1とstrand2が接続している場合
    connected_strands = {}
    
    for b in bonds_dic:
        
        if b[1] == "HB" and bonds_dic[b] < 0.0:
            particle_id1 = b[0][0]
            particle_id2 = b[0][1]
            strand_id1 = particle2strand[particle_id1]
            strand_id2 = particle2strand[particle_id2]
            # strand idが異なる場合
            if strand_id2 != strand_id1:

                # print("Wow")
               
                if strand_id1 in connected_strands and not (strand_id2 in connected_strands[strand_id1]):
                    connected_strands[strand_id1].append(strand_id2)
                else:
                    connected_strands[strand_id1] = [strand_id2]

                if strand_id2 in connected_strands and not (strand_id1 in connected_strands[strand_id2]):
                    connected_strands[strand_id2].append(strand_id1)
                else:
                    connected_strands[strand_id2] = [strand_id1]
               
    return  connected_strands

def get_particle_strands_data(target_dir):
    strands2particle, particle2strand = make_initial_strands_data(target_dir)
    strands2particle, particle2strand = get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    return strands2particle, particle2strand 

def get_connected_strands_data(target_dir):
    strands2particle, particle2strand = make_initial_strands_data(target_dir)
    # bondsファイルを取得する
    connected_strands = get_connection_strands2(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    return connected_strands

def get_nb_seq_of_structure(dir):
    strands = get_connected_strands_data(dir)
    # {10: [26], 11: [10], 6: [17], 17: [13], 13: [18], 18: [13], 26: [10], 1: [27], 27: [1], 22: [28], 28: [22], 5: [28], 21: [29], 29: [21]}
    max_size, max_component = get_nb_node(strands)
    return max_size, max_component

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("usage : python3 get_big_structure_each.py type_of_l target")

    type_of_l = sys.argv[1]
    target = sys.argv[2]
    base_dir = "/home/user/SA-EDS/"

    dirs = glob.glob(base_dir + target + "/" + type_of_l + "*")

    file_path=f"home/user/SA-EDS/dataset/{type_of_l}_data_{target}_each.pkl"

    f = open(file_path, "rb")
    data = pickle.load(f)

    for dir10 in dirs:
        dir = glob.glob(dir10 + "/*")

        for dir1 in dir:
            print("dir1", dir1)
            temp = dir1.split("/")[-1].split("_")[-2]

            try:
                get_nb_seq_of_structure(dir1)
                max_size, max_component = get_nb_seq_of_structure(dir1)                
                data[(dir1.split("/")[-1], temp)]["max_nb_strand"] = max_size
                data[(dir1.split("/")[-1], temp)]["max_nb_strand_component"] = max_component
            except Exception as e:
                print(f"Error in get_nb_seq_of_structure for {dir1}: {e}")
        
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)
