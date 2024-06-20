import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
from common import get_target_file as gtf
from common import check_seq as cs
from common import check_dir as cd
from classify_seq import make_input_seq as mis
from classify_seq import get_structure_from_pil as gsfp
from classify_seq import get_structure_from_req as gsfq
import load_results as lr
import csv
import pickle
import glob

def get_target_x(target_dir):
    seq_file_name = gtf.get_seq(target_dir)
    seq = seq_file_name.split("/")[-1][3]

    seq_lst = []

    if cd.is_random(target_dir):
        seq_lst = gsfq.seq2structure(target_dir)
    if cd.is_fromQD(target_dir):
        seq_lst = gsfq.seq2structure(target_dir)
    if cd.is_fromKED(target_dir):
        seq_lst = gsfq.seq2structure(target_dir)
    if cd.is_fromKED(target_dir):
        seq_lst = gsfq.seq2structure(target_dir)
    else:
        # seq_lst = gsfq.seq2structure(target_dir)
        # seq_lst = gsfp.seq2structure(seq)
        seq_lst = gsfq.seq2structure(target_dir)

    # print(seq_lst)
    return seq_lst

def get_x_dic(csv_path, dirs):
    
    x_dic = {}

    # print(csv_path)

    x_dic_l1 = mis.seq_dic(csv_path)
    x_dic_l1_sub = x_dic_l1

    # print(x_dic_l1)
    # print(dirs)


    for d in dirs:
        # print("dir : ", d)
        structurs = get_target_x(d)
        for str in structurs:
            # print(cs.lst2str(str))
            x_dic_l1_sub[cs.lst2str(str)] = 1
        x_dic[d] = x_dic_l1_sub
        x_dic_l1_sub =mis.seq_dic(csv_path)
    
    return x_dic

def get_x_domain_dic(dirs):

    dic = {}
    
    for dir in dirs:
        x_domain_dic = gsfq.seq2domain(dir)
        dic[dir] = x_domain_dic

    return dic

def make_x_dir_initial():
    if len(sys.argv) != 3:
        print("usage : python3 get_x2.py <type of l> <target>")
        print("ex)")
        print("type of l : L1")
        print("target : initial")

    dirs = glob.glob("/home/user/SA-EDS/" + sys.argv[2] + "/" + sys.argv[1] + "*")
    dir2domainseq_dic = {}

    type_of_l = sys.argv[1]
    for dir in dirs:
        target = dir
        csv_path = "/home/user/SA-EDS/conf/input_seq_" + type_of_l +".csv"
        dirs = lr.load_10_dir(path=target)
        seq_dic = get_x_dic(csv_path, dirs)
        domain_dic = get_x_domain_dic(dirs)
        # print(seq_dic)
        # 統合した辞書を作成
        merged_dict = {}

        for key in seq_dic.keys():
            # 新しい辞書を作成してdic1とdic2の値を統合
            new_key = "/".join(key.split("/")[:-2])
            temp = key.split("_")[-2]
            merged_dict[(temp, new_key)] = {"sequence": domain_dic[key], "domain": seq_dic[key]}

        dir2domainseq_dic.update(merged_dict)

    result_path = "/home/user/SA-EDS/dataset/x_initial_" + type_of_l + ".pkl"
    print(result_path)
    with open(result_path, "wb") as tf:
        pickle.dump(dir2domainseq_dic,tf)


def make_x_dir_QD():
    if len(sys.argv) != 3:
        print("usage : python3 get_x2.py <type of l> <target>")
        print("ex)")
        print("type of l : L1")
        print("target : ../input/results/QD_1/L1_1")

    # dirs = glob.glob(sys.argv[2] + "/" + sys.argv[1] + "*")
    dirs = glob.glob(sys.argv[2] + "/*")

    dir2domainseq_dic = {}

    for dir in dirs:
        type_of_l = sys.argv[1]
        target = dir

        csv_path = "../input/input_seq_" + type_of_l +".csv"
        dirs = lr.load_10_dir(path=target)
        seq_dic = get_x_dic(csv_path, dirs)
        domain_dic = get_x_domain_dic(dirs)
        # 統合した辞書を作成
        merged_dict = {}

        for key in seq_dic.keys():
            # 新しい辞書を作成してdic1とdic2の値を統合
            new_key = "/".join(key.split("/")[:-2])
            temp = key.split("_")[-2]
            merged_dict[(temp, new_key)] = {"sequence": domain_dic[key], "domain": seq_dic[key]}

        dir2domainseq_dic.update(merged_dict)
    target=sys.argv[2].split("/")[3]
    result_path = f"../data/dic/x_{target}_" + type_of_l + ".pkl"
    with open(result_path, "wb") as tf:
        pickle.dump(dir2domainseq_dic,tf)

if __name__ == '__main__':
    make_x_dir_initial()
    # make_x_dir_QD()
