import pickle
import numpy as np
import csv
import os
import glob
import sys

if sys.argv[1] == "L1":
    from qd_l1 import L1Individual
if sys.argv[1] == "L2":
    from qd_l2 import L2Individual
if sys.argv[1] == "L3":
    from qd_l3 import L3Individual

# final.pからスコアの高いindを取り出す

def make_req(type_of_l, filename, lst, target, result_path):
    temp_f = open("/home/user/SA-EDS/conf/requirement_" + type_of_l + ".txt", "r")
    if os.path.isdir(result_path + "r" + filename) == False:
        os.mkdir(result_path + "r" + filename) 
    # if os.path.isdir(result_path + "r" + target + "/r" + filename ) == False:
    #     os.mkdir(result_path + "r" + target + "/r" + filename ) 
    new_f = open(result_path + "r" + filename + "/req_r" + filename + ".txt", "w") # ここ変える
    tmp_theme = ""
    for temp_l in temp_f:
        if temp_l[0] == '#':
            tmp_theme = temp_l[:-1]
            if temp_l[:-1] == "# structure":
                new_f.write("# structure\n")
                new_f.write("number_of_types = " + str(len(lst)) + "\n")
                for l in lst:
                    new_f.write(l + "\n")

        if "# structure" == tmp_theme:
            pass
        else:
            new_f.write(temp_l[:-1] + "\n")
            

def req(comp):
    new_comp = []
    for i, c in enumerate(comp):
        new_comp.append("s" + str(i) + " = " + c + " @initial 1.0 M")
    return new_comp

def Ind2complexes(lst, type_of_l):
    f = open("/home/user/SA-EDS/conf/input_seq_" + type_of_l + ".csv", "r")
    r = csv.reader(f)

    comp = []

    seq_lst = []
    for l in r:
        for e in l:
            seq_lst.append(e)
    
    for i, e in enumerate(lst):
        if e == 1:
            # print(seq_lst[i])
            comp.append(seq_lst[i])
    
    return comp


# https://gitlab.com/leo.cazenille/qdpy/-/blob/master/qdpy/containers.py
def readFinal(path, type_of_l, target, result_path):
    with open(path, "rb") as f:
        # data = pickle.load(f)
        data = pickle.load(f)  
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    ind_dic = {}
    for i, ind in enumerate(grid):
        # if ind.features[1] >= 2 and ind.features[1] <= 6: # ここ変える
        lst = ind.indexes
        comp = Ind2complexes(lst, type_of_l)
        comp = req(comp)
        if not ind.name in ind_dic:
            ind_dic[ind.name] = ind
        elif ind_dic[ind.name].fitness[0] < ind.fitness[0]:
            ind_dic[ind.name] = ind

    for i, ind in enumerate(ind_dic):
        make_req(type_of_l=type_of_l, filename=f"{ind.name}", lst=comp, target=target, result_path=result_path) # ここ変える

    make_req(type_of_l=type_of_l, filename=f"{ind.name}", lst=comp, target=target, result_path=result_path) # ここ変える

def write_csv(path, type_of_l, target):
    with open(path, "rb") as f:
        # data = pickle.load(f)
        data = pickle.load(f)  

    grid = data['container']
    with open("/home/user/SA-EDS/results/qd_output_{type_of_l}_{target}.csv", "w") as f:
        f.write("fitness,scaled_energy,b/c,temp,strand_set\n")  # ヘッダーを書く
        for i, ind in enumerate(grid):
            lst = ind.indexes
            comp = Ind2complexes(lst, type_of_l)
            # comp = req(comp)
            fitness = ind.fitness[0]
            scaled_energy = ind.features[0]
            bc = ind.features[1]
            temp = ind.temp
            strand_set = ",".join([strand.name for strand, conc in ind.strands])
            

            comp = ",".join(comp)
            
            f.write(f"{fitness},{scaled_energy},{bc},{temp},{comp} \n")
 

                # make_req(type_of_l=type_of_l, filename=ind.name, lst=comp, target=target)  # 必要に応じてコメントアウトを解除


def req(comp):
    new_comp = []
    for i, c in enumerate(comp):
        new_comp.append("s" + str(i) + " = " + c + " @initial 1.0 M")
    return new_comp

def Ind2complexes(lst, type_of_l):
    f = open("/home/user/SA-EDS/conf/input_seq_" + type_of_l + ".csv", "r")
    r = csv.reader(f)

    comp = []

    seq_lst = []
    # for l in r:
    #     for e in l:
    #         seq_lst.append(e)
    
    for l in r:
        seq_lst.append(" ".join(l))

    
    for i, e in enumerate(lst):
        if e == 1:
            # print(seq_lst[i])
            comp.append(seq_lst[i])



    return comp



def readFinals(path, type_of_l):
    # script/optimize_tutorial/results/optimizationresults_20230724060622/final
    finals = glob.glob(path + "final/*")
    for final in finals:
        readFinal(final, type_of_l)


def main(target, type_of_l, result_path):
    # path="results/optimizationresults_" + target + "/final.p" # ここ変える
    write_csv(target, type_of_l, target)
    readFinal(target, type_of_l, target, result_path)
    # readFinal(path, type_of_l, target) # ここ変える

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("usage : python3 get_ind.py <target name> <type of l> <req_num>")
        print("usage : python3 get_ind.py int_initial L1 1")
    print(sys.argv)
    type_of_l = sys.argv[1]
    target = sys.argv[2]
    req_num = sys.argv[3]
    
    final_path = f"/home/user/SA-EDS/results/{target}_{type_of_l}/final.p"
    result_path = f"/home/user/SA-EDS/conf/req_{type_of_l}_{req_num}/"
    print(result_path)
    main(final_path, type_of_l, result_path)
