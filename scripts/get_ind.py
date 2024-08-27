import pickle
import numpy as np
import csv
import os
import glob
import sys
from optimize_l3_original_ind_temp_eigen import L3Individual

# if sys.argv[2] == "L1":
#     from optimize_l1_original_ind_temp_eigen import L1Individual
# if sys.argv[2] == "L2":
#     from optimize_l2_original_ind_temp_eigen import L2Individual
# if sys.argv[2] == "L3":
#     from optimize_l3_original_ind_temp_eigen import L3Individual

# final.pからスコアの高いindを取り出す

def make_req(type_of_l, filename, lst, target, result_path):
    temp_f = open("/home/user/SA-EDS/conf/requirement_" + type_of_l + ".txt", "r")
    print(result_path + "r" + filename)
    if os.path.isdir(result_path + "r" + filename) == False:
        os.mkdir(result_path + "r" + filename) 
    # if os.path.isdir(result_path + "r" + target + "/r" + filename ) == False:
    #     os.mkdir(result_path + "r" + target + "/r" + filename ) 
    print(result_path + "r" + filename + "/req_r" + filename + ".txt", "w") 
    new_f = open(result_path + "r" + filename + "/req_r" + filename + ".txt", "w") # ここ変える
    tmp_theme = ""
    for temp_l in temp_f:
        if temp_l[0] == '#':
            tmp_theme = temp_l[:-1]
            if temp_l[:-1] == "# structure":
                new_f.write("# structure\n")
                new_f.write("number_of_types = " + str(len(lst)) + "\n")
                for l in lst:
                    print(type(l), l)
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
        data = pickle.load(f)
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    for ind in grid:
        if 1:
        # if ind.features[1] >= 2 and ind.features[1] <= 6: # ここ変える
            lst = ind.indexes
            comp = Ind2complexes(lst, type_of_l)
            comp = req(comp)
            print(f"fitness : {ind.fitness[0]} scaled energy : {ind.features[0]}, b/c : {ind.features[1]}, temp : {ind.temp}") # feature0 : よこ feature1 : たて
            print(f"strand set : {[str(strand) for strand, conc in ind.strands]}")
            print()
            make_req(type_of_l=type_of_l, filename=ind.name, lst=comp, target=target, result_path=result_path) # ここ変える

            # print(ind.fitness)
            # print(ind.features)

def write_csv(path, type_of_l, target):
    with open(path, "rb") as f:
        data = pickle.load(f)

    grid = data['container']
    with open("/home/user/SA-EDS/results/qd_output_L3.csv", "w") as f:
        f.write("fitness,scaled_energy,b/c,temp,strand_set\n")  # ヘッダーを書く
        for ind in grid:
            lst = ind.indexes
            comp = Ind2complexes(lst, type_of_l)
            # comp = req(comp)
            fitness = ind.fitness[0]
            scaled_energy = ind.features[0]
            bc = ind.features[1]
            temp = ind.temp
            print([strand.name for strand, conc in ind.strands])
            strand_set = ",".join([strand.name for strand, conc in ind.strands])
            
            print(f"{ind.features[0]} {ind.features[1]} {ind.features[2]}")
            print(f"fitness : {fitness} scaled energy : {scaled_energy}, b/c : {bc}, temp : {temp}")
            print(f"strand set : {strand_set}")
            print()

            comp = ",".join(comp)
            
            f.write(f"{fitness},{scaled_energy},{bc},{temp},{comp} \n")

                # make_req(type_of_l=type_of_l, filename=ind.name, lst=comp, target=target)  # 必要に応じてコメントアウトを解除


def req(comp):
    print(comp)
    new_comp = []
    for i, c in enumerate(comp):
        new_comp.append("s" + str(i) + " = " + c + " @initial 1.0 M")
    return new_comp

def Ind2complexes(lst, type_of_l):
    f = open("/home/user/SA-EDS/conf/input_seq_" + type_of_l + ".csv", "r")
    r = csv.reader(f)

    print(sum(lst))

    comp = []

    seq_lst = []
    # for l in r:
    #     for e in l:
    #         seq_lst.append(e)
    
    for l in r:
        seq_lst.append(" ".join(l))

    
    for i, e in enumerate(lst):
        if e == 1:
            print(seq_lst[i])
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
    # if len(sys.argv) != 3:
    #     print("usage : python3 get_ind.py <target name> <type of l>")
    # target = sys.argv[1]
    # type_of_l = sys.argv[2]
    target = "/home/user/SA-EDS/results/int_L3/final.p"
    # target = "/home/user/SA-EDS/results/final.p"
    type_of_l = "L3"
    result_path = "/home/user/SA-EDS/conf/req_L3_1/"
    main(target, type_of_l, result_path)
