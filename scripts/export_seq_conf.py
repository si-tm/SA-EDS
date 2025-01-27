import sys
import glob
import common
import common.get_target_file as gtf
import classify_seq.get_structure_from_req as gr
import csv

# reqから以下を作成
# seq_dic = {'a': ['AGAAT', 'AGAAT', 'AGGGG', 'CTCTC', 'ATCGC', 'CGTGT', 'GAAGC', 'CACCG', 'TATGA', 'CATCG', 'GTGAC', 'CTACC', 'AACCG', 'CGGAA', 'TCAGA', 'CGACC', 'CAACG', 'CATCT', 'GAAGT', 'ACCGT'], 'b': ['TAACC', 'TAACC', 'ATGAA', 'AAACT', 'GTCCG', 'CCTAA', 'GTGGT', 'TTTGC', 'GTTCC', 'CTACC', 'AGATT', 'AGTTT', 'CTTTT', 'TCCGA', 'CGGGG', 'CAAGC', 'GCTTG', 'TTGGT', 'GGGGG', 'CCCTT'], 'c': [], 'd': [], 'e': [], 'f': []}
# first_indexes = [[2, 1, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7], [2, 1, 11, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7]]

# reqを収集
def get_req_lst(type_of_l, target):
    req_lst = []
    dir_lst = glob.glob(f"/home/user/SA-EDS/{target}/{type_of_l}_*")
    for dir in dir_lst:
        print(dir)
        d_lst = glob.glob(f"{dir}/*")
        for d in d_lst:
            req = gtf.get_req(d)
            req_lst.append(req)
            break
    return req_lst

def get_target_lst(type_of_l, target):
    target_lst = []
    req_lst = glob.glob(f"/home/user/SA-EDS/{target}/r*/*")
    return req_lst

# それぞれのreqから、seqとindexを取得
def get_seq(target_lst, type_of_l):
    seq_dic = {
        'a': [],
        'b': [],
        'c': [],
        'd': [],
        'e': [],
        'f': []
    }


    for target in target_lst:
        f = open(target, "r")
        domains = {}
        for l in f:
            if l[:6] == "domain":
                # domain a = CGGCCAGTAA
                domains[l.split(" ")[1]] = l.split(" ")[3].split('\n')[0]
        tmp_dic = domains
        for key in tmp_dic:
            seq_dic[key].append(tmp_dic[key])
    print(seq_dic)


def get_indexes(target_lst,type_of_l):

    domain_dic = {}
    f = open(f"/home/user/SA-EDS/conf/input_seq_{type_of_l}.csv")
    csv_f = csv.reader(f)
    for i, seq in enumerate(csv_f):
        seq_str = "".join(seq)
        domain_dic[seq_str] = i

    first_indexes = []

    for target in target_lst:
        # [['a', 'a', 'a', 'a*'], ['a', 'a', 'a', 'b'], ['a', 'a', 'a*', 'a'], ['a', 'a', 'b', 'b*'], ['b', 'b*', 'b*', 'a'], ['b*', 'b', 'a', 'a']]
        f = open(target, "r")
        structures = []
        for l in f:
            if l[0] == "s":
                # "="と"@initial"の間を出力
                lst = l.split(" ")
                # structure = lst[lst.index("=") + 1: lst.index("@")]
                if "@initial" in lst:
                    structure = lst[lst.index("=") + 1: lst.index("@initial")]
                    print(structure)
                elif "@ initial" in lst:
                    structure = lst[lst.index("=") + 1: lst.index("@ initial")]
                    print(structure)
                elif "@" in lst:
                    print("wow")
                    structure = lst[lst.index("=") + 1: lst.index("@")]
                else:
                    structure = lst[lst.index("=") + 1: lst.index("@initial")]
                structures.append(structure)
        tmp_d_lst = []
        for d in structures:
            d_str = "".join(d)
            tmp_d_lst.append(domain_dic[d_str])
        first_indexes.append(tmp_d_lst)  

    print(first_indexes)      

if __name__ == '__main__':
    arg_lst = [item for item in sys.argv if item]
    if len( arg_lst ) != 3:
        print("usage: python3 script/export_seq.py <type_of_l> <req_dir>")
        # python3 script/export_seq.py L3 conf req_L3_1
    
    type_of_l = arg_lst[1]
    target = f"conf/{arg_lst[2]}"

    print(type_of_l)
    print(target)

    target_lst = get_target_lst(type_of_l, target)
    get_seq(target_lst, type_of_l)
    get_indexes(target_lst, type_of_l)

    