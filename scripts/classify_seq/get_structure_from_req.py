import sys
# script/から実行
sys.path.append("common")
sys.path.append(".")
import common.get_target_file as gtf

class req():
    def __init__(self, target_dir):
        self.req = gtf.get_req(target_dir)
        self.seq2req_dic = {}  
        self.req2seq_dic = {}
        # self.get_structure_req()
    
    def get_structure_seq(self, target_dir):
        self.req = gtf.get_req(target_dir)
        f = open(self.req, "r")
        structures = []
        for l in f:
            if l[0] == "s":
                # "="と"@initial"の間を出力
                lst = l.split(" ")
                # structure = lst[lst.index("=") + 1: lst.index("@")]
                if "@ initial" in lst:
                    structure = lst[lst.index("=") + 1: lst.index("@ initial")]
                    print(structure)
                elif "@initial" in lst:
                    structure = lst[lst.index("=") + 1: lst.index("@initial")]
                else:
                    print("wow")
                    structure = lst[lst.index("=") + 1: lst.index("@")]
                structures.append(structure)
        return structures
    
    def get_domain_seq(self, target_dir):
        self.req = gtf.get_req(target_dir)
        f = open(self.req, "r")
        domains = {}
        for l in f:
            if l[:6] == "domain":
                # domain a = CGGCCAGTAA
                domains[l.split(" ")[1]] = l.split(" ")[3].split('\n')[0]
        return domains

def seq2structure(target_dir):
    req_data = req(target_dir)
    return req_data.get_structure_seq(target_dir)

def seq2domain(target_dir):
    req_data = req(target_dir)
    return req_data.get_domain_seq(target_dir)



def main():
    str = seq2structure("../input/results/oxdna_random_1/L1/d-15-9-2-14-1-0/L1_d-15-9-2-14-1-0_0/L1_d-15-9-2-14-1-0_0/")
    print(str)

if __name__ == '__main__':
  main()