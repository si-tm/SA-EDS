import pickle
import networkx as nx
import matplotlib.pyplot as plt
import time
import numpy.linalg
import networkx as nx
import numpy as np
from scipy.sparse import csgraph
import nupack 
from nupack import *
from nupack import Model as nupackModel
import sys

"""
type of l : L1, L2, L3
indexes : [0, 1, 1, 0 ... 1]
structure : ["a a", "a a*","a b","a b*","a* a","a* a*","a* b","a* b*","b a","b a*","b b","b b*","b* a","b* a*","b* b","b* b*"]
temperature : 277

G = getGraph_deltaG(dic, seq, temp,type_of_l="L3")
e = calcEigenvalues(G)
"""

def ind2eigen(type_of_l, indexes, structure, temperature, domains):
    dic = {}
    for i, string in enumerate(structure):
        dic[string] = indexes[i]
    G = getGraph_deltaG(dic, domains, temperature,type_of_l)
    e = calcEigenvalues(G)
    return e
    


def getGraph_deltaG(dic, seq, temp, type_of_l="L1"):

    # get domain
    if type_of_l == "L1" or type_of_l == "L2":
        a = Domain(str(seq['a']), name='Domain a', material='dna')
        b = Domain(str(seq['b']), name='Domain b', material='dna')
        domain_name = {
            "a" : a,
            "b" : b,
            "a*" : ~a,
            "b*" : ~b
        }
    elif type_of_l == "L3":
        domain_name = {
            "a" : a,
            "b" : b,
            "c" : c,
            "d" : d,
            "e" : e,
            "f" : f,
            "a*" : ~a,
            "b*" : ~b,
            "c*" : ~c,
            "d*" : ~d,
            "e*" : ~e,
            "f*" : ~f
        }
    else:
        print("input type of l")
        return
    

    # generate graph
    G = nx.Graph()
    nodes = set()
    strands = {}
    num = 0

    # get strands
    for key in dic:
        if dic[key] == 1:
            # generate strand
            seq = ""
            for domain in key.split(' '):
                seq += domain_name[domain].to_string()
            strands[f"s{num}"] = Strand(seq, name=f's{num}')
            num += 1

    # calc delta G
    s_conc = {}
    for str_name in strands:
        # s_conc[str_name] = 1e-8
        s_conc[strands[str_name]] = 1e-8

    t1 = Tube(strands=s_conc, complexes=SetSpec(max_size=2), name='t1')
    model1 = nupackModel(material='dna', kelvin=float(temp))
    tube_results = tube_analysis(tubes=[t1], model=model1)

    # add edges
    for complex in tube_results.complexes:
        # print(complex.strands[0].name, complex.strands[1].name, tube_results.complexes[complex].free_energy)
        # print(complex.strands)
        if tube_results.complexes[complex].free_energy*(-1) < 0:
            continue
        if len(complex.strands) == 2:
            G.add_edge(complex.strands[0].name, complex.strands[1].name, weight=tube_results.complexes[complex].free_energy*(-1))
        elif len(complex.strands) == 1:
            G.add_edge(complex.strands[0].name, complex.strands[0].name, weight=tube_results.complexes[complex].free_energy*(-1))
        
    # plot
    # options = {
    #     "font_size": 12,
    #     "node_size": 3000,
    #     "node_color": "white",
    #     "edgecolors": "black",
    #     "linewidths": 5,
    #     "width": 5,
    # }

    # nx.draw_networkx(G, **options)

    # # Set margins for the axes so that nodes aren't clipped
    # ax = plt.gca()
    # ax.margins(0.20)
    # plt.axis("off")
    # plt.show()

    return G

def calcEigenvalues(G):
    # どのノードとも繋がっていないノードを見つけて削除します
    # isolated_nodes = [node for node, degree in dict(G.degree()).items() if degree == 0]
    # G.remove_nodes_from(isolated_nodes)

    L = nx.normalized_laplacian_matrix(G)
    
    e = numpy.linalg.eigvals(L.toarray())
    e_sort = numpy.sort(e)
    # print(len(e_sort), G)
    if len(e_sort) == 1:
        return e_sort[0]
    return e_sort[1]
   
# def main():
#     e = ind2eigen(type_of_l="L1",
#               indexes=[1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#               structure=["a a", "a a*","a b","a b*","a* a","a* a*","a* b","a* b*","b a","b a*","b b","b b*","b* a","b* a*","b* b","b* b*"],
#               temperature=277,
#               domains={'a' : 'CGGCCAGTAA', 'b' : 'GCCGGTGAAC'})
#     print(e)

def main():
    if len(sys.argv) != 3:
        print("usage : python3 connectivity L{1-3} target")
    type_of_l = sys.argv[1]
    target = sys.argv[2]
    file_path=f"home/user/SA-EDS/dataset/x_{target}_{type_of_l}.pkl"
    f = open(file_path, "rb")
    data = pickle.load(f)
    for (temp, dir) in data:
        print(temp, dir)
        dic = data[(temp, dir)]["domain"]
        seq = data[(temp, dir)]["sequence"]
        print(data[(temp, dir)])
        G = getGraph_deltaG(dic, seq, temp,type_of_l=type_of_l)
        e = calcEigenvalues(G)
        data[(temp, dir)]["eigenValue_2"] = e
        print(data[(temp, dir)])
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)
   
if __name__ == '__main__':
    main()
