#!/usr/bin/env python
# coding: utf-8

from qdpy import algorithms, containers, plots
import math
import numpy as np
from keras.models import load_model

import tensorflow as tf
print("TensorFlow version:", tf.__version__)
#tf.get_logger().setLevel('ERROR')

import nupack
from nupack import *
from nupack import Model as nupackModel
from tensorflow.keras import Model
from scipy.stats import linregress
import numpy as np

from qdpy.phenotype import Individual, GenIndividuals
from typing import Iterable, Any
import numpy as np
from numpy import random

from qdpy import algorithms, containers, plots, tools
from qdpy.containers import Container
from qdpy.algorithms.evolution import Evolution
from qdpy.base import DomainLike, registry
# from qdpy import registry
import math
from functools import partial
import functools
import csv
import datetime

import pickle
import connectivity as cn

"""
use b/c as y_axis
""" 

type_of_l = "L1"
target = "int_initial"

class L1Individual(Individual):
    strands: Iterable[tuple[Strand,float]]
    complexes: ComplexSet
    evaluation: Any
    # input/structure_seq/input_seq_L1.csv
    # [0, 1, 1, 0 ... 1] : a a*, a b, ... b* b*
    indexes: [bool] 
    structure : [str] 
    domains : {}
    nb_comb : int
    temp_lst : [int] = [277, 298, 308, 318, 328, 338, 348, 358]
    # temp_lst : [int] = [277, 308, 328, 358]
    temp : int
    input : [int]
    # export_seq.py にて生成
    seq_dic = {'a': ['CGGCCAGTAA', 'CGGCCAGTAA', 'CGTGCACGTT', 'AACACTGTTC', 'TACATCAATG', 'AGAATCGATT', 'GGCCCGCGGG', 'CAGGTGTTGA', 'GCACCGTGCT', 'CTTGCGCCCG', 'CTTTTCACCC', 'CCAGGCTGGG', 'GGCTTCCGAC', 'ATCGGTCGAT', 'CAAATATTTG', 'CCTCTCGACA', 'CAGCCTGCTG', 'CGCTCAGCGT', 'CACGTCGACG', 'GGGAGCTCCT'], 'b': ['GCCGGTGAAC', 'GCCGGTGAAC', 'AGGCTCCAGC', 'GACTAACGCA', 'GGCCTAGGCT', 'TCTACTGTGC', 'TGGTCGTAAA', 'TGGATTATCC', 'CGTACTGGCG', 'TTAGCGACCT', 'AGCTTAAGCG', 'ATCCTTTTCC', 'AAGCGCTTCG', 'CAGGCACCTA', 'CGTAAACGGG', 'AACCTGGGTT', 'TCGCTCTTAT', 'GGATACCGGG', 'TATTGGTGTA', 'GTGCGTAAAT'], 'c': [], 'd': [], 'e': [], 'f': []}
    domains : {}
    first_indexes = [[0, 1, 3, 9, 11, 4, 14], [0, 1, 3, 9, 11, 4, 14], [0, 1, 11, 4, 6, 5, 12], [0, 1, 3, 9, 11, 4, 5], [1, 8, 10, 11, 6, 7, 12], [0, 1, 3, 9, 5], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7], [2, 1, 11, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7]]
    eigenvalue_2 : float

    def __init__(self, specs:dict[str,Any]={}, indexes: [bool] = [], a_string = "GTTACTTGGA", b_string = "GGTTTTTTGC", nb_comb = 16, **kwargs : Any) -> None:
        seq_num = random.randint(0, 19) 
        self.structure = self.init_structure()
        self.specs = specs
        self.nb_comb = nb_comb
        # self.indexes = [0]*self.nb_comb
        self.gen_indexes(seq_num)
        # self.temp = 277
        self.temp = random.choice(self.temp_lst)
        self.input = self.indexes + [self.temp]
        # strands from indexes
        self.domains = self.init_domains(seq_num)
        self.strands = self.indexes2strands(indexes)
        if "max_size" not in self.specs:
            self.specs["max_size"] = 3
        self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**specs))
        self.evaluation = None
        # print("L1", self.indexes, self.structure, self.temp, self.domains)
        self.eigenvalue_2 = cn.ind2eigen(type_of_l="L1", indexes=self.indexes, structure=self.structure, temperature=self.temp, domains=self.domains)
        self.input = self.indexes + [self.temp] + [self.eigenvalue_2]
        super().__init__(iterable=self.strands,**kwargs)
        #self.tube = Tube(dict(self.strands), complexes=SetSpec(**specs), name =f"tube_{self.name}")

    def gen_indexes(self, seq_num):
        self.indexes = [0]*self.nb_comb
        for index in self.first_indexes[seq_num]:
            self.indexes[index] = 1

    def init_domains(self, seq_num):
        domain_name_lst = ['a', 'b']
        # domain_seq_lst = ["GTTCCAGCACCTTCACT", "TCAGCCGGTGGACTGAG", "GGCAACGTCCTGTTACT", "GAACCCGGGAAAGAC", "TTGAAGAGTAGAGCACA", "CCTAGAGAGGCGCACAT"]
        domain_seq_lst = []
        # seq_num = random.randint(0, 19)
        for domain_name in domain_name_lst:
            domain_seq_lst.append(self.seq_dic[domain_name][seq_num])
        domains = {}
        for i, domain_name in enumerate(domain_name_lst):
            domains[domain_name] = Domain(domain_seq_lst[i], name=domain_name)
            domains[domain_name + "*"] = ~Domain(domain_seq_lst[i], name=domain_name)
        return domains
    
    def init_structure(self):
        # input/structure_seq/input_seq_L2.csv
        f = open("/home/user/SA-EDS/conf/input_seq_L1.csv", "r")
        lst = csv.reader(f)
        structure = []
        for l in lst:
            tmp = ""
            for elem in l:
                tmp += elem
                tmp += " "
            structure.append(tmp[:-1])
        # print(structure)
        return structure

    def reinit(self):
        strands_set = set([s for s,_ in self.strands])
        strands_lst = list(strands_set)
        # self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**self.specs))
        self.complexes = ComplexSet(strands=strands_lst, complexes=SetSpec(**self.specs))
        # self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**self.specs))
        self.evaluation = None
        self.tube = Tube(dict(self.strands), complexes=SetSpec(**self.specs), name =f"tube")
        self.eigenvalue_2 = cn.ind2eigen(type_of_l="L1", indexes=self.indexes, structure=self.structure, temperature=self.temp, domains=self.domains)
        self.input = self.indexes + [self.temp] + [self.eigenvalue_2]
    
    def indexes2strands(self, indexes, concentration=1e-10):
        tmp_index = 0
        strands = []
        for num, index in enumerate(indexes):
            if index == 1:
                tmp_domain = Domain("", name="s"+str(tmp_index))
                tmp_strand = Strand(str(tmp_domain), name="Strand tmp")
                for domain in self.structure[num].split(" "):
                    tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.domains[domain]), name=domain)
                    # if domain == "a":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.a), name="x")
                    # elif domain == "b":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.b), name="x")
                    # elif domain == "a*":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.a))), name="x")
                    # elif domain == "b*":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.b))), name="x")
                    tmp_domain = Domain(str(tmp_strand), name="s"+str(tmp_index))                
                self.strands.append((tmp_strand, concentration))
                tmp_index += 1
        return strands
    
    def domainlist_to_string(self, a):
        if isinstance(a,Domain):
            return a.to_string()
        if isinstance(a,DomainList):
            return ''.join([b.to_string() for b in a])

    def __eq__(self, other): 
        if not isinstance(other, L1Individual):
            # don't attempt to compare against unrelated types
            return False
        return self.input == other.input

    def eliminateStrand(self):
        nb_strands = [i for i, a in enumerate(self.indexes) if a == 1]
        self.indexes[random.choice(nb_strands)] = 0
        self.reinit()

    def addStrand(self):
        nb_strands = [i for i, a in enumerate(self.indexes) if a == 0]
        self.indexes[random.choice(nb_strands)] = 1
        self.reinit()

    def mutateStrand(self):
        nb_1_strands = [i for i, a in enumerate(self.indexes) if a == 1]
        self.indexes[random.choice(nb_1_strands)] = 0
        nb_0_strands = [i for i, a in enumerate(self.indexes) if a == 0]
        self.indexes[random.choice(nb_0_strands)] = 1
        self.reinit()

    def changeTemp(self):
        self.temp = random.choice(self.temp_lst)
        self.reinit()

@registry.register
class L1GenIndividuals(GenIndividuals):
    def __next__(self):
        return L1Individual()

@registry.register
class L1Evo(Evolution):
    strands_number_domain: DomainLike
    sel_pb: float
    init_pb: float
    mut_pb: float

    def __init__(self, container: Container, budget: int,
            strands_number_domain: DomainLike = (3.,6.), concentration_domain: DomainLike = (1e-10,1e-6),
            sel_pb: float = 0.5, init_pb: float = 0.5, mut_pb: float = 0.2, dimension: int = 2, nupack_specs = {"max_size": 3}, **kwargs):
        self.strands_number_domain = strands_number_domain
        self.concentration_domain = concentration_domain
        self.sel_pb = sel_pb
        self.init_pb = init_pb
        self.mut_pb = mut_pb
        self.nupack_specs = nupack_specs
    
        select_or_initialise = partial(
                tools.sel_or_init,
                sel_fn = tools.sel_random,
                sel_pb = 0.9, #変更
                init_fn = self.my_init_fn,
                init_pb = 0.1) #変更
        def vary_fn(ind):
            # print("nb_strands", np.sum(ind.indexes))
            remove_pb = 0.2 if np.sum(ind.indexes) > strands_number_domain[0] else 0.0
            add_pb = 0.4 if np.sum(ind.indexes) < strands_number_domain[1] else 0.0
            temp_pb = 0.3
            rand_pb = 1 - remove_pb - add_pb - temp_pb             
            mut_type = np.random.choice(4,p=[remove_pb, add_pb, temp_pb, rand_pb])
            if mut_type == 0:
                ind.eliminateStrand()
            elif mut_type == 1:
                ind.addStrand()
            elif mut_type == 2:                 
                ind.changeTemp()
            else:
                # 1の位置を一つ入れ替える
                ind.mutateStrand()
            return ind
        
        
        super().__init__(container, budget, dimension=dimension, # type: ignore
                select_or_initialise=select_or_initialise, vary=vary_fn, base_ind_gen=L1GenIndividuals(), **kwargs) # type: ignore

        
    def my_init_fn(self, base_ind):
        # print("called init fn")
        # generate random indexes
        indexes = [0]*base_ind.nb_comb
        nb_strands = random.choice(range(int(self.strands_number_domain[0]), int(self.strands_number_domain[1])))
        while nb_strands > 0:
            index = random.choice(range(0, base_ind.nb_comb - 1))
            if indexes[index] == 0:
                indexes[index] = 1
                nb_strands -= 1
        # generate strands from indexes
        base_ind.indexes = indexes
        base_ind.indexes2strands(indexes)
        base_ind.specs = self.nupack_specs
        base_ind.reinit()
        return base_ind

    def _internal_ask(self, base_ind):
        return super()._internal_ask(base_ind)    

def nupack_val(ind, nupackmodel = nupackModel(material='dna'), fitness_scale = (1e-8), complexes_scale=200):
    nupackmodel = nupackModel(material='dna', kelvin=ind.temp)
    tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
    score = 0
    complexes = tube_results.complexes
    for t, v in tube_results.tubes.items():
        #print("tube:",t,v)
        for c in t.complexes:
            if v[c] == 0:
                continue
            score += complexes[c].free_energy*v[c]
    return - np.tanh(score/fitness_scale), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1))

def import_sigmoid_model():
    with open(f'/home/user/SA-EDS/saved_model/bagging_model_{type_of_l}_{target}_sigmoid.pkl', 'rb') as f:
        regr_loaded = pickle.load(f)
    return regr_loaded

def predict_sigmoid(input, min_val=-3, max_val=5):
    regr = import_sigmoid_model()
    input_reshaped = np.array(input).reshape(1, -1)
    sigmoid = regr.predict(input_reshaped)[0]
    return (sigmoid - min_val)/(max_val - min_val)

def set_eval(ind, averageModel, scale=10.0 ):
    strands = ind.indexes
    score = 2*math.atan(averageModel.predict([ind.input])[0]/scale)/math.pi
    energy, nb_comp = nupack_val(ind)
    fit0 = ind.temp_lst.index(ind.temp)
    fit1 = energy
    fit2 = predict_sigmoid(ind.input)
    # print("fit1 : ", fit1)
    features = (fit0, fit1, fit2)
    # print(features, score)
    return (score,), features

def run_qdpy(dirpath="/home/user/SA-EDS/results/int_initial_L1"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
        
    grid = containers.Grid(
        # shape=(16,16,8), 
        shape=(9, 16,16), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),),  
        # features_domain=((0., 1.), (1, 1728))) #軸 deviation, number of strands 
        features_domain=( (0., 7.), (0., 1.), (0., 1.))) #軸 free energy, sigmoid, temp
        # features_domain=( (277., 358.), (0., 1.), (0., 1.))) #軸 free energy, sigmoid, temp
    
    algo = L1Evo(
        grid, 
        budget=100000, 
        # budget=1000,
        # budget=3000, 
        batch_size=100,
        optimisation_task="maximization")
    
    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.AlgorithmLogger(algo)

    with open('/home/user/SA-EDS/saved_model/bagging_model_L1_int_initial.pkl', 'rb') as f:
        
    # with open(f'/home/user/SA-EDS/saved_model/bagging_model_L1_initial_sigmo.pkl', 'rb') as f:
        regr_loaded = pickle.load(f)

    eval_fn = functools.partial(set_eval,averageModel=regr_loaded)
    logger.log_base_path = dirpath 
    logger.final_filename = "final.p"
    best = algo.optimise(eval_fn)
    print(algo.summary())

    # Plot the results
    
    # gridを入れる
    plots.default_plots_grid(logger) 
    # plotGridSubplots
    # xlabel, ylabel
    print("All results are available in the '%s' pickle file." % logger.final_filename)

def main():
    run_qdpy()
    
if __name__ == "__main__":
    main()
