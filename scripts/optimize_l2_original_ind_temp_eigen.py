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
import random

import connectivity as cn

class L2Individual(Individual):
    strands: Iterable[tuple[Strand,float]]
    complexes: ComplexSet
    evaluation: Any
    # input/structure_seq/input_seq_L2.csv
    # [0, 1, 1, 0 ... 1] : a a*, a b, ... b* b*
    indexes: [bool] 
    structure : [str] 
    # a : nupack.core.Domain
    # b : nupack.core.Domain
    a : Strand
    b : Strand
    temp_lst : [int] = [277, 298, 308, 318, 328, 338, 348, 358]
    temp : int
    input : [int]
    # export_seq.py にて生成
    seq_dic = {'a': ['AGAAT', 'AGAAT', 'AGGGG', 'CTCTC', 'ATCGC', 'CGTGT', 'GAAGC', 'CACCG', 'TATGA', 'CATCG', 'GTGAC', 'CTACC', 'AACCG', 'CGGAA', 'TCAGA', 'CGACC', 'CAACG', 'CATCT', 'GAAGT', 'ACCGT'], 'b': ['TAACC', 'TAACC', 'ATGAA', 'AAACT', 'GTCCG', 'CCTAA', 'GTGGT', 'TTTGC', 'GTTCC', 'CTACC', 'AGATT', 'AGTTT', 'CTTTT', 'TCCGA', 'CGGGG', 'CAAGC', 'GCTTG', 'TTGGT', 'GGGGG', 'CCCTT'], 'c': [], 'd': [], 'e': [], 'f': []}
    domains : {}
    first_indexes = [[2, 1, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7], [2, 1, 11, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7]]
    eigenvalue_2 : float

    def __init__(self, specs:dict[str,Any]={}, indexes: [bool] = [], a_string = "CCCTT", b_string = "GAGTC", nb_comb = 256, **kwargs : Any) -> None:
        seq_num = random.randint(0, 19)
        a_string = self.seq_dic['a'][seq_num]
        b_string = self.seq_dic['b'][seq_num]
        self.domains = {}
        self.domains['a'] = a_string
        self.domains['b'] = b_string
        self.structure = self.init_structure()
        self.specs = specs
        self.nb_comb = nb_comb
        self.temp = 277
        # self.indexes = [0]*self.nb_comb
        self.gen_indexes(seq_num)
        self.temp = 277
        # strands from indexes
        self.a = Domain(a_string, name='a')
        self.b = Domain(b_string, name='b')
        self.strands = self.indexes2strands(indexes)
        if "max_size" not in self.specs:
            self.specs["max_size"] = 3
        self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**specs))
        self.evaluation = None
        self.eigenvalue_2 = cn.ind2eigen(type_of_l="L2", indexes=self.indexes, structure=self.structure, temperature=self.temp, domains=self.domains)
        self.input = self.indexes + [self.temp] + [self.eigenvalue_2]
        super().__init__(iterable=self.strands,**kwargs)
        #self.tube = Tube(dict(self.strands), complexes=SetSpec(**specs), name =f"tube_{self.name}")

    def gen_indexes(self, seq_num):
        self.indexes = [0]*self.nb_comb
        for index in self.first_indexes[seq_num]:
            self.indexes[index] = 1

    
    def init_structure(self):
        # input/structure_seq/input_seq_L2.csv
        f = open("../../input/structure_seq/input_seq_L2.csv", "r")
        lst = csv.reader(f)
        structure = []
        for l in lst:
            for elem in l:
                structure.append(elem)
        return structure

    def reinit(self):
        strands_set = set([s for s,_ in self.strands])
        strands_lst = list(strands_set)
        # self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**self.specs))
        self.complexes = ComplexSet(strands=strands_lst, complexes=SetSpec(**self.specs))
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
                    if domain == "a":
                        tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.a), name="x")
                    elif domain == "b":
                        tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.b), name="x")
                    elif domain == "a*":
                        tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.a))), name="x")
                    elif domain == "b*":
                        tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.b))), name="x")
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
        if not isinstance(other, L2Individual):
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
class L2GenIndividuals(GenIndividuals):
    def __next__(self):
        return L2Individual()

@registry.register
class L2Evo(Evolution):
    strands_number_domain: DomainLike
    sel_pb: float
    init_pb: float
    mut_pb: float

    def __init__(self, container: Container, budget: int,
            strands_number_domain: DomainLike = (2.,6.), concentration_domain: DomainLike = (1e-10,1e-6),
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
            remove_pb = 0.2 if ind.indexes.count(1) > strands_number_domain[0] else 0.0
            add_pb = 0.4 if len(ind.strands) < strands_number_domain[1] else 0.0
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
                select_or_initialise=select_or_initialise, vary=vary_fn, base_ind_gen=L2GenIndividuals(), **kwargs) # type: ignore

        
    def my_init_fn(self, base_ind):
        # print("called init fn")
        # generate random indexes
        indexes = [0]*base_ind.nb_comb
        nb_strands = random.choice(range(int(self.strands_number_domain[0]), int(self.strands_number_domain[1])))
        while nb_strands > 0:
            index = random.choice(range(0, 15))
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

def getModel(path="../../saved_model/L2_ave_230530"):
    model = load_model(path)
    return model

def seq_ratio_val(ind):
    base_dic = {
        'A' : 0,
        'T' : 0,
        'G' : 0,
        'C' : 0,
    }
    base_ratio = {
        'A' : 0,
        'T' : 0,
        'G' : 0,
        'C' : 0,
    }

    for strand in ind.strands:
        for base in strand[0].to_string():
            base_dic[base] += 1
    base_sum = base_dic['A'] + base_dic['T'] + base_dic['G'] + base_dic['C']
    base_ratio['A'] = float(base_dic['A'])/float(base_sum)
    base_ratio['T'] = float(base_dic['T'])/float(base_sum)
    base_ratio['G'] = float(base_dic['G'])/float(base_sum)
    base_ratio['C'] = float(base_dic['C'])/float(base_sum)

    return base_ratio



# def eval_fn(ind, nupackmodel = Model(material='dna'), fitness_scale = (1e-4)/2,max_strands = 15, complexes_scale=270):
#     """An example evaluation function. It takes an individual as input, and returns the pair ``(fitness, features)``, where ``fitness`` and ``features`` are sequences of scores."""

#     tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
#     score = 0
#     complexes = tube_results.complexes
#     for t, v in tube_results.tubes.items():
#         #print("tube:",t,v)
#         for c in t.complexes:
#             score += complexes[c].free_energy*v[c]
#     features = ((len(ind.strands)-1)/(max_strands-1), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1)))
#     ind.fitness.values = (- np.tanh(score/fitness_scale),)
#     ind.features.values = features
#     #print("DEBUG",ind.fitness, ind.features)
#     return (- np.tanh(score/fitness_scale),), features

def nupack_val(ind, nupackmodel = nupackModel(material='dna'), fitness_scale = (1e-8)/2, complexes_scale=200):
    nupackmodel = nupackModel(material='dna', kelvin=ind.temp)
    tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
    score = 0
    complexes = tube_results.complexes
    for t, v in tube_results.tubes.items():
        #print("tube:",t,v)
        for c in t.complexes:
            if v[c] == 0:
                continue
            # print(complexes[c].free_energy, v[c])
            score += complexes[c].free_energy*v[c]
    return - np.tanh(score/fitness_scale), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1))

def set_eval(ind, averageModel,scale=10.0 ):
    # fit0 = 2*math.atan(deviationModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    score = 2*math.atan(averageModel.predict([ind.input])[0]/scale)/math.pi
    energy, nb_comp = nupack_val(ind)
    fit0 = energy
    # fit1 = np.sum(strands)
    fit1 = nb_comp
    fit1 = (ind.temp - 277)/(358-277)
    # fit1 = energy
    base_ratio = seq_ratio_val(ind)
    fit1 = (ind.temp - 277.0)/(358.0-277.0)
    fit1 = base_ratio['C']
    fit1 = base_ratio['C']/0.6
    fit1 = (base_ratio['C'] - 0.05)/0.5

    features = (fit0, fit1)
    return (score,), features

def run_qdpy(dirpath="test"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
        
    grid = containers.Grid(
        shape=(16,16), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),),  
        # features_domain=((0., 1.), (1, 256))) #軸 deviation, number of strands 
        features_domain=((0., 1.), (0., 1.))) #軸 deviation, free energy
    
    algo = L2Evo(
        grid, 
        # budget=6000, 
        budget=10000, 
        batch_size=100,
        optimisation_task="maximization")
    
    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.AlgorithmLogger(algo)

    # averageModel = getModel('../../saved_model/l2_ave_230613')
    # deviationModel = getModel('../../saved_model/l2_dev_230613')
    
    # averageModel = getModel('../../saved_model/l2_ave_20230829')
    # deviationModel = getModel('../../saved_model/l2_dev_20230830')

    # averageModel = getModel('../../saved_model/l2_ave_20230903')
    with open('../saved_model/bagging_model_L2_initial.pkl', 'rb') as f:
        regr_loaded = pickle.load(f)

    averageModel = getModel('../../saved_model/l2_ave_20230903')
    deviationModel = getModel('../../saved_model/l2_dev_20230903')

    eval_fn = functools.partial(set_eval,averageModel=regr_loaded )
    best = algo.optimise(eval_fn)
    print(algo.summary())

    # Plot the results
    now = datetime.datetime.now()
    timestamp = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
    logger.final_filename = dirpath + f"/qdpy_log_L2_{timestamp}.p"
    print(logger.final_filename)
    plots.default_plots_grid(logger)
    print("All results are available in the '%s' pickle file." % logger.final_filename)

def main():
    run_qdpy()
    
if __name__ == "__main__":
    main()
