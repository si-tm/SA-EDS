#!/usr/bin/env python
# coding: utf-8
print("")
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
from typing import Iterable, Any, MutableSequence
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

import pickle
import random

import connectivity as cn

class L1Individual(Individual):
    strands: Iterable[tuple[Strand,float]]
    complexes: ComplexSet
    evaluation: Any
    eigenvalue_2 : float
    # input/structure_seq/input_seq_L1.csv
    # [0, 1, 1, 0 ... 1] + [temperature] + [eigen]: a a*, a b, ... b* b*
    indexes: [bool] 
    structure : [str] = ["a a", "a a*","a b","a b*","a* a","a* a*","a* b","a* b*","b a","b a*","b b","b b*","b* a","b* a*","b* b","b* b*"]
    # a : nupack.core.Domain
    # b : nupack.core.Domain
    a : Strand
    b : Strand
    temp_lst : [int] = [277, 298, 308, 318, 328, 338, 348, 358]
    temp : int
    input : [int]
    # export_seq.py にて生成
    seq_dic = {'a': ['CGGCCAGTAA', 'CGGCCAGTAA', 'CGTGCACGTT', 'AACACTGTTC', 'TACATCAATG', 'AGAATCGATT', 'GGCCCGCGGG', 'CAGGTGTTGA', 'GCACCGTGCT', 'CTTGCGCCCG', 'CTTTTCACCC', 'CCAGGCTGGG', 'GGCTTCCGAC', 'ATCGGTCGAT', 'CAAATATTTG', 'CCTCTCGACA', 'CAGCCTGCTG', 'CGCTCAGCGT', 'CACGTCGACG', 'GGGAGCTCCT'], 'b': ['GCCGGTGAAC', 'GCCGGTGAAC', 'AGGCTCCAGC', 'GACTAACGCA', 'GGCCTAGGCT', 'TCTACTGTGC', 'TGGTCGTAAA', 'TGGATTATCC', 'CGTACTGGCG', 'TTAGCGACCT', 'AGCTTAAGCG', 'ATCCTTTTCC', 'AAGCGCTTCG', 'CAGGCACCTA', 'CGTAAACGGG', 'AACCTGGGTT', 'TCGCTCTTAT', 'GGATACCGGG', 'TATTGGTGTA', 'GTGCGTAAAT'], 'c': [], 'd': [], 'e': [], 'f': []}
    domains : {}
    first_indexes = [[0, 1, 3, 9, 11, 4, 14], [0, 1, 3, 9, 11, 4, 14], [0, 1, 11, 4, 6, 5, 12], [0, 1, 3, 9, 11, 4, 5], [1, 8, 10, 11, 6, 7, 12], [0, 1, 3, 9, 5], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7], [2, 1, 11, 4, 6, 7], [2, 1, 4, 6, 7], [2, 1, 11, 6, 7, 14], [8, 11, 4, 12, 13], [8, 11, 4, 12, 14, 13], [8, 12, 13], [8, 12, 14, 13], [1, 8, 9, 11, 14, 13], [2, 6, 7], [2, 1, 6, 7]]

    def __init__(self, specs:dict[str,Any]={}, indexes: [bool] = [], a_string = "GTTACTTGGA", b_string = "GGTTTTTTGC", **kwargs : Any) -> None:
        seq_num = random.randint(0, 19)
        a_string = self.seq_dic['a'][seq_num]
        b_string = self.seq_dic['b'][seq_num]
        self.domains = {}
        self.domains['a'] = a_string
        self.domains['b'] = b_string
        self.specs = specs
        # self.indexes = [0]*16
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
        # print("L1", self.indexes, self.structure, self.temp, self.domains)
        self.eigenvalue_2 = cn.ind2eigen(type_of_l="L1", indexes=self.indexes, structure=self.structure, temperature=self.temp, domains=self.domains)
        self.input = self.indexes + [self.temp] + [self.eigenvalue_2]
        super().__init__(iterable=self.strands,**kwargs)
        #self.tube = Tube(dict(self.strands), complexes=SetSpec(**specs), name =f"tube_{self.name}")

    def gen_indexes(self, seq_num):
        self.indexes = [0]*16
        for index in self.first_indexes[seq_num]:
            self.indexes[index] = 1

    def reinit(self):
        # Check for and remove duplicate strands
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
        added_domains = set()  # 重複チェックのためのセット
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
                if tmp_domain.to_string() not in added_domains:  # 重複していない場合のみ追加                
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
        # print(self.indexes, self.strands)
        nb_strands = [i for i, a in enumerate(self.indexes) if a == 1]
        self.indexes[random.choice(nb_strands)] = 0
        # print(self.indexes)
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
                sel_pb = 1,
                init_fn = self.my_init_fn,
                init_pb = 0)
        def vary_fn(ind):
            # print(strands_number_domain[0], len(ind.strands))
            remove_pb = 0.3 if np.sum(ind.indexes) > strands_number_domain[0] else 0.0
            add_pb = 0.4 if np.sum(ind.indexes) < strands_number_domain[1] else 0.0
            temp_pb = 0.2
            rand_pb = 1 - remove_pb - add_pb - temp_pb
            mut_type = np.random.choice(4,p=[remove_pb, add_pb, temp_pb, rand_pb])
            if mut_type == 0:
                ind.eliminateStrand()
            elif mut_type == 1:
                ind.addStrand()
            elif mut_type == 2:
                ind.changeTemp()
            return ind
        
        
        super().__init__(container, budget, dimension=dimension, # type: ignore
                select_or_initialise=select_or_initialise, vary=vary_fn, base_ind_gen=L1GenIndividuals(), **kwargs) # type: ignore

        
    def my_init_fn(self, base_ind):
        # print("called init fn")
        # generate random indexes
        indexes = [0]*16
        temp = 277
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

def getModel(path="../../saved_model/l1_ave_230530"):
    model = load_model(path)
    return model


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


def nupack_val(ind, nupackmodel = nupackModel(material='dna'), fitness_scale = (1e-8)/2, complexes_scale=100):
    nupackmodel = nupackModel(material='dna', kelvin=ind.temp)
    tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
    score = 0
    complexes = tube_results.complexes
    for t, v in tube_results.tubes.items():
        #print("tube:",t,v)
        for c in t.complexes:
            score += complexes[c].free_energy*v[c]
    return - np.tanh(score/fitness_scale), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1))

# def set_eval(ind, averageModel, deviationModel, scale=10.0, ):
def set_eval(ind, averageModel, scale=10.0, ):
    base_ratio = seq_ratio_val(ind)
    # print(ind.input)
    # [1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    energy, nb_comp = nupack_val(ind)
    # score = 2*math.atan(averageModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    # fit0 = 2*math.atan(deviationModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    # fit0 = 2*math.atan(averageModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    score = 2*math.atan(averageModel.predict([ind.input])[0]/scale)/math.pi
    fit0 = energy
    fit1 = nb_comp
    fit1 = (ind.temp - 277.0)/(358.0-277.0)
    fit1 = (base_ratio['C'] - 0.1)/0.4
    features = (fit0, float(fit1))
    # print(ind.input,score, features)

    # print(features, energy)
    return (score,), features

def my_callback(algo,ind, added, xattr):
    container = algo.container
    index = container._index_grid_features(ind.features)
    # print("TELL", ind.fitness, ind.features, index, added, xattr)
    if added:
        value = container.quality_array[index]
        if np.isnan(value):
            print("CHECK CONTAINER", container._size)
            index2 = super(containers.Grid,container)._add_internal(ind, True, False)
            print("FORCE ADD", container._size,index2)
            can_be_added, index3 = container._add_to_collection(container.items,ind)
            print("MORE DIRECT ADD", can_be_added,index3)
    

def run_qdpy(dirpath="/home/user/SA-EDS/results"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
        
    grid = containers.Grid(
        shape=(16,16), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),),  
        features_domain=((0, 1), (0, 1))) #軸 ave, number of strands 
    
    algo = L1Evo(
        grid, 
        # budget=10000, 
        budget=10000, 
        batch_size=100,
        optimisation_task="maximization")
    
    algo.add_callback("tell", my_callback)
    # print(isinstance(grid.items,MutableSequence))
    
    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.AlgorithmLogger(algo)

    # f = open('../saved_model/randomforest_model.pkl', 'rb')
    # regr_loaded = pickle.load(f)
    # f.close()

    with open('/home/user/SA-EDS/saved_model/bagging_model_L1_initial.pkl', 'rb') as f:
        regr_loaded = pickle.load(f)

    eval_fn = functools.partial(set_eval,averageModel=regr_loaded)
    logger.log_base_path = dirpath 
    logger.final_filename = "final.p"
    best = algo.optimise(eval_fn)
    print(algo.summary())

    # Plot the results
    
    plots.default_plots_grid(logger) # plotされたものがどこかわからない final.pがどうしても実行時のディレクトリになってしまう
    # print("All results are available in the '%s' pickle file." % logger.final_filename)
    print("All results are available in the '%s' pickle file." % logger.final_filename)

def main():
    run_qdpy()
    
if __name__ == "__main__":
    main()
