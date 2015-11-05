# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:38:17 2015

@author: jordan
"""
import CreateCells as cc
import pandas as pd
from numpy.random import choice

filename = '/mnt/Data/GA_Results/GA_2015_Nov_03_1.h5'
ab = cc.Create_Cell('ABfile')

def read_hdf(filename):
    #set empty dataframes for fitness, results, gbars
    fitness = pd.DataFrame()
    results = pd.DataFrame()
    gbars = pd.DataFrame()
    #get the amount of nodes (generations)
    with pd.HDFStore(filename) as store:
        num_gen = store._handle.root._v_nchildren
    #for each node....
    for i in range(num_gen):
        #get the name of the current node and subfolder
        name = '/gen%d/fitness' %i
        name2 = '/gen%d/feature' %i
        name3 = '/gen%d/gbars' %i
        #read each subfolder into a variable
        hdf = pd.read_hdf(filename, name)
        hdf2 = pd.read_hdf(filename, name2)
        hdf3 = pd.read_hdf(filename, name3)
        #rename the index since we are appending all generations together
        hdf = hdf.set_index([['gen%d_%d' %(i,x)for x in range(len(hdf))]])
        hdf2 = hdf2.set_index([['gen%d_%d' %(i,x)for x in range(len(hdf2))]])
        hdf3 = hdf3.set_index([['gen%d_%d' %(i,x)for x in range(len(hdf3))]])
        #append to the data frame
        fitness = fitness.append(hdf)
        results = results.append(hdf2)
        gbars = gbars.append(hdf3)
    return fitness, results, gbars

f, r, g = read_hdf(filename)

def good_models(fitness, threshold):
    f2 = fitness.copy()
    f2['tot'] = f2.sum(axis=1)
    good = list(f2[(f2.tot >= threshold)].index)
    print "There are {} good models with a threshold of {}".format(
                                                           len(good), threshold)
    return good

def plot_good(cell, good_indexes, sample=None, title=True):

    if sample:
        index = sorted(choice(good_indexes, size=sample, replace=False))

    elif len(good_indexes) > 15:
        x = raw_input("There are %d 'good' indexes, run anyway? (y/n) "
                                                        % len(good_indexes))
        if x.lower() == 'y':
            index = good_indexes
        else:
            print "Please rerun the function with a sample = # argument!"
            return
    else:
        index = good_indexes

    for good in index:
        gbars = g.loc[good].to_dict()
        for k, v in gbars.items():
            p = k.partition('_')
            cell._optimize[p[0]][p[2]] = v
        cell.ga_change_gbar(cell._optimize)
        out = cell.run()
        if title:
            cell.plot('vs', 'va', title=good)
        else:
            cell.plot('vs', 'va')