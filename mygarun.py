# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 08:07:25 2015

@author: jordan
"""

import CreateCells as cc
import myga
import inspyred
import argparse
import os, errno
from random import Random
from time import time

def main(prng=None):
    prng=None
    if prng is None:
        prng = Random()
        prng.seed(time())

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', help='Infile for building a Cell')
    parser.add_argument('dir', nargs='?',
                        default=os.path.expanduser('~/Documents/GA_Results/'),
                        help='Directory for the GA output file')
    parser.add_argument('-p', dest="pop", type=int,
                        default=2, help='population size')
    parser.add_argument('-g', dest="gen", type=int,
                        default=2, help='max generations')
    args = parser.parse_args()

    if args.infile:
        myfile = args.infile
    else:
        myfile = raw_input('Please provide the parameter file location!: ')

    try:
        os.makedirs(args.dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(args.dir):
            pass
        else:
            raise

    problem = cc.Create_Cell(myfile)

    ga = myga.Ga(problem, args.dir)

    ea = inspyred.ec.emo.NSGA2(prng)

    ea.variator = [inspyred.ec.variators.blend_crossover,
                   inspyred.ec.variators.gaussian_mutation]

    ea.terminator = inspyred.ec.terminators.generation_termination

    final_pop = ea.evolve(generator = ga.generator,
                          evaluator = ga.evaluator,
                          bounder = ga.bounder,
                          maximize = True,
                          pop_size = args.pop,
                          max_generations = args.gen)

    return ea

if __name__=='__main__':
    main()
