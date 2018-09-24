#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import libs
import os
from mpi4py import MPI

from math import sin, cos

from gaft import GAEngine
from gaft.components import BinaryIndividual, Population
from gaft.operators import RouletteWheelSelection, UniformCrossover, FlipBitMutation

# Analysis plugin base class.
from gaft.plugin_interfaces.analysis import OnTheFlyAnalysis
# Built-in best fitness analysis.
from gaft.analysis.fitness_store import FitnessStore

indv_template = BinaryIndividual(ranges=[(0, 10)], eps=0.001)
population = Population(indv_template=indv_template, size=50)
population.init()  # Initialize population with individuals.

# Use built-in operators here.
selection = RouletteWheelSelection()
crossover = UniformCrossover(pc=0.8, pe=0.5)
mutation = FlipBitMutation(pm=0.1)

engine = GAEngine(population=population, selection=selection,
                  crossover=crossover, mutation=mutation,
                  analysis=[FitnessStore])
                  
@engine.fitness_register
def fitness(indv):
    x, = indv.solution
    return x + 10*sin(5*x) + 7*cos(4*x)
    
# alternate
#@engine.fitness_register
#@engine.minimize
#def fitness(indv):
#    x, = indv.solution
#    return x + 10*sin(5*x) + 7*cos(4*x)
    
@engine.analysis_register
class ConsoleOutput(OnTheFlyAnalysis):
    master_only = True
    interval = 1
    def register_step(self, g, population, engine):
        best_indv = population.best_indv(engine.fitness)
        msg = 'Generation: {}, best fitness: {:.3f}'.format(g, engine.fmax)
        engine.logger.info(msg)
        
if '__main__' == __name__:
    comm = MPI.COMM_WORLD
    print("Hello! I'm rank %d running on node %s from %d running in total." % (comm.rank, os.uname()[1], comm.size))
    comm.Barrier() # wait for sync here
    engine.run(ng=100)