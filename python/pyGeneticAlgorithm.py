import numpy as np
import random
import time

def fitfunc(func):
    def f(self,*args,**kwargs):
        self._fitness = np.power(2,np.sum(self._population,axis=1,dtype=float))
        func(self)
    return f

class GeneticAlgorithm(object):
    def __init__(self,crossrate,mutate_rate,pop_size,dna_size,iteration,sel_range,sel_size,maxsec):
        self._cross_rate = crossrate
        self._mutate_rate = mutate_rate
        self._pop_size = pop_size
        self._dna_size = dna_size
        self._iteration = iteration
        self._sel_range = sel_range
        self._sel_size = sel_size
        self._maxsec = maxsec
        self._population = None
        self._population_next = None
        self._fitness = None
        self._chance = None
        self._best = None

    def initpop(self):
        self._population = np.zeros(shape=(self._pop_size, self._dna_size))
        for i in range(self._pop_size):
            for j in range(self._dna_size):
                self._population[i,j] = self._sel_range[random.randint(0,self._sel_size-1)]
        self._fitness = np.zeros(shape=(self._pop_size))
        self._chance = np.zeros(shape=(self._pop_size))
        self._population_next = np.zeros((self._pop_size,self._dna_size))

    def cross(self):
        for i in range(self._pop_size):
            if random.random() <= self._cross_rate:
                num1_rand = random.randint(0,self._pop_size-1)
                num2_rand = random.randint(0,self._pop_size-1)
                pos_rand = random.randint(0,self._dna_size-1)
                temp = self._population[num1_rand,pos_rand:self._dna_size]
                self._population[num1_rand,pos_rand:self._dna_size] = self._population[num2_rand,pos_rand:self._dna_size]
                self._population[num2_rand, pos_rand:self._dna_size] = temp

    def mutation(self):
        for i in range(self._pop_size):
            if random.random() <= self._mutate_rate:
                pos_rand = random.randint(0, self._dna_size-1)
                sel_rand = random.randint(0, self._sel_size-1)
                self._population[i,pos_rand] = self._sel_range[sel_rand]
    @fitfunc
    def fresh(self):
        sumfit = np.sum(self._fitness)
        self._chance[0] = self._fitness[0]/sumfit
        for i in range(1,self._pop_size,1):
            self._chance[i] = self._chance[i-1] + self._fitness[i]/sumfit

    def select(self):
        xh = np.argmax(self._fitness,axis=0)
        self._best = self._population[np.argmax(self._fitness,axis=0)]
        for i in range(self._pop_size):
            chance = random.random()
            if chance < self._chance[0]:
                self._population_next[i,:] = self._population[0,:]
            for j in range(self._pop_size-1):
                if chance >= self._chance[j] and chance < self._chance[j+1]:
                    self._population_next[i,:] = self._population[j+1,:]
        self._population = self._population_next
        self._chance = np.zeros(shape=(self._pop_size))
        self._fitness = np.zeros(shape=(self._pop_size))

    def generate(self):
        self.initpop()
        self.fresh()
        timest = time.time()
        for i in range(self._iteration):
            timeend = time.time()
            if (timeend-timest)>self._maxsec:
                break
            self.select()
            self.cross()
            self.mutation()
            self.fresh()

    @property
    def best(self):
        return self._best
