# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:13:35 2014

@author: Serdar Ozsezen (c)

License: BSD-3-Clause

Simple demonstration of Monte Carlo and 
genetic algorithm on energy of a system 
with four spins. 

"""
from numpy.random import rand, randint
from math import exp
from copy import deepcopy as equals

class System:
    """

    Defines the system with it's properties (spins and energy). 
    
    """
    def __init__(self,spins=None):
        if spins is None:
            self.spins = []
        else:
            self.spins = spins
            
    def __str__(self):
        d = {1: "Up", 0: "Down"}
        return "("+ d[self.spins[0]] + ", " + d[self.spins[1]] + ", " + \
        d[self.spins[2]] + ", " + d[self.spins[3]] + ")" + \
        " Energy: " + str(self.energy()) + "R"

    def __getitem__(self,key):
        return self.spins[key]

    def __setitem__(self,idx,val):
        self.spins[idx] = val

    def energy(self):
        """
        Upward = 1 (0.5R)
        Downward = 0 (0.5R)
        Up, up = 1, 1 (0.1R)
        Down, down = 0, 0 (0.1R)
        Up, down or down, up = 1, 0 or 0, 1 (-0.3R)
        """
        energy = 2.0 # Sum without secondary interactions
        for i in range(len(self.spins)):
            if i+1 < len(self.spins):
                if self.spins[i] == self.spins[i+1]:
                    energy += 0.1
                elif self.spins[i] != self.spins[i+1]:
                    energy -= 0.3
        return energy

def randomize():
    """
    Make a random list of ones and zeros for creating a random system.
    """
    random_list = []
    for i in range(4):
        if rand() > 0.5: # If random number is higher than 0.5 take 1 (Up)
            random_list.append(1)
        else: # If random number is less than or equal to 0.5 take 0 (Down)
            random_list.append(0)
    return random_list

def random_generation():
    """
    
    Randomly creates three systems.

    """
    print "Three random systems generated: "
    print "First: " + str(System(randomize()))
    print "Second: " + str(System(randomize()))
    print "Third: " + str(System(randomize()))
    
def monte_carlo(steps=10):
    """
    
    Creates a random system and tries to decrease the energy with 
    Monte Carlo/Metropolis method with given steps. 

    """
    print "Step number: " + str(steps-1)
    sys1 = System(randomize()) # Create a random system
    sys2 = equals(sys1)
    sys2[randint(0,4)] = randint(0,2) # Make a random perturbation in system
    print "First :" + str(sys1)
    print "Second :" + str(sys2)
    for i in range(1,steps):
        if sys1.energy() > sys2.energy(): # Take first generated structure
            sys1 = equals(sys2)           # with lower energy
            print "Step " + str(i) + ": " + str(sys1)
        else:
            if -exp(sys2.energy()-sys1.energy()) >= rand(): # Metropolis criterion
                sys1 = equals(sys2)                         # -exp(deltaE/RT) >= random(0,1)
                print "Step " + str(i) + ": " + str(sys1)
            else:
                sys2[randint(0,4)] = randint(0,2) # Perturb the system again
                print "Step " + str(i) + ": " + str(sys1)

def genetic_algorithm():
    
    # Genetic operators
    def mutate(system):
        """
        
        Point mutation to system.

        """
        idx = randint(0,4) # Pick one of four spins randomly
        val = randint(0,2) # Randomly take up or down
        sys = equals(system)
        sys[idx] = val # Assign randomly taken spin to randomly chosen index
        return sys
        
    def cross_over(system1,system2):
        """
        
        Cross-over between two systems. 
        
        """
        sys1 = equals(system1)
        sys2 = equals(system2)
        sys1[randint(0,4)] = sys2[randint(0,4)]
        sys2[randint(0,4)] = sys1[randint(0,4)]
        return (sys1,sys2)
    
    population = [System(randomize()), System(randomize()), 
                  System(randomize()), System(randomize())]
                  
    pop_energies = [population[0].energy(),population[1].energy(),
                    population[2].energy(),population[3].energy()]
                  
    print "First randomly created population:\n"
    for system in population: print system
    print

    parents = []
    # Choosing the parents according to fittnes function
    # If random number is equal or above 0.4 take the system with minimum 
    # energy, if it is less than 0.4 take a random one with higher energy.
    # Choice is biased to choosing system with minimum energy. 
    for i in range(2):
        if rand() >= 0.4:
            parents.append(population[pop_energies.index(min(pop_energies))])
        else:
            parents.append(population[randint(0,len(population))])
        population.pop(pop_energies.index(min(pop_energies)))
        pop_energies.pop(pop_energies.index(min(pop_energies)))
    
    print "Parents:\n"
    for i in parents: print i
    print
    
    children = []
    # Create children from parents by crossing over parents' spins.
    for i in range(4):
        children.append(cross_over(parents[0],parents[1])[0])
    
    print "Children:\n"
    for i in children: print i
    print
    
    # Mutate two random children 
    ran1 = randint(0,4)
    ran2 = randint(0,4)
    children[ran1] = equals(mutate(children[ran1]))    
    children[ran2] = equals(mutate(children[ran2]))
    # Do cross-over for two random children
    ran3 = randint(0,4)
    ran4 = randint(0,4)
    children[ran3] = equals(cross_over(children[ran3],children[ran4])[0])
    children[ran4] = equals(cross_over(children[ran3],children[ran4])[1])
    
    print "Final population after two cross-overs and two mutations:\n"
    for i in children: print i
    print
    
    
    
    
if __name__ == "__main__":
    print "Random Generation"
    print "-------------------------------"
    random_generation()
    print 
    print "Monte Carlo/Metropolis Method"
    print "-------------------------------"
    monte_carlo()
    print
    print "Genetic Algorithm"
    print "-------------------------------"
    genetic_algorithm()