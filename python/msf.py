# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 13:21:58 2014

@author: Serdar Ozsezen (c)

Lisence: BSD-3-Clause

An example of mean square fluctuation calculation
for a set of given MD snapshots in a zip file.

"""

import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from numpy import array, mean, sum, square, sqrt, zeros, flipud
import zipfile

def read_snapshots():
    """

    () -> list(file)
    
    Returns the snapshot pdb files as a list. 
    
    """
    snapshots = []    
    zipped = zipfile.ZipFile("tim_md.zip","r")
    for pdb in zipped.namelist()[1:]:
        snapshots.append(zipped.open(pdb))
    return snapshots

def fetch_alphacarbons(pdb_file):
    """
    
    (file) -> numpy.array
    
    Grabs alpha carbon coordinates and returns them as array.
    
    """
    carbons = []
    for line in pdb_file:
        if "ATOM  " in line and " CA " in line:
            split = line.split()
            carbons.append([float(split[5]),float(split[6]),float(split[7])])
    return array(carbons)

def mean_pos(coords):
    """
    
    (numpy.array) -> numpy.array
    
    Gives the mean positions of the atoms.
    
    """
    return mean(coords, axis=0)
    
def D_r(coords):
    """
    
    (numpy.array) -> numpy.array
    
    Returns deltaR_i for one atom.

    """
    return coords - mean_pos(coords)
    
def D_r2(coords):
    """
    
    (numpy.array) -> numpy.array    
    
    Returns the <deltaR_i^2> for one atom.    
    
    """
    return sum(square(D_r(coords)))
    
def array_alphacarbons():
    """
    
    () -> numpy.array
    
    Returns the coordinates of all alphacarbon atoms in all snapshots.
    
    """
    alphacarbons = []
    for model in read_snapshots():
        alphacarbons.append(fetch_alphacarbons(model))
    return array(alphacarbons)
    
def msf():
    """
    
    () -> numpy.array
    
    Returns the mean square fluctuations as an array. 
    
    """
    alphacarbons = []
    mean_sq_fluct = []
    for model in read_snapshots():
        alphacarbons.append(fetch_alphacarbons(model))
    cas = array(alphacarbons)
    for i in range(len(cas[0])):
        mean_sq_fluct.append(D_r2(cas[:,i])/10)
    return array(mean_sq_fluct)

def cross_correlation():
    """
    
    () -> numpy.ndarray
    
    Returns the cross correlation matrix as numpy.ndarray object. 
    
    """
    cas = array_alphacarbons()
    rng = len(cas[0])
    cross = zeros((rng,rng))
    fluct = msf()
    for i in range(rng):
        for j in range(rng):
            cross[i,j] = sum(D_r(cas[:,i])*D_r(cas[:,j]))*0.1\
            /sqrt(fluct[i])/sqrt(fluct[j])
    return cross
    
if __name__ == "__main__":
    
    # Calculate mean square fluctuations
    mean_square_fluctuations = msf()
    # Calculate cross correlations
    cross_correlation_map = flipud(cross_correlation())
    # Draw mean square fluctuations
    fig1, ax1 = plt.subplots()
    ax1.set_title("Mean square fluctuations for alphacarbons of TIM",
                  fontweight="bold",fontsize=20)
    ax1.plot(mean_square_fluctuations)
    loc = plticker.MultipleLocator(base=20)
    ax1.tick_params(labelsize=17)
    ax1.xaxis.set_major_locator(loc)    
    ax1.set_ylabel(u"Mean Square Fluctuations $(\u00c5^2)$", fontsize=20)
    ax1.set_xlabel(u"Residue Index", fontsize=20)
    ax1.grid()
    # Draw cross correlation map
    fig2, ax2 = plt.subplots()
    ax2.set_title("Cross correlation map for TIM",
              fontweight="bold",fontsize=20)
    img = ax2.imshow(cross_correlation_map,
                     extent=(0,len(cross_correlation_map),
                             0,len(cross_correlation_map)))
    ax2.tick_params(labelsize=17)
    loc2 = plticker.MultipleLocator(base=40)
    ax2.xaxis.set_major_locator(loc2)
    ax2.yaxis.set_major_locator(loc2)
    ax2.set_xlabel("Residue Index",fontsize=20)
    ax2.set_ylabel("Residue Index",fontsize=20)
    fig2.colorbar(img)
    