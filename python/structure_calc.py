# -*- coding: utf-8 -*-
"""
Created on Tue Mar 08 14:41:42 2014

@author: Serdar Ozsezen (c)

Lisence: BSD-3-Clause

Various structure related calculations 

"""
import matplotlib.pyplot as plt
from numpy import array, dot, zeros, mean
from numpy.linalg import norm
from math import acos, atan2, degrees, pow, sqrt

""" MA is the list of masses of atoms """

MA = {"C":float(12),"O":float(16),"N":float(14),"S":float(32)}

def read_pdb(pdb_name):
    """ 
    
    (string) -> file
    
    Reads the PDB file from disk. PDB file should be in the same 
    working directory with the script.
    
    """
    return open(pdb_name,"r")
    
def fetch_atoms(pdb_file,chain_name):
    """ 
    
    (file,string)-> dictionary(int: tuple(string,string,numpy.array))
    
    Grabs all the atoms from PDB file with given chain name. Saves all the
    atoms according to their numbers in PDB file with their coordinates,
    names and residue names as shown in example:
    
    65: (3,"ALA","CA", array([5.4, 2.3, 20.6]))
    
    """
    atoms = {}
    chain = " " + chain_name + " "    
    for line in pdb_file:
        if "ATOM  " in line and chain in line:
            split = line.split()
            atoms[int(split[1])] = (int(split[5]),split[3],split[2],\
            array([float(split[6]),float(split[7]),float(split[8])]))
    return atoms

def fetch_res_names(atoms):
    """

    (file,string) -> list(string)

    Returns the residue names as a list of strings. 
    
    """
    names = []
    for val in atoms:
        if atoms[val][2] == "CA":
            names.append(atoms[val][1])
    return names

def fetch_backbone(atoms):
    """
    
    (dictionary) -> tuple(list(numpy.array),list(tuple))
    
    Grabs the coordinates of backbone atoms (N, CA, C) and assigns the 
    coordinates to list elements.
    
    """
    backbone = []
    backbone_names = []
    for val in atoms:
        if atoms[val][2] == "N" or atoms[val][2] == "CA" or atoms[val][2] == "C":
            backbone.append(atoms[val][3])
            backbone_names.append((atoms[val][0],atoms[val][1],atoms[val][2]))
    return (backbone, backbone_names)
    
def fetch_alphacarbons(atoms):
    """
    (dictionary) -> tuple(list(numpy.array),list(tuple))
    
    Grabs the coordinates of alpha carbons and assigns the coordinates to 
    list elements.

    """
    carbons = []
    carbons_names = []
    for val in atoms:
        if atoms[val][2] == "CA":
            carbons.append(atoms[val][3])
            carbons_names.append((atoms[val][0],atoms[val][1],atoms[val][2]))
    return (carbons,carbons_names)

def centroids(atoms):
    """
    
    (dictionary) -> list(numpy.array)
    
    Finds the coordinates of centroids of the side chains. 
    
    """
    residues = fetch_res_names(atoms)
    side_chain ={}
    centroid_list =[]
    for i in range(len(residues)):
        for val in atoms:
            if atoms[val][1] == residues[i] and atoms[val][0] == i+1 and \
            atoms[val][2] != "N" and atoms[val][2] != "CA" and \
            atoms[val][2] != "C" and atoms[val][2] != "O":
                side_chain[val] = atoms[val]
        centroid_list.append(center_of_mass(side_chain))
        side_chain.clear()                        
    return centroid_list
    
def pseudo_length(atoms):
    """
    
    (dictionary) -> list(float)
    
    Calculates the pseudobond lengths between side chain centroids and 
    alpha carbons.
    
    """
    alphacarbons = fetch_alphacarbons(atoms)[0]
    centroid_list = centroids(atoms)
    bonds = []    
    for i in range(len(alphacarbons)):
        bonds.append(norm(centroid_list[i]-alphacarbons[i]))
    return bonds

def bond_lengths(backbone):
    """
    
    (list(numpy.array)) -> list(float)
        
    Calculates the bond lengths from coordinates of the atoms and returns them 
    as a list of floats. 
    
    """
    bondlengths = []   
    for i in range(len(backbone)):
        if i+1 < len(backbone):
            bondlengths.append(norm(backbone[i+1]-backbone[i]))
    return bondlengths

def bond_angles(backbone):
    """
    
    (list(numpy.array)) -> list(float)
    
    Calculates the bond angles in degrees and returns them as a list of floats.
    
    """
    bondangles = []
    for i in range(len(backbone)):
        if i+2 < len(backbone):
            v1 = backbone[i] - backbone[i+1]
            v2 = backbone[i+2] - backbone[i+1]
            bondangles.append(degrees(acos(dot(v2,v1)/norm(v1)/norm(v2))))
    return bondangles
    
def cross(vector1, vector2):
    """
    
    (numpy.array,numpy.array) -> numpy.array
    
    Returns the cross products of two vectors as a NumPy array. 
    
    """
    return array([vector1[1]*vector2[2]-vector1[2]*vector2[1],\
    vector1[2]*vector2[0]-vector1[0]*vector2[2],vector1[0]*vector2[1]\
    -vector1[1]*vector2[0]])
    
def dihedral_angles(backbone):
    """
    
    (list(numpy.array)) -> list(float)
    
    Calculates and returns the torsional angles as a list of floats. 
    
    """
    dihedralangles = []
    for i in range(len(backbone)):
        if i+3 < len(backbone):
            v1 = backbone[i+1] - backbone[i]
            v2 = backbone[i+2] - backbone[i+1]
            v3 = backbone[i+3] - backbone[i+2]
            angle = atan2(dot(cross(cross(v1,v2),cross(v2,v3)),(v2/norm(v2))),\
            dot(cross(v1,v2),cross(v2,v3)))
            dihedralangles.append(degrees(angle))
    return dihedralangles

def side_chain_angles(atoms):
    """
    
    (dictionary) -> list(float)
    
    Calculates the angles between alpha carbons and side chain centroids. 
    
    
    """
    angles = []
    carbons = fetch_alphacarbons(atoms)[0]
    sidechains = centroids(atoms)
    for i in range(len(carbons)):
        if i+1 < len(carbons):
            v1 = carbons[i] - carbons[i+1]
            v2 = sidechains[i+1] - carbons[i+1]
            angle = degrees(acos(dot(v1,v2)/norm(v1)/norm(v2)))
        angles.append(angle)
    return angles
            
    
    
def center_of_mass(atoms):
    """

    (dictionary) -> numpy.array    
    
    Calculates the center of mass coordinates of a chain or side chain.
    
    """
    coords_mass = array([0,0,0]); total_mass = 0.0;
    for val in atoms:
        if "C" in atoms[val][2]:
            coords_mass += atoms[val][3]*MA["C"]
            total_mass += MA["C"]
        elif "O" in atoms[val][2]:
            coords_mass += atoms[val][3]*MA["O"]
            total_mass += MA["O"]
        elif "N" in atoms[val][2]:
            coords_mass += atoms[val][3]*MA["N"]
            total_mass += MA["N"]
        elif "S" in atoms[val][2]:
            coords_mass += atoms[val][3]*MA["S"]
            total_mass += MA["S"]
    if total_mass == 0.0:
        center_of_mass = array([nan, nan, nan])
    else:
        center_of_mass = array([coords_mass[0]/total_mass,\
        coords_mass[1]/total_mass,coords_mass[2]/total_mass])
    return center_of_mass

def radius_of_gyration(atoms,center_of_mass):
    """

    (dictionary,array) -> float
    
    Calculates the radius of gyration by using formula below:
    
    Rg^2 = sum(|r_i - Rc|^2)/N
    
    where r_i is the coordinates of ith atom, Rc is center of mass and N is the
    number of atoms in the chain. 
    
    """
    diff = zeros(shape=(len(atoms)))
    for val in atoms:
        diff[val-1] = pow(norm(atoms[val][3] - center_of_mass),2)
    return sqrt(mean(diff))

def R_end(alphacarbons):
    """
    
    (dictionary) -> numpy.array

    Calculates the end to end vector.

    """
    return alphacarbons[0] - alphacarbons[-1]
   
def print_bond_angles(pdb_code, names, angles):
    """

    (string,list,list) -> NoneType
    
    Prints the bond angles by giving the names of relevant atoms.
    
    """
    print("Bond angles of " + pdb_code + ":\n")
    print("{:>17s}{:>26s}".format("Atoms","Angles"))
    print("-"*50)
    s = " "
    for i in range(len(names)):
        if i+2 < len(names):
            string = str(names[i][0])+s+names[i][1]+s+names[i][2]+s+"-"+s+\
            str(names[i+1][0])+s+names[i+1][1]+s+names[i+1][2]+s+"-"+s+\
            str(names[i+2][0])+s+names[i+2][1]+s+names[i+2][2]
            angle = angles[i]
            print("{:<37s}{a:5.2f}".format(string,a=angle))
    print "\n"
    
def print_torsion_angles(pdb_code, names, angles):
    """
    
    (string,list,list) -> NoneType
    
    Prints the torsional angles by giving the names of relevant atoms.
    
    """
    print("Torsion angles of " + pdb_code + ":\n")
    print("{:>21s}{:>35s}".format("Atoms","Angles"))
    print("-"*60)
    s = " "
    for i in range(len(names)):
        if i+3 < len(names):
            string = str(names[i][0])+s+names[i][1]+s+names[i][2]+s+"-"+s+\
            str(names[i+1][0])+s+names[i+1][1]+s+names[i+1][2]+s+"-"+s+\
            str(names[i+2][0])+s+names[i+2][1]+s+names[i+2][2]+s+"-"+s+\
            str(names[i+3][0])+s+names[i+3][1]+s+names[i+3][2]
            angle = angles[i]
            print("{:<50s}{a:5.2f}".format(string,a=angle))
    print "\n"
    
def print_sidechain_angles(pdb_code, names, angles):
    """
    
    (string,list,list) -> NoneType
    
    Prints the angles between alphacarbons and side chain centroids.
    
    """
    print("Angles between alpha-carbons and side chain centroids "+\
    "for " + pdb_code + ":")
    print "(SCC = Side Chain Centroid)\n"
    print("{:>25s}{:>21s}".format("Interaction Centers","Angles"))
    print("-"*50)
    s = " "
    for i in range(len(names)):
        if i+1 < len(names):
            string  = str(names[i][0])+s+names[i][1]+s+names[i][2]+s+"-"+s+\
            str(names[i+1][0])+s+names[i+1][1]+s+names[i+1][2]+s+"-"+s+\
            str(names[i+1][0])+s+names[i+1][1]+s+"SCC"
            angle = angles[i]
            print("{:>34s}      {a:5.2f}".format(string,a=angle))
    print "\n"

def print_sidechain_lengths(pdb_code, names, lengths):
    """
    
    (string,list,numpy.array) -> NoneType
    
    Prints the lengths of the pseudobonds between side chain centroids and 
    relevant alphacarbon atoms. 
    
    """
    s = " "
    print("Pseudobond length between alphacarbons and "\
    "side chain centroids for "+ pdb_code + " are:")
    print "(SCC = Side Chain Centroid)\n"
    print("{:>11s}{:>14s}".format("Interaction Centers",u"Lengths(\u00c5)"))
    print("-"*32)
    for i in range(len(names)):
        string = string = str(names[i][0])+s+names[i][1]+s+names[i][2]+s+"-"+s+\
        str(names[i][0])+s+names[i][1]+s+"SCC"
        length = lengths[i]
        print("{:<24s} {l:5.3f}".format(string,l=length))
    print "\n"
    
def print_bond_lengths(pdb_code, names, lengths):
    """
    
    (string,list,numpy.array) -> NoneType
    
    Prints the torsional angles by giving the names of relevant atoms.
    
    """
    s = " "
    print("Bond lengths of " + pdb_code + ":\n")
    print("{:>11s}{:>21s}".format("Atoms",u"Lengths(\u00c5)"))
    print("-"*32)
    for i in range(len(names)):
        if i+1 < len(names):
            string = string = str(names[i][0])+s+names[i][1]+s+names[i][2]+s+"-"+s+\
            str(names[i+1][0])+s+names[i+1][1]+s+names[i+1][2]
            length = lengths[i]
            print("{:<24s}{l:5.3f}".format(string,l=length))
    print "\n"
    
def calculate(pdb_name,chain_name):
    """
    (string,string) -> dictionary
    
    Calculates all relevant properties (available in this script) for 
    given chain.
    
    """
    
    PDB = read_pdb(pdb_name+".pdb")
    atoms_PDB = fetch_atoms(PDB,chain_name)
    bb_PDB, names_bb_PDB = fetch_backbone(atoms_PDB)
    ca_PDB, names_ca_PDB = fetch_alphacarbons(atoms_PDB)
    bond_PDB = bond_lengths(bb_PDB)
    ba_PDB = bond_angles(bb_PDB)
    da_PDB = dihedral_angles(bb_PDB)
    sc_bonds_PDB = pseudo_length(atoms_PDB)
    sc_angles_PDB = side_chain_angles(atoms_PDB)
    bond_red_PDB = bond_lengths(ca_PDB)    
    ba_red_PDB = bond_angles(ca_PDB)
    da_red_PDB = dihedral_angles(ca_PDB)
    com_PDB = center_of_mass(atoms_PDB)
    S = radius_of_gyration(atoms_PDB,com_PDB)    
    Re = R_end(ca_PDB)
    return {"backbone":bb_PDB, "atom names": names_bb_PDB,
            "alphacarbons":ca_PDB, "residues": names_ca_PDB,
            "bond lengths": bond_PDB, "bond angles": ba_PDB,
            "torsion angles": da_PDB, "side chain bond lengths": sc_bonds_PDB,
            "side chain bond angles": sc_angles_PDB, 
            "virtual bond lengths": bond_red_PDB,
            "virtual bond angles": ba_red_PDB,
            "CA torsion angles": da_red_PDB,
            "center of mass": com_PDB, "radius of gyration": S,
            "end to end vector": Re}
            
def print_all(pdb,chain_name):
    """
    
    (string,string) -> NoneType
    
    Calculates and prints all relevant properties (available in this script) for 
    given chain.

    """
    res = calculate(pdb,chain_name)
    print "*** Calculations for all atom chain model for "+ pdb +" ***\n\n"
    print_bond_lengths(pdb,res["atom names"],res["bond lengths"])
    print_bond_angles(pdb,res["atom names"],res["bond angles"])
    print_torsion_angles(pdb,res["atom names"],res["torsion angles"])
    print "*** Calculations for virtual bond model for " + pdb + " ***\n\n"    
    print_sidechain_lengths(pdb,res["residues"],res["side chain bond lengths"])
    print_bond_lengths(pdb, res["residues"],res["virtual bond lengths"])    
    print_sidechain_angles(pdb,res["residues"],res["side chain bond angles"])
    print_bond_angles(pdb+" reduced model", res["residues"],res["virtual bond angles"])
    print_torsion_angles(pdb+" reduced model",res["residues"],res["CA torsion angles"])
    Re = res["end to end vector"]
    print "End to end vector: " + str(Re[0])+"i "+str(Re[1])+"j "+ \
    str(Re[2])+"k"
    print "Radius of gyration: " + str(res["radius of gyration"])

def draw(pdb,results):
    """
    
    (string,dictionary) -> NoneType
    
    Draws the graphics for calculation results.
    
    """
    da = results["torsion angles"]
    ba = results["bond angles"]
    bl = results["bond lengths"]
    da_red = results["CA torsion angles"]
    ba_red = results["virtual bond angles"]
    ba_sc = results["side chain bond angles"]
    bl_sc = results["side chain bond lengths"]
    bl_red = results["virtual bond lengths"]
    
    plt.figure(1)
    plt.subplot(211)
    plt.title("Torsion angles and bond angles of "+pdb,fontweight="bold",fontsize=15)
    plt.plot(da,"ro")
    plt.xlabel("Index",fontsize=15); plt.ylabel("Torsion Angles (Degrees)",fontsize=15)
    plt.grid()
    plt.subplot(212)
    plt.plot(ba,"go")
    plt.xlabel("Index",fontsize=15); plt.ylabel("Bond Angles (Degrees)",fontsize=15)
    plt.grid()
    
    plt.figure(2)
    plt.title("Bond lengths of backbone atoms of "+pdb,fontweight="bold",fontsize=15)
    plt.plot(bl,"ro")
    plt.xlabel("Index",fontsize=15); plt.ylabel(u"Bond Lengths (\u00c5)",fontsize=15)
    plt.grid()
    
    plt.figure(3)
    plt.subplot(211)
    plt.title("Torsion angles and bond angles of "+pdb+" for reduced model"\
    ,fontweight="bold",fontsize=15)
    plt.plot(da_red,"ro")
    plt.xlabel("Index",fontsize=15); plt.ylabel("Torsion Angles (Degrees)",fontsize=15)
    plt.grid()
    plt.subplot(212)
    plt.plot(ba_red, "go")
    plt.xlabel("Index",fontsize=15); plt.ylabel("Bond Angles (Degrees)",fontsize=15)
    plt.grid()
    
    plt.figure(4)
    plt.title("Side chain bond angles of "+pdb+" for reduced model"\
    ,fontweight="bold",fontsize=15)
    plt.plot(ba_sc,"ro")
    plt.xlabel("Index",fontsize=15); plt.ylabel("Bond Angles (Degrees)",fontsize=15)
    plt.grid()
    
    plt.figure(5)
    plt.subplot(211)
    plt.title("Bond lengths of "+pdb+" for reduced model"\
    ,fontweight="bold",fontsize=15)
    plt.plot(bl_sc,"ro")
    plt.xlabel("Index",fontsize=15); plt.ylabel(u"Bond Lengths (\u00c5)",fontsize=15)
    plt.grid()
    plt.subplot(212)
    plt.plot(bl_red,"go")
    plt.xlabel("Index",fontsize=15); plt.ylabel(u"Bond Lengths (\u00c5)",fontsize=15)
    plt.grid()
    
if __name__ == "__main__":
    res1 = calculate("1AKE","A")
    
    draw("1AKE", res1)
    