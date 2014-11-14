# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 10:34:37 2014

@author: Serdar

For parsing XML files of sequence identity clusters of PDB.org

"""

from urllib import urlretrieve
from xml.dom import minidom


def cluster_lists(PDB,chain):
    """

    (string,string) -> dict
    
    Grabs the chain names from PDB sequence identity clusters and returns them 
    as a dictionary, sorted according to sequence identity cluster percentage.     

    """
    clusters = { 100:[],90:[],70:[],50:[],40:[],30:[] }
    keys = sorted(clusters.keys(),reverse=True)
    lst = []
    for percentage in clusters:
        urlretrieve("http://www.pdb.org/pdb/rest/sequenceCluster?cluster=" + 
        str(percentage) + "&structureId="+PDB[0:4]+"." + chain,"cluster" +
        str(percentage) + ".xml")
        xmldoc = minidom.parse("cluster" + str(percentage) + ".xml")
        for i in xmldoc.getElementsByTagName('pdbChain'):
            clusters[percentage].append(i.attributes['name'].value)
    
    for i in range(len(keys)):
        if keys[i] != 100:
            lst += clusters[keys[i-1]]
            clusters[keys[i]] = list(set(clusters[keys[i]]) - set(lst))
        
    return clusters

def write_to_file(file_name,clusters):
    """

    (string,dict) -> NoneType
    
    Writes the names of the chains in the clusters into a text file. Starts 
    from highest sequence identity. 
    
    """
    percentages = [100,90,70,50,40,30]   
    with open(file_name,"w") as pdb_list:
        for i in percentages:
            print "----------------------"
            print "Sequence identity: " + str(i) + "\n"
            for pdb in clusters[i]:
                print pdb
            pdb_list.writelines("\n".join(clusters[i]))
            pdb_list.write("\n")
    
    
if __name__ == "__main__":
    a = cluster_lists("8ohm","A")
    write_to_file("clusters_8ohm.txt",a)