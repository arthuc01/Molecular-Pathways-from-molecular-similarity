# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:52:04 2015

@author: chxcja

Can we reconstruct a biosynthetic pathway simply from the structures. 
It turns out the answer is no - but you can get large contiguous parts of the
pathway - but there needs to be chemical insight into identifying unrealistic 
chemical or metabolic transformations.

This model uses some of the structures from the pseudomonic acid pathway
see cdx file
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw
import networkx as nx
G=nx.Graph()
 
mols = Chem.SmilesMolSupplier('pseudomonicAcid-pathway-smiles.txt')
print len(mols)

#set up nodes from each molecule
G.add_nodes_from([0,1,2,3,4,5,6])

#fps = [FingerprintMols.FingerprintMol(x) for x in mols]

#Morgan fingerprints work better - radius doesn't seem to make much diff.
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in mols]

for mol in range(len(mols)-1):
    #fig = Draw.MolToMPL(mols[mol])
    max2Start=0
    max2Temp=0
    maxStart=0
    maxTemp=0
    for mol2 in range(mol, len(mols)):

        if mol != mol2:
            temp = DataStructs.FingerprintSimilarity(fps[mol],fps[mol2])
            if temp>maxTemp: 
                max2Start = maxTemp
                max2Temp = maxStart
                maxTemp =  temp
                maxStart = mol2
            print mol, mol2, temp
    #Create node between most similar molecules
    G.add_edge(mol, maxStart,weight = maxTemp, len=(1-maxTemp))
    #G.add_edge(mol, max2Start,weight = max2Temp, len=(1-max2Temp))
    print "**", mol, maxStart, maxTemp
    
import matplotlib.pyplot as plt
nx.draw(G, prog = "neato", with_labels=True)

print(nx.shortest_path(G,source=0,target=6))
