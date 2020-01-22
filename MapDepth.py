#!/usr/bin/env python3

import sys
import re
import os


"""
This script allow map the depth on a pdb structure.
The depth correspond to the distance between each CA atoms and the average phosphate plane
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2020/01"
__copyright__ = "CC_by_SA"
__dependencies__ = ""


def GetMembCenter(pdb):
	"""
	This function allow to calculate the membrane center
	:param pdb: a PDB file of a protein/membrane complex
	:return: the geometrical center of the membrane
	"""

	z_coord = []

	with open(pdb) as inputfile:
		for line in inputfile:
			if "ATOM" in line:
				if line[12:16].replace(" ","") == "P": # only the phosphate are taking into account
					z_coord.append(float(line[46:54]))
	return np.mean(z_coord)


def GetPhosPlane(pdb,memb):
	"""
	Calculation of the average upper phosphate plane
	:param pdb:  a PDB file of a protein/membrane complex
	:param memb: the center of geometry calculated by the function GetMembCenter
	:return: the average upper phosphate plane
	"""

	z_phosphate_upper_plane = []

	with open(pdb) as inputfile:
		for line in inputfile:
			if "ATOM" in line:
				if line[12:16].replace(" ","") == "P":
					if float(line[46:54]) > memb:
						z_phosphate_upper_plane.append(float(line[46:54]))

	return np.mean(z_phosphate_upper_plane)


def GetDists(Phos,pdb):
	"""
	Calculatoon of the CA - avg phosphate plane distance
	:param Phos: avg phosphate plane distance
	:param pdb: a PDB file of a protein/membrane complex
	:return: a dictionary with key= residue ID; value= dist.
	"""
	distDico = {}
	
	with open(pdb) as inputfile:
		for line in inputfile:
				if line[12:16].replace(" ","") == "CA":
					z_courant =  float(line[46:54]) - Phos 
					distDico[line[22:26]] = z_courant
	return distDico

def Mapped(pdb,dico):
	"""
	Map the dist. on a pdb file
	:param pdb: a PDB file of a protein/membrane complex
	:param dico: the dico created by the function GetDists
	:return: write a pdb file with the depth as the B-factor
	"""

	output = open("DepthMaped.pdb","w")
	with open(pdb) as inputfile:
		for line in inputfile:
			if "PROA" in line:
				resid = line[22:26]
				extend = line[46:].replace(line[60:66],"%6.2f"%dico[resid])
				line = "%s%s"%(line[0:46],extend)
				line = line.replace("HSD","HIS")
				line = line.replace("HSE","HIS")
				output.write(line)

if __name__ == '__main__':

	memb = GetMembCenter(sys.argv[1])
	Phos = GetPhosPlane(sys.argv[1],memb)	
	dists = GetDists(Phos,sys.argv[1])
	Mapped(sys.argv[1],dists)