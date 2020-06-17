#!/usr/bin/env python3

import sys
import argparse
import numpy as np


"""
This script allow map the depth on a pdb structure.
The depth correspond to the distance between each CA atoms and the average phosphate plane
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2020/01"
__copyright__ = "CC_by_SA"
__dependencies__ = ""



def GetArgs():
	"""
	Get argument for the calculation using flags
	:return: argument: log file and pdb
	"""

	parser = argparse.ArgumentParser(description='Map the depth on a PDB file. The depth, here, correspond to\
	the distance between the upper phosphate (COM*) plane and the protein (COM).')

	parser.add_argument('-pdb', metavar="pdb", type=str, help= "pdb file")
	parser.add_argument('-Up', action='store_true', help="Does it bound to the Upperleaflet?")
	parser.add_argument('-Low',action='store_true', help="Does it bound to the lowerleaflet?")

	args = parser.parse_args()

	return args


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
					z_coord.append(abs(float(line[46:54])))

	return np.mean(z_coord)


def GetProtCenter(pdb):
	"""
	This function allow to calculate the Prot. center
	:param pdb: a PDB file of a protein/membrane complex
	:return: the geometrical center of the prot.
	"""

	z_coord = []

	with open(pdb) as inputfile:
		for line in inputfile:
			if "ATOM" in line:
				if line[12:16].replace(" ","") == "CA": # only the CA are taking into account
					print(line)
					z_coord.append(abs(float(line[46:54])))
	return np.mean(z_coord)

def GetPhosPlane(pdb,memb,up,low):
	"""
	Calculation of the average upper phosphate plane
	:param pdb:  a PDB file of a protein/membrane complex
	:param memb: the center of geometry calculated by the function GetMembCenter
	:param Prot: the center of geometry calculated by the function GetprotCenter
	:return: the average upper phosphate plane
	"""

	z_phosphate_upper_plane = []

	with open(pdb) as inputfile:
		for line in inputfile:
			if "ATOM" in line:
				if line[12:16].replace(" ","") == "P":

					if up :
						if float(line[46:54]) > memb:
							z_phosphate_upper_plane.append(abs(float(line[46:54])))
					if low:
						if float(line[46:54]) < memb:
							z_phosphate_upper_plane.append(abs(float(line[46:54])))

	return np.mean(z_phosphate_upper_plane)


def GetDists(Phos,pdb):
	"""
	Calculatoon of the CA - avg phosphate plane distance
	:param Phos: avg phosphate plane distance
	:return: a dictionary with key= residue ID; value= dist.
	"""
	distDico = {}

	with open(pdb) as inputfile:
		for line in inputfile:
				if line[12:16].replace(" ","") == "CA":
					z_courant =  (abs(float(line[46:54])) - Phos)
					distDico[int(line[22:26])] = z_courant
	return distDico

def Mapped(pdb,dico):
	"""
	Map the dist. on a pdb file
	:param pdb: a PDB file of a protein/membrane complex
	:param dico: the dico created by the function GetDists
	:return: write a pdb file with the depth as the B-factor
	"""

	pdbout = pdb.replace(".pdb","_mapped.pdb")
	output = open(pdbout,"w")

	with open(pdb) as inputfile:
		for line in inputfile:
			if "PROA" in line:
				resid = int(line[22:26])
				extend = line[46:].replace(line[60:66],"%6.2f"%dico[resid])
				line = "%s%s"%(line[0:46],extend)
				line = line.replace("HSD","HIS")
				line = line.replace("HSE","HIS")
				output.write(line)

def writeDist(dico,pdbname):
	"""
	Write the calculated depth in a file
	:param dico: the dico created by the function GetDists
	:param pdbname: the name of the input pdb
	:return: Write the depth in a txt file
	"""

	txt_file_name = pdbname.replace(".pdb","_depth.txt")

	output = open(txt_file_name,"w")

	for key in dico.keys():
		output.write("{0}\t{1}\n".format(key,dico[key]))

	output.close()

if __name__ == '__main__':
	args = GetArgs()  # collect arguments

	if args.Low == 0 and args.Up == 0:
		sys.exit('Please specify the leaflet? (-Up or -Low)')

	elif args.Low == 1 and args.Up == 1:
		sys.exit('The calculation cannot be done on both leaflets')

	memb = GetMembCenter(args.pdb)
	#prot = GetProtCenter(sys.argv[1])
	Phos = GetPhosPlane(args.pdb,memb,args.Up,args.Low)
	dists = GetDists(Phos,args.pdb)
	Mapped(args.pdb,dists)
	writeDist(dists,args.pdb)
