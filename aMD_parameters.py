#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Propose parameters for accelerated molecular dynamics based of a classical molecular dynamics
aMD documentation: https://www.ks.uiuc.edu/Research/namd/2.9/ug/node63.html
A NAMD log file of cMD and a pdb is necessary.
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2020/05"

import pandas as pd
import numpy as np
import argparse
import sys

def GetArgs():
    """
    Get argument for the calculation using flags
    :return: argument: log file and pdb
    """

    parser = argparse.ArgumentParser(description='Propose parameters for accelerated molecular dynamics based of a \
    classical MD')

    parser.add_argument('-log', metavar = "log", type = str, help = "log file from a classical MD")
    parser.add_argument('-pdb', metavar="pdb", type=str, help= "pdb file, no water")
    parser.add_argument('-first', metavar="first",default=0, type=int, help="step where to start the average")

    args = parser.parse_args()

    return args


def CheckPDBExtension(pdb):
    """
    Check extension of PDB file
    :param pdb: the argument given with the -pdb flag
    :return: 0 if the extension is not '.pdb'
    """
    if pdb[-4:] == ".pdb":
        return 1
    else:
        return 0

def CheckPDBAtom(pdb):
    """
    Check if the PDB file contains ATOM line
    :param pdb: a pdb file
    :return: 0 if no atom in pdb file else return 1
    """
    flag = 0
    with open(pdb) as input_file:
        for line in input_file:
            if line[0:4] == "ATOM":
                flag = 1
    return flag


def CheckLOGExtesnion(log):
    """
    Check extension of the log file
    :param log: the argument given with the -log flag
    :return: 0 if the extension is not '.log'
    """
    if log[-4:] == ".log":
        return 1
    else:
        return 0

def CheckLOGEnergy(log):
    """
    Check extension of the log file
    :param log: the argument given with the -log flag
    :return: 0 if the log does not contain energies
    """
    flag = 0
    with open(log) as input_file:
        for line in input_file:
            if line[0:7] == "ENERGY:":
                flag = 1
    return flag



def CheckArguments(arguments):
    """
    Check if the correct argument has been given
    :param arguments: list of argument given
    :return: error message if something wrong (FILE extension or incorrect data)
    """
    if not CheckPDBExtension(arguments.pdb):
        sys.exit('ERROR: argument for -pdb is not a PDB file')
    if not CheckPDBAtom(arguments.pdb):
        sys.exit('ERROR: no atom in the PDB file')
    if not CheckLOGExtesnion(arguments.log):
        sys.exit('ERROR: argument for -log is not a LOG file')
    if not CheckLOGEnergy(arguments.log):
        sys.exit('ERROR: no energies in the LOG file')

def GetNRJ(NamdLogFile):
    """
    Get total and dihedral energies for NAMD log file
    :param NamdLogFile: NAMD log file
    :return: numpy matrix containing the energies
    """


    with open(NamdLogFile) as inputfile:
        for line in inputfile:
            if line[0:7] == "ENERGY:":
                line = line.split()

                if line[1] == "0":
                    matrix = [int(line[1]),float(line[11]),float(line[4])] #[step,total,dihe]

                else:
                    tmp = [int(line[1]),float(line[11]),float(line[4])] #[step,total,dihe]
                    matrix = np.vstack((matrix,tmp))

    names = ["step","Total","Dihedral"]
    matrix = pd.DataFrame(matrix,columns = names)

    return matrix

def GetFinalStep(desire_time,timestep, first_step):
    """
    Get the step corresponding to the desire time
    e.g 5000000 if the desire time is 10ns and the timestep is 2fs

    :param desire_time: mk the average during this time
    :param timestep: the NAMD timestep
    :param first_step: desire first step
    :return: the step corresponding to the desire time
    """
    desire_time = desire_time*1000000
    final_step = desire_time/timestep + first_step

    return final_step


def GetSystemInfo(pdb):
    """
    Get nb of atoms and nb of residues
    :param pdb: a pdb file
    :return: nb of atoms and nb of residues
    """
    residues = []
    nbatoms = 0

    with open(pdb) as inputfile:
        for line in inputfile:
            if line[0:4] == "ATOM":
                nbatoms += 1
                if line[22:26] not in residues:
                    residues.append(line[22:26])

    return nbatoms,len(residues)

def CalculateParam(NrjMtx,timestep,desire_time, first_step):
    """
    Get avg. dhedral and total energy
    :param NrjMtx: energies matrix
    :param timestep: NAMD timestep
    :param desire_time: mk the average during this time
    :return: avg. dhedral and total energies
    """
    final_step = GetFinalStep(desire_time,timestep,first_step)

    SelectedValues = NrjMtx.loc[first_step <= NrjMtx['step']]
    SelectedValues = SelectedValues.loc[SelectedValues['step'] <= final_step]

    return np.mean(SelectedValues["Total"]), np.mean(SelectedValues["Dihedral"]) #AVG TOTAl, AVG DIHE

def DihedralParam(avgDihe,nbRes):
    """
    Calculate E and alpha
    :param avgDihe: avg. dhedral energy
    :param nbRes: number of residues
    :return: E and alpha
    """
    EDihedral = avgDihe + (3.5*nbRes)
    alphaDihedral = (3.5*nbRes)/5

    return EDihedral,alphaDihedral

def TotalParam(avgTotal,nbatoms):
    """
    Calculate E and alpha for total energies
    :param avgTotal: avg. total energy
    :param nbatoms: number of atoms
    :return: E and alpha
    """
    ETotal = avgTotal + (0.175*nbatoms)
    alphaTotal = 0.175*nbatoms

    return ETotal,alphaTotal

def PrintCommandLine(avgTotal,avgDihedral,pdb):
    """
    Print param. for boost
    :param avgTotal: avg. total energy
    :param avgDihedral: avg. dihe. energy
    :param pdb: pdb file
    :return: Print param. for boost
    """
    nbatoms, nbres = GetSystemInfo(pdb)
    Etotal, alphaTotal = TotalParam(avgTotal, nbatoms)
    EDihedral, alphaDihedral = DihedralParam(avgDihedral, nbres)


    print("\nOnly dihedral:\n")
    print"""
    accelMD on
    accelMDdihe on
    accelMDE %f
    accelMDalpha %f
    accelMDOutFreq 5000
    """%(EDihedral,alphaDihedral)


    print("\nDual Boost:\n")
    print"""
    accelMD on
    accelMDdihe on
    accelMDE %f
    accelMDalpha %f
    accelMDdual on
    accelMDTE %f
    accelMDTalpha %f
    accelMDOutFreq 5000
    """%(EDihedral,alphaDihedral,Etotal,alphaTotal)



if __name__ == '__main__':
    args = GetArgs() #collect arguments
    CheckArguments(args) # Check arguments
    NrjMtx = GetNRJ(args.log) # collect energies from NAMD log file
    avgTotal, avgDihedral = CalculateParam(NrjMtx,2,10,args.first) # Caculate the parameters for aMD
    PrintCommandLine(avgTotal, avgDihedral,args.pdb) # Print the command lines for NAMD

