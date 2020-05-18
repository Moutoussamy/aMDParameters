#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Propose parmaters for accelerated molecular dynamics based of a classical molecular dynamics
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2020/05"

import sys
import pandas as pd
import numpy as np



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

def GetFinalStep(desire_time,timestep):
    """
    Get the step corresponding to the desire time
    e.g 5000000 if the desire time is 10ns and the timestep is 2fs

    :param desire_time: mk the average during this time
    :param timestep: the NAMD timestep
    :return: the step corresponding to the desire time
    """
    desire_time = desire_time*1000000
    final_step = desire_time/timestep

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

def CalculateParam(NrjMtx,timestep,desire_time):
    """
    Get avg. dhedral and total energy
    :param NrjMtx: energies matrix
    :param timestep: NAMD timestep
    :param desire_time: mk the average during this time
    :return: avg. dhedral and total energy
    """
    final_step = GetFinalStep(desire_time,timestep)
    SelectedValues = NrjMtx.loc[NrjMtx['step'] <= final_step]

    return np.mean(SelectedValues["Total"]), np.mean(SelectedValues["Dihedral"]) #AVG TOTAl, AVG DIHE

if __name__ == '__main__':

    NrjMtx = GetNRJ(sys.argv[1]) # collect energies from NAMD log file
    CalculateParam(NrjMtx,2,10)
    GetSystemInfo(sys.argv[2])

