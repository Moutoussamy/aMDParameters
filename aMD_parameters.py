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



if __name__ == '__main__':
    GetNRJ(sys.argv[1]) # collect energies from NAMD log file