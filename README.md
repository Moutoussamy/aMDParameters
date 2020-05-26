# aMDParameters
Give parameters for accelerated molecular dynamics (https://www.ks.uiuc.edu/Research/namd/2.9/ug/node63.html).
Calculate the parameters for accelerated molecular dynamics based on a classical MD simulation (cMD).


- Inputs:
    - a PDB file of your system (only the protein, no water)
    - a log file from a classical MD

The input parameters Edihed, adihed, Etotal, atotal is calculated from  10 ns cMD simulation.

- Usage:
python aMD_parameters.py -log [NAMDLogFile] -pdb [PDBFile] -first [FISRT STEP, Default = 0]

- Example output:

Only dihedral:

    accelMD on
    accelMDdihe on
    accelMDE 12541.139773
    accelMDalpha 199.500000
    

Dual Boost:

    accelMD on
    accelMDdihe on
    accelMDE 12541.139773
    accelMDalpha 199.500000
    accelMDdual on
    accelMDTE -70549.059176
    accelMDTalpha 397.075000
