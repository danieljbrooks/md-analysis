"""
Simple script showing an example MSD calculation.

Contact: daniel.brooks@alumni.caltech.edu
"""

import mdtraj as md
import sys

from msd import compute_msd, fit_d

#Load the trajectory using mdtraj.
TRJ_FILE = "short.lammpstrj"
GEO_FILE = "structure.pdb"

#Compute and print the lithium MSD. 
li_ind = t.topology.select("element Li")
[msd_ave, msd_atomic] = compute_msd(t.xyz[:,li_ind,:])
print("Li MSD(t): ", msd_ave)

#Compute and print the diffusion coefficient.
D = fit_d(msd_ave, t.time)
print("D_Li (cm^2/s): ", D)
