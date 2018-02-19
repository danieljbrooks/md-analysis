"""
Simple script showing an example MSD calculation.

Contact: daniel.brooks@alumni.caltech.edu

"""

#Import the required packages and include scripts in path.
import mdtraj as md
import sys
sys.path.append("/ul/brooks/scripts/")
from msd import compute_msd

#Load the trajectory using mdtraj.
t = md.load("short.lammpstrj", top="structure.pdb")

#Compute and print the lithium MSD. 
li_ind = t.topology.select("element Li")
[ave, atomic] = compute_msd(t.xyz[:,li_ind,:])
print("Average MSD: ", ave)

