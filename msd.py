"""
Compute mean-squared displacement (MSD) averaged over entire trajectory
for diffusion simulations.

Contact: daniel.brooks@alumni.caltech.edu
"""

import numpy as np

def compute_msd(xyz, step_jump=1, start_frame=0):
  """
  Computes the MSD of a set of atoms in a trajectory.
  
  The MSD is computed from the displacements from all times in the 
  trajectory separated by times specific in the compute_MSD_frames vector.
 
  Parameters
  ----------
  xyz : numpy array
      Array containing xyz coordinates of atoms over which to compute MSD. 
      Array shape is of the form: xyz.shape = [n_frames, n_atoms, 3].
  step_jump: integer
      The MSD is computed for times separated by step_jump timesteps.
      Default value is 1. 
  start_frame: integer
      Ignore this many frames at the start of the trajectory.
      Default value is 0.

  Returns
  -------
  MSD_vector : numpy array
      The computed MSD at each of the frames specified in compute_MSD_frames.
  atom_MSD_vector : numpy array
      The atomic MSD at each of the frames specified in compute_MSD_frames.
  """
  
  #Determine the size of the trajectory. 
  end_frame = xyz.shape[0]
  n_atoms = xyz.shape[1]

  #Determine the time separations for MSD computation.
  min_step = 0
  max_step = end_frame - start_frame
  compute_MSD_frames = np.arange(min_step, max_step, step_jump) #MSD times.
  num_frames = compute_MSD_frames.shape[0]
  
  #Initialize vectors that store MSD. 
  MSD_vector = np.zeros(num_frames)
  atom_MSD_vector = np.zeros((n_atoms, num_frames))

  #Efficiently compute MSD over the specified time separations.
  mp = 0
  for MSD_step in compute_MSD_frames:
    
    #Specify all positions separated by a time MSD_step (vectorized).
    mat1 = xyz[np.arange(start_frame, end_frame - MSD_step), :,:]
    mat2 = xyz[np.arange(start_frame + MSD_step, end_frame), :,:]
    
    #Compute the squared displacements.
    MSD_mat = np.square(mat1-mat2)

    #Compute the averaged 3D MSD.
    ave_MSD = (np.sum(MSD_mat) / (MSD_mat.size*1.0)) * 3.0
    MSD_vector[mp] = ave_MSD

    #Compute the MSD for each atom.
    for a in range (0, n_atoms):
      cur_atom_MSD_mat = MSD_mat[:,a,:]
      ave_atom_MSD = (np.sum(cur_atom_MSD_mat)/(1.0*cur_atom_MSD_mat.size)) * 3.0
      atom_MSD_vector[a,mp] = ave_atom_MSD
    
    #Consider the next MSD time separation.
    mp = mp + 1
  
  #Return the average and atomic MSD vectors.
  return (MSD_vector, atom_MSD_vector)

def fit_d(MSD, times, frac_range=[0.25, 0.50]):
  """
  Fits D linearly to the given MSD over the fraction of the trajectory
  specified by fitting_range.

  Parameters
  ----------
  MSD : numpy array
      Contains MSD at a series times (nm^2).
  times : numpy array
      Times of each MSD measurement (ps).
  frac_range : [float lower, float upper]
      Optional: Perform fit between [lower, upper] fractions of trajectory.

  Returns
  -------
  D : float (cm^2/s)
  """

  #Define the fitting range.
  nf = len(MSD)
  frame_bound = [round(nf*frac_range[0]), round(nf*frac_range[0])]

  #Extract the relevant frames.
  MSD_fit = MSD[frame_bound[0]:frame_bound[1]]
  times_fit = times[frame_bound[0]:frame_bound[1]]

  #Perform the fit.
  coeff = np.polyfit(times_fit, MSD_fit, 1)
  D_nm2_ps = coeff[0]/6.0

  #Specify unit conversion factors to cm^2/s.
  nm2_2_cm2 = pow(10,14)
  ps_2_s = pow(10,12)
  unit_coeff = nm2_2_cm2/ps_2_s

  #Perform unit conversion and return D.
  D_cm2_s = D_nm2_ps * unit_coeff
  return D_cm2_s
