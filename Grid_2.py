import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import product
import seaborn as sns
import scipy as sc
import math
import pandas as pd
from glob import glob
from itertools import groupby
import pickle
# %matplotlib notebook
from PIL import Image
from tqdm import tqdm
import sys
from fnmatch import fnmatch
from matplotlib.colors import  LinearSegmentedColormap
import os

import Data_Tagger as dt

# TODO:

# Scale grid with box size
# refresh grids after second run
# combine xtcs


class GridBox(dt.Base_Obj):
    
   
    
    def __init__(self, structure_file, traj_file, init_spacing,
                 ligand_selection=None, mask=None, tag='', temperature=300.65):
        
        
        """
        Creates a GridBox object. This object creates 
    
    
        For the creation of an mda Universe:
        
        gro (str/list): structure file or list of structures
        xtc (str/list): trajectory file or list of trajectories
        (Doesn't necessairly have to be these file types)
        
        init_spacing (float or int): to form the grid, as many gripoints that can fit into the structure file
                        will be put in with this spacing. The spacing is fixed across the trajcetory, so
                        this code might not work properly if you have large fluctuations in box size
                        
        ligand_selection (str or callable): If a str, it will be passed to mda.Universe.select_atoms().
                                            Refer to the MDA manual on selection syntax. The atom selection will
                                            then be used as AtomGroup.residues. If you don't like that, pass a func instead
                                            If a function is passed it should return an MDA object with the attribute  
                                            atoms
                                            
        mask (str or callable): Pass a str or function that selects the which gridpoints to ignore.
                                This is done by turning the gridpoints into dummy atoms in the mda Universe.
                                The gripoints have the resname GRIDPOINT. GP selected by this, will not have 
                                a tally of whether a ligand atom is in proximity to the point. 
                                If one only wants a count of the GP explored around a protein, the could pass a str like:
                                \" resname GRIDPOINT and not around 5 protein \" 
                                
                                The function must return a string that is passed to u.select_atoms(). This allows
                                a string to be constructed based on properties found within that system; say excluding
                                the gridpoints below a relative coord
                                
                                
                                
        tag: An arbitrary value to store on the GridBox object. I use it to name GridBox's, 
             giving each obj a unique str or tuple to identify it.
                    
        temperature (float): The temperature to be used when the free energy at each GP is calculated based on the 
                                probability of a ligand atom occuping that GP.

        """
        
        
        self.structure_file = structure_file
        self.traj_file = traj_file
        self.init_spacing = init_spacing
        self.temperature = temperature
        self.ligand_selection = ligand_selection
        self.mask = mask
        # tag is just for object identification
        self.tag = tag
        # Attr to store tag information. This makes the inheritence work more smoothly
        self.attr_vals = {}
        
        
        
        # Create universe if xtc is not a list
        # else assume it is a list and pair gro with xtc (make sure order is correct)
        # universes will be made later
        if not isinstance(self.traj_file, list):
            self.u = mda.Universe(self.structure_file, self.traj_file)
            self.ligs = self._get_ligands(self.u, ligand_selection)
        else:
            self.structure_file.sort()
            self.traj_file.sort()
            self._gro_traj_pairs = list(zip(self.structure_file, self.traj_file))
            
            # just for setup purposes
            self.u = mda.Universe(self._gro_traj_pairs[0][0], self._gro_traj_pairs[0][1])
        
        
        
        self.npoints = self._get_npoints(self.u, self.init_spacing)      
        
        self._setup_grid()
        
        
        if self.mask != None:
            self._construct_mask()
            
        
    def _add_count(self, gp):
        # Work around for the issue where the grid does not scale with the box
        try:
            self.grid[gp[0], gp[1], gp[2]][-1] += 1
        except IndexError:
            pass
    
    # identify how to respond to the ligand_selection input
    def _get_ligands(self, u, ligand_selection):
        
        if ligand_selection == None:
            raise Exception('The ligand_selection kwarg must be parsed. It must be either a str or func')
            
        elif callable(ligand_selection):
            return ligand_selection(u)
        
        elif isinstance(ligand_selection, str):
            return u.select_atoms(ligand_selection).residues
        
        else:
            raise Exception('ligand selection was neither a str or func')
        
    def _construct_position_func(self, peratom, verbose):
        
        # determine whether to iterate through the atoms or the residues, taking the position or COM as needed
        if peratom == True:
            iter_group = self.ligs.atoms
            if verbose:
                print(f'iter group has {len(iter_group)} atoms')
            # wrapped in func so that the attr is callable
            position_func = lambda l: l.positions
       
        else:
            
            iter_group = self.ligs
            if verbose:
                print(f'itergroup has {len(iter_group)} residues')
            
            
            def position_func(ligs):
                pos = [lig.atoms.center_of_mass() for lig in ligs]
                return np.array(pos)
            
        return position_func, iter_group
                    
        
   
    
    def run(self, peratom=False, start=0, end=-1, step=1, grid_traj_step=0, verbose=False,):
        
        """
        Calculates the how frequently the grid points are occupied by the ligand_selection.
        A number of outputs are stored on the object (see below)
        
        Args, kwargs
        peratom (bool): If False, the center of mass of the ligand selection is used to consider the 
                        occupancy of a gridpoint
        start (int): Frame index of the trajectory to begin the calculation at.
        
        end (int): Frame index of the trajectory to end the calculation at. -1 for last frame of traj
        
        step (int): NUmber of frames of the trajectory to skip the calculation of. 
        
        grid_traj_step (int): Saves a copy of the grid every grid_traj_step steps. Useful for studying
                            convergence. This can very quickly consume to much memory if this value is
                            too low. If this value is not 0, the results will be stored as a list of np. arrays
        
        Returns: None
        The results are stored on the grid object. The results are stored in a number of different ways.
        Most are 3d np.arrays, where each element represents the corresponding gridpoint: 
                                                                                                                   
        self.coords_grid: 4d matrix. The first 3d each represent the grid, with the 4th D, its real coords
        self.counts: 3d matrix. elements are the raw counts for each gridpoint.
        self.zcounts: Same as counts except nanvalues are 0
        self.p_grid: 3d: values are the probabilty of gridpoint being occupied
        self.fe_grid: 3d: free energy value of the gridpoint being occupied by the ligand.  
        
        
        """
        
        
        
        
        if isinstance(self.traj_file, list):
            self._combined_xtc_analysis(peratom=peratom, start=start, end=end, step=step)
            return None
        
        
        
        
        self.n_explored = []
        self.times = []
        
        # Saves the grid every grid_traj_step 
        if grid_traj_step > 0:
            self.grid_traj = []
        
        
        self.peratom = peratom
        
        position_func, iter_group = self._construct_position_func(self.peratom, verbose)
        
        
        if grid_traj_step > 0:
            
            def grid_step():
                if ts.frame % grid_traj_step == 0 or ts.frame == self.u.trajectory[-1].frame:
                    self.grid_traj.append(np.copy(self.grid))
        else:
            def grid_step():
                pass
        
        
        # Loop through frames 
        print('Starting grid calculation ... ')
        for ts in tqdm(self.u.trajectory[start:end:step]):
            
            
            positions = position_func(iter_group)
            
            gps = self._find_closest_grid_point(positions)
            
            for gp in gps:
                self._add_count(gp)
            
        
            self.times.append(ts.time)
            
            
            # I have commented this out as it is unnecessary and time consuming if I am not
            # actually needing the n_epxlored vs time
            # nan is counted as nonzero, so must be changed to zero for the count
#             nonzero_grid = np.where((np.isnan(self.grid[:,:,:,-1])) | (self.grid[:,:,:,-1] == 0),
#                                     0, self.grid[:,:,:,-1])
#             self.n_explored.append(np.count_nonzero(nonzero_grid))
            
            grid_step()
            
        
        # calc probabilities, energies, remove zeroes for vis, etc.
        print('Postprocessing... calculating final probabilites, energies, and a grid for visualisaiton ')
        
        if grid_traj_step == 0:
        # save outputs to self
            self.coords_grid, self.counts, self.zcounts, self.p_grid, self.fe_grid = self._post_processing(self.grid)
        else:
            self.coords_grid = []
            self.counts = []
            self.zcounts = []
            self.p_grid = []
            self.fe_grid = []        
            
            for grid in self.grid_traj:
                coords_grid, counts, zcounts, p_grid, fe_grid = self._post_processing(grid)
                
                self.coords_grid.append(coords_grid)
                self.counts.append(counts)
                self.zcounts.append(zcounts)
                self.p_grid.append(p_grid)
                self.fe_grid.append(fe_grid)
    
    
    
    # modified version of run() but for combinding results across multiple traj
    def _combined_xtc_analysis(self, peratom=False, start=0, end=-1, step=1,):
        
        self.n_explored = [0]  
        self.times = []
        
        self.peratom = peratom
        

        
        # Loop through frames 
        print('Starting grid calculation ... ')
        
        for gro, xtc in self._gro_traj_pairs:
            
            u = mda.Universe(gro, xtc)
            
            self.ligs = self._get_ligands(u, self.ligand_selection)
            # determine whether to iterate through the atoms or the residues, taking the position or COM as needed
            
            if self.peratom == True:
                iter_group = self.ligs.atoms
                # wrapped in func so that the attr is callable
                position_func = lambda l: l.positions
           
            else:
                
                iter_group = self.ligs
                
                def position_func(ligs):
                    pos = []
                    for lig in ligs:
                        pos.append(lig.atoms.center_of_mass())
                    return np.array(pos)
        
        
            # Find sim length and make sure that the n_explored will not give index error
            sim_len = len(u.trajectory)
            diff = sim_len - len(self.n_explored)
            if diff > 0:
                self.n_explored = self.n_explored + [0] * diff
                
        
            print(f'Working on {xtc} ...')
        
            for ts in u.trajectory[start:end:step]:
                
                
                positions = position_func(iter_group)
                
                gps = self._find_closest_grid_point(positions)
                
                for gp in gps:
                    self._add_count(gp)
                
    #             for lig in iter_group:
    #                 pos = position_func(lig)
    #                 gp = self._find_closest_grid_point(pos)
    #                 self._add_count(gp)
                if verbose:
                    if np.isnan(gp):
                        print(f'gp is maksed. time: {ts.time}')
                self.times.append(ts.time)
                
                # nan is counted as nonzero, so must be changed to zero for the count
                nonzero_grid = np.where((np.isnan(self.grid[:,:,:,-1])) | (self.grid[:,:,:,-1] == 0),
                                        0, self.grid[:,:,:,-1])
                self.n_explored[ts.frame] += np.count_nonzero(nonzero_grid)
            
        # calc probabilities, energies, remove zeroes for vis, etc.
        print('Postprocessing... calculating final probabilites, energies, and a grid for visualisaiton ')
        self._post_processing(self.grid)
    
    

        
        
        
    
    
    
    # Grid points that are found in the a mda selection will be converted from 0 to nan
    # This stops them from counts acruing
    # i.e. all points found by passed function will not be counted
    def _construct_mask(self,):
        
        # get real coords from gp indices
        # create list of gp coords
        # create dummy universe and add GRIDPOINT res to it and give it gp coords
        # merge with the self.u
        # Create a selection of the GRIDPOINT atoms / gp. THese will be the only gp that are counted
        
        self.gp_coords_list = self.gp_list * self.init_spacing
        ngps = len(self.gp_coords_list)
        udum = mda.Universe.empty(ngps, n_residues=ngps,  trajectory=True )
        udum.add_TopologyAttr('type', ['GRIDPOINT'] * ngps)
        udum.add_TopologyAttr('resname', ['GRIDPOINT'] * ngps)
        udum.add_TopologyAttr('name', ['GRIDPOINT'] * ngps)
        udum.add_TopologyAttr('mass', [1.0] * ngps)
        udum.atoms.positions = self.gp_coords_list
        
        
        
        self.umerge = mda.Merge(self.u.atoms, udum.atoms)
         
        print('Finding masked grid points ...')
        

        if callable(self.mask):
            self.selected = self.umerge.select_atoms(self.mask(self.umerge))
            print(f'selection found  {len(self.selected)} atoms')
            self.masked_coords = self.selected.positions
        else:
            self.masked_coords = self.umerge.select_atoms(self.mask).positions
            

        self.masked_gp = (self.masked_coords / self.init_spacing).astype(int)
        print(f'Found {len(self.masked_coords)} grid points to be excluded')
        print('Excluding ... ')
        
        print(len(self.masked_coords))
        # adding to nan will still remain a nan
        # This avoids putting an if statement in the traj loop 

        for gp in self.masked_gp:
            self.grid[gp[0], gp[1], gp[2]][-1] = np.nan
        
        print('Congratulations, the mask has been setup')
 
    # Deprecated version as basing grid size on furthest atoms was problematic if 
    # water and memb has been removed
#     def _get_npoints(self, u, spacing):
        
#         mns = u.atoms.positions.min(axis=0)
#         mxs = u.atoms.positions.max(axis=0)
        
#         boxdim = mxs - mns
#         np_real = boxdim / spacing
        
#         # Round every num up, to avoid any index errors
#         diff = np_real - np_real.astype(int)
#         add_arr = np.where(diff < 0.5, 1, 0)
#         np_preround = np_real + add_arr
#         npoints = np.around(np_preround).astype(int)
        
#         return list(npoints)
    
    
    def _get_npoints(self, u, spacing):
        
        boxdim = self.u.trajectory[0].dimensions[:3]
        
        np_real = boxdim / spacing
        
        # Round every num up, to avoid any index errors
        diff = np_real - np_real.astype(int)
        add_arr = np.where(diff < 0.5, 1, 0)
        np_preround = np_real + add_arr
        npoints = np.around(np_preround).astype(int)
        
        return list(npoints)
    
    
    

    # Given a point in space and a particular grid, which gridpoint is that point closest to?
    def _find_closest_grid_point(self, point ):
        
        # find how many grid points there are between point and 0
        ng_real = point / self.init_spacing
        
        # round down to nearest grid point
        ng_int = ng_real.astype(int)
    
        # find difference to know whether point is closer to the lower or upper gridpoint
        remainder = ng_real - ng_int
        # add 0 if lower, or grid_width if upper
        add_arr = np.where(remainder > 0.5, 1, 0)
        closest_grid_index = ng_int + add_arr
        
        return closest_grid_index
    
    def _setup_grid(self,):
        
        # setup grid where index is equal to its value
        # gplist has shape (npoints multiplied by eachother) by 3
        # It is  a list of indices essentially
        self.gp_list = np.array(list(product(np.arange(0, self.npoints[0], 1),
                np.arange(0, self.npoints[1], 1),
                np.arange(0, self.npoints[2], 1))))
        
        # init the count list
        counts_list = np.zeros([len(self.gp_list), 1])
        
        
        # add the counts to the grip points. Each i in the list will be [x,y,z,count]
        gp_counts_list = np.hstack((self.gp_list, counts_list))
        
        # reshape into a 3d grid, where the index is equal to its value in the grid
        self.grid = gp_counts_list.reshape(self.npoints + [4])
        
        return self.grid
    
    def _post_processing(self, grid ):
        
        
        coords_grid = grid[:,:,:,:3] * self.init_spacing
        counts = grid[:,:,:,3]
        counts[counts == 0] = np.nan
        
        zcounts = np.copy(counts)
        zcounts = np.nan_to_num(zcounts, nan=0)
        p_grid = zcounts / np.sum(zcounts)
        fe_grid =  (- sc.constants.k * self.temperature * np.log(p_grid)) / 1000 * sc.constants.Avogadro
        fe_grid[fe_grid == np.inf] = np.nan
    
        return coords_grid, counts, zcounts, p_grid, fe_grid
    
    
    
    def save_nexplored(self, outputname):
        f = open(outputname)
        
        #header info
        f.write('# Number of gridpoints explored vs time')
        f.write('# Produced with the following files:')
        f.write(f'# {self.structure_file} {self.traj_file}')
        f.write(f'# {self.tag} ')
        
        for t, n in zip(self.times, self.n_explored):
            f.write(f'{t} {n}')
        f.close()
        
    def output_grid(self, filename, attr, exclude_val):
        
        
        #flatten grids into 1d array
#         coords_list = self.coords_grid.reshape(np.prod(self.coords_grid.shape[:3]), 3)
        coords_list = self.gp_list
        property_list = getattr(self, attr).flatten()
                              
        print(coords_list.shape)
        print(property_list.shape)
        
        # Adjust values for viewablility. I dunno why the values are weird. Check up on this later
        property_list = property_list * (10 ** 23)
        # A to nm
        coords_list = coords_list * 10
        
        self.coords_with_prop = []
        self.props = []
        for i, gp in enumerate(coords_list):
            if np.isnan(property_list[i]):
                continue
            else:
                self.coords_with_prop.append(gp)
                self.props.append(property_list[i])
        
        print(len(self.coords_with_prop))
        print(len(self.props))
        
        u = mda.Universe.empty(
            len(self.coords_with_prop),
            trajectory=True        
        )
        
        
        u.add_TopologyAttr('name', ['H']*len(self.coords_with_prop))                 
        u.atoms.positions = self.coords_with_prop
                                          
        u.add_TopologyAttr('tempfactors')
        print(self.props[:5])
        u.atoms.tempfactors = self.props
        u.atoms.write(filename)
        
        
# make it not system specific
#  Should be able to find an object based on a tuple
class GridManagerOLD(dt.Base_Manager):
    
    def __init__(self, scheme, dirs, gro_str, xtc_str, resolution, prot_present=False, res_selection=None, ligand_selection=None):
        
        """ The way of detecting what the ligand selection should be is shoddy but does work. I've created
        a new version of the gridmanager that does it a better way but is a little more restricted in the systems it can upload
        """


        self.objs = []
        self.attr = []
        self.file_attr = {}
        
        if prot_present:
            with open('isoform_natoms.txt', 'r') as f:
                iso_natoms = {}
                for line in f.readlines():
                    i, na = line.split()
                    iso_natoms[i] = na
                
        
 
        for d in dirs:
        
            print(f'uploading: {d}')
            # Find the necessary structure and traj files for this grid 
            if d[-1] == scheme.slash:
                end = ""
            else:
                end = scheme.slash
            
            gro = d + end + gro_str
            xtc = d + end + xtc_str
            
            
            # Determine how the ligand will be selected
            # if protein is in file
            attr_vals = scheme._pair_attr_vals(d+'.')
            if prot_present:
                iso = attr_vals['iso']
                na_ind = iso_natoms[iso]
                ligand_selection = f'index {na_ind}:10000'
            
            # If there is no protein, assume all atoms are ligand atoms
            # assumes there are less than 1000 lig atoms
            # !!!!! Not a good solution if there is more than prot and ligand in the system !!!!!!
            else:
                ligand_selection = 'index 0:10000'
            
            if res_selection != None:
                ligand_selection=res_selection
            
            
            if ligand_selection != None:
                ligand_selection = ligand_selection
            
            # Create the grid
            try: 
                grid = GridBox(gro, xtc, resolution,
                               ligand_selection=ligand_selection,
                               tag=tuple(attr_vals.values()) )
            
                # The grid box obj now inherits from dt.Base_obj, so I can call the _label method
                # to give the grid obj the attr specified in the scheme
                grid._label(**attr_vals)
                self.objs.append(grid)
            except FileNotFoundError:
                print(f'File not found for {d}')
                continue
        
    
    
    def run(self, **kwargs):
        
        for grid in self.objs:
            print(grid)
            grid.run(**kwargs)


class GridManager(dt.Base_Manager):
        
    def __init__(self, scheme, dirs, gro_str, xtc_str, resolution, ):



        self.objs = []
        self.attr = []
        self.file_attr = {}
        
        if 'ligand' not in scheme.attr_list:
            raise Exception('The gridmanager requires that there be an attribute called \'ligand\' in the attribute list of the scheme. \n This tells the GridManager how to make the ligand selection for grid calculation')



        for d in dirs:
        
            print(f'uploading: {d}')
            
            
            gro = os.path.join(d, gro_str)
            xtc = os.path.join(d, xtc_str)
            
            # Get attributes about this system from the file name according to scheme object
            attr_vals = scheme._pair_attr_vals(d+'.')
            
            ligand_selection = f'resname {attr_vals["ligand"]}'
            
            # Create the grid
            try: 
                grid = GridBox(gro, xtc, resolution,
                               ligand_selection=ligand_selection,
                               tag=tuple(attr_vals.values()) )
            
                # The grid box obj now inherits from dt.Base_obj, so I can call the _label method
                # to give the grid obj the attr specified in the scheme
                grid._label(**attr_vals)
                self.objs.append(grid)
                print(f"The ligand selection will be: \'{ligand_selection}\'")


            except FileNotFoundError:
                print(f'File not found for {d}')
                continue
        
    
    
    def run(self, **kwargs):
        
        for grid in self.objs:
            print(grid)
            grid.run(**kwargs)




        
def assess_convergence(grid_objs, selection1, selection2, steps, sim_length=4500, **kwargs):
    """
    This function calculates the convergence of a group of trajectories by measuring the average energy of
    a two selected regions of the grid across time. If the relative difference between these two 
    regions is consistent, then the simulation is converged (at least in respect to those regions).
    For my system this will be a average FE of a slice of the bulk and the average FE of the cavity.
    
    Args: 
    grid_objs (list): list of GridBox objects which should already be initialised. 
                        each grid should be a single replicate if you want ot calculuate convergence across time
    
    selection1 (tuple): a tuple of 6 integers which correspond to the start and end selection 
                        ranges for dimensions x,y, and z of the grid
    
    selection2 (tuple): same as above
                    
    steps (float): The number of frames across trajectories to assess
    
    Kwargs: will be passed to GridBox.run()
    
    """
    
    # initialise raw results array. stores av fe of each selection, for every step, for every grid
    # +1 for the last frame that is always inlcuded (now adding 100 to try and avoid index error. graphing should ignore extra nan anyway)
    n_datapoints = int(sim_length / steps) + 100
    energies = np.full(( n_datapoints, len(grid_objs), 2), np.nan)
    
    
    for g, gridobj in enumerate(grid_objs):
        print(gridobj.tag, g)
        # Run the grid, saving every number of steps
        gridobj.run(grid_traj_step=steps, **kwargs)
        
        # For each of the steps, make the two selections
        # and calculate the average FE of each
        for i in range(len(gridobj.grid_traj)):
        
            
            s1, s2 = selection1, selection2
            
            avfe1 = np.nanmean(gridobj.fe_grid[i][s1[0]:s1[1], s1[2]:s1[3], s1[4]:s1[5],])
            energies[i, g, 0] = avfe1
            
            avfe2 = np.nanmean(gridobj.fe_grid[i][s2[0]:s2[1], s2[2]:s2[3], s2[4]:s2[5],])
            energies[i, g, 1] = avfe2
            
            print(i, g, avfe1, avfe2)
            
    
             
    return energies
        

    ################################
############################
# initialise raw results array. stores av fe of each selection, for every step, for every grid
    # 
    n_datapoints = int(sim_length / steps)
    energies = np.full(( n_datapoints, len(grid_objs.objs), 2), np.nan)
    
    
    for g, gridobj in enumerate(grid_objs):
        print(gridobj.tag, g)
        # Run the grid, saving every number of steps
        gridobj.run(grid_traj_step=steps, **kwargs)
        
        # For each of the steps, make the two selections
        # and calculate the average FE of each
        for i in range(len(gridobj.grid_traj)):
        
            
            s1, s2 = selection1, selection2
            
            avfe1 = np.nanmean(gridobj.fe_grid[i][s1[0]:s1[1], s1[2]:s1[3], s1[4]:s1[5],])
            energies[i, g, 0] = avfe1
            
            avfe2 = np.nanmean(gridobj.fe_grid[i][s2[0]:s2[1], s2[2]:s2[3], s2[4]:s2[5],])
            energies[i, g, 1] = avfe2
            
            print(i, g, avfe1, avfe2)
            
    
             
    return energies
####################################
##################################
        

    
    
if __name__ == '__main__':

    glob_str = sys.argv[1].replace('\\', '')
    outputfile = sys.argv[2]


    floodrun_scheme = dt.Scheme(['_', 'iso', 'memb', 'ligand', '_', 'rep'], split_list=['rep'])


    all_grids = GridManagerOLD(floodrun_scheme, sorted(glob(glob_str)),
                 'mdeq_Prot-DRUG-only.gro', 'trjcat_centered-Prot_Prot-DRUG-only.xtc', 
                1, prot_present=True)
    
    all_grids.objs.sort(key = lambda o: (o.iso, o.ligand))

    selection1 = (0, 91, 0, 91, 85, 90 )
    selection2 = (30, 70, 30,70, 40,70)
    skipnframes = 250

    for key, grids in groupby(all_grids.objs, key = lambda o: (o.iso[:3], o.ligand)):
        print(key)
        grids=list(grids)
        energies = assess_convergence(grids, selection1, selection2, skipnframes)

        # np.save(f'Convergence-test/Convergence_ediff_{grids[0].iso[:3]}_{grids[1].ligand}.npy', energies, )
        np.save(outputfile, energies,)
    completed_run_str = f"Finished calculating the convergence energies for {glob_str}"
    print('')
    print('#'*20)
    print(completed_run_str)
    print('#'*20)
    print('')


