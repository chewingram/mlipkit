import numpy as np
from ase.io import read
from ase.build import make_supercell
from pathlib import Path

from my_utils.utils_stdep import launch_stdep
from mlipkit.MTP import MPT_model

temperature = 200 
ucell = read(f'T200K_unitcell.json')
scell_mat = np.eye(3) * 4

pot_files = {'potential_file': '/gpfs/scratch/ehpc14/ulie583683/Work/ML/MoS2/MTP/mlip/mtp_6/training/pot.mtp'}
mlip_bin = '/Users/samuel/Work/codes/mlip-2/build1/mlp'
bin_pref = 'mpirun -n 6'

mlip_model = MTP_model(root_dir='pre_trained_model', 
                       name='MTP_model', 
                       pre_trained=True, 
                       pre_trained_pot_filepaths=pot_files)

launch_stdep(root_dir = './',
                ucell=ucell,
                make_supercell = True,
                scell_mat=scell_mat,
                scell=None,
                temperature = temperature,
                preexisting_ifcs = False,
                preexisting_ifcs_path=None,
                max_freq = 15,
                quantum = True,
                first_order = False,
                displ_threshold_firstorder = 0.0001,
                max_iterations_first_order = 20,
                rc2 = 20,
                niters = 20,
                nconfs = [4, 8, 10, 20, 40, 80, 100, 200]
                remove_infiles = True,
                mlip_bin = mlip_bin,
                bin_prefix = bin_pref,
                mlip_model=mlip_model)
