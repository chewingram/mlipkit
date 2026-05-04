import ase
from ase.io import read, write
from ase.io.trajectory import TrajectoryWriter, TrajectoryReader, Trajectory
from ase.md.langevin import Langevin
from ase.md.nptberendsen import NPTBerendsen
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import numpy as np
from mattersim.forcefield import MatterSimCalculator
from copy import deepcopy as cp
from pathlib import Path

calc = MatterSimCalculator(load_path='/gpfs/scratch/ehpc14/ulie583683/mattersim_pth/mattersim-v1.0.0-1M.pth', device='cpu')
conf = read('init_structure.traj')
conf.calc = calc

T1 = 2000
T2 = 11000
nsteps = 7000
step = (T2-T1)/nsteps
ts = 1 # in fs
temps = np.arange(T1, T2, step)
temps = np.concatenate((temps, np.ones(10000)*T2))

MaxwellBoltzmannDistribution(conf, temperature_K=T1)
conf.calc = calc

trajectory = []
trajpath = Path('./Trajectory.traj')
logpath = Path('./md_log')

if trajpath.is_file():
    trajpath.unlink()

if logpath.is_file():
    logpath.unlink()
    
for i, T in enumerate(temps):
    traj = Trajectory(trajpath.absolute(), 'a')
    simulation = NPTBerendsen(conf,
                          timestep=ts * ase.units.fs,
                          temperature_K=T,
                          pressure_au=0,
                          taut=100*units.fs,
                          taup=1000*units.fs,
                          compressibility_au=4.57e-5 / units.bar,
                          logfile=logpath.absolute()
                         )
    simulation.run(1)
    traj.write(cp(conf))
    traj.close
    trajectory.append(cp(conf))
    
write('Trajectory.traj', trajectory)
    
                          
                          
            
		      
