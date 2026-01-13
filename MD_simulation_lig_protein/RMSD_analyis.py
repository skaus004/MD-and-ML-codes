# Cell 2: Imports
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys

# Suppress non-critical warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Cell 3: Load the trajectory and topology files dynamically
ID = sys.argv[1]
print(ID)
topology_file = f'./{ID}/IF_chai4_amber-{ID}-minimized.pdb'
trajectory_file = f'./{ID}/IF_chai4_amber-{ID}-output.dcd'

u = mda.Universe(topology_file, trajectory_file)
# Cell 4: Select atoms
backbone = u.select_atoms("protein and backbone")
unk_atoms = u.select_atoms("resname UNK")

if len(unk_atoms) == 0:
    raise ValueError("No atoms with resname UNK found.")
# Cell 5: Align trajectory to backbone
# Create a reference universe using the first frame
ref = mda.Universe(topology_file)
aligner = align.AlignTraj(u, ref, select="protein and backbone", in_memory=True)
aligner.run()
# Cell 6: RMSD calculation for aligned UNK
reference_positions = unk_atoms.positions.copy()
rmsd_values = []

for ts in u.trajectory:
    current = unk_atoms.positions
    diff = current - reference_positions
    rmsd = np.sqrt((diff ** 2).sum(axis=1).mean())
    rmsd_values.append(rmsd)

time = np.arange(len(rmsd_values)) * u.trajectory.dt / 1000
import numpy as np

# Exclude the first 10 RMSD values
rmsd_excluded = rmsd_values[100:]

# Calculate mean and standard deviation
rmsd_mean = np.mean(rmsd_excluded)
rmsd_std = np.std(rmsd_excluded)

print(f"Mean RMSD (excluding first 100 frames): {rmsd_mean:.4f}")
print(f"Standard Deviation: {rmsd_std:.4f}")

# Time array in nanoseconds
time = np.arange(len(rmsd_values)) * u.trajectory.dt / 1000  # Assuming dt is in ps, so time is in ns
