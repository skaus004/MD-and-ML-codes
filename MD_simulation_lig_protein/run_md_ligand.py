#!/usr/bin/env python

from openmm import *
from openmm.app import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology as OffTopology
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem
import argparse
import sys
import os
from tempfile import NamedTemporaryFile


def run_simulation(pdb_filename, ligand_filename, sim_time_ns):
    """
    Performs a molecular dynamics simulation on a given PDB file.

    Args:
        pdb_filename (str): The path to the input PDB file.
        ligand_filename (str): The path to the ligand file (SDF or MOL2 format).
        sim_time_ns (float): The total simulation time in nanoseconds.
    """

    # Set the md_prefix from the pdb_filename
    md_prefix = os.path.splitext(os.path.basename(pdb_filename))[0] + '-' \
        + os.path.splitext(os.path.basename(ligand_filename))[0] + '-'
    sim_time = sim_time_ns * nanosecond
    step_size = 2 * femtosecond
    steps = round(sim_time / step_size)

    # Load the ligand molecule based on the postfix of the ligand_filename
    if ligand_filename.endswith('.sdf'):
        ligand_molecule = Molecule.from_file(ligand_filename)
    elif ligand_filename.endswith('.mol2'):
        ligand_molecule_rdkit = Chem.MolFromMol2File(ligand_filename)
        if ligand_molecule_rdkit is None:
            raise ValueError(f"Could not read the ligand file: {ligand_filename}")
        ligand_molecule = Molecule.from_rdkit(ligand_molecule_rdkit)
    else:
        raise ValueError(f"Unsupported ligand file format: {ligand_filename}")
    gaff_generator = GAFFTemplateGenerator(molecules=ligand_molecule)
    # Extract the ligand tology and positions and convert to OpenMM
    ligand_topology_off = OffTopology.from_molecules([ligand_molecule])
    ligand_topology_openmm = ligand_topology_off.to_openmm()
    ligand_positions = ligand_molecule.conformers[0].to_openmm()

    # Load the PDB file for the protein
    pdb = PDBFile(pdb_filename)

    
    # Specify the forcefield
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    forcefield.registerTemplateGenerator(gaff_generator.generator)

    # Modeller for building the system
    modeller = Modeller(pdb.topology, pdb.positions)

    # Add the ligand to the system
    modeller.add(ligand_topology_openmm, ligand_positions)

    # Delete water if present
    modeller.deleteWater()
    # Add hydrogens to the protein
    modeller.addHydrogens(forcefield)

    # Add solvent and ions
    modeller.addSolvent(forcefield, padding=1.0 * nanometer)

    # Create the system
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1.0 * nanometer, constraints=HBonds)

    #Integrator
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, step_size)

    # Simulation
    platform = Platform.getPlatform('CUDA')
    properties = {'Precision': 'mixed'}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimize energy
    print("Minimizing energy")
    simulation.minimizeEnergy()

    # Write the minimized structure to a PDB file
    PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(),
                        open(md_prefix + 'minimized.pdb', 'w'), keepIds=True)
    
    
    
    # Reporters
    simulation.reporters.append(DCDReporter(md_prefix + 'output.dcd', 2000))
    simulation.reporters.append(StateDataReporter(sys.stdout, 10000, step=True,
                                                potentialEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(StateDataReporter(md_prefix + "md_log.txt", 2000, step=True,
                                                potentialEnergy=True, temperature=True, volume=True))

    # NVT equilibration
    print("Running NVT")
    simulation.step(20000)

    # NPT production
    system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin))
    simulation.context.reinitialize(preserveState=True)

    print("Running NPT")
    simulation.step(steps)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run OpenMM MD simulation.')
    parser.add_argument('pdb_file', type=str, help='Input PDB file.')
    parser.add_argument('ligand_file', type=str, help='Input ligand file (SDF or MOL2).')
    parser.add_argument('sim_time', type=float, help='Simulation time in nanoseconds.')

    args = parser.parse_args()

    run_simulation(args.pdb_file, args.ligand_file, args.sim_time)