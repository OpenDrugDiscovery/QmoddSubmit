import mdtraj
import openmm
import datamol as dm
import numpy as np

from loguru import logger
from functools import partial
from typing import Union, Optional

import openmm.unit as unit
from openmm import app, State
from openmm.unit import Quantity
from openff.toolkit.utils.exceptions import (
    UnassignedValenceParameterException,
    UnassignedProperTorsionParameterException,
)
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology

from qmodd_submit.utils import flatten



def generate_conformers(
    smi: str,
    ff_engine: ForceField,
    n_conformers: int = 50,
    n_starting_points: int = 10,
    n_conformers_per_starting_point: int = 100,
    n_conformers_after_rmsd_filter: int = 25,
) -> Union[Molecule, None]:
    """
    Generate the conformations for a molecule given a smiles string.

    Parameters
    ----------
    smi: str
        The smiles string of the molecule
    ff_engine: str
        The force field engine to use for the simulation
    n_conformers: int
        The total number of conformers to generate
    n_starting_points: int
        The number of starting points to use for the MD simulation
    n_conformers_per_starting_point: int 
        The number of conformers to generate from each starting point
    n_conformers_after_rmsd_filter: int
        The number of conformers to keep after the RMSD filter

    Returns
    -------
    mol: openff.toolkit.topology.Molecule
        The molecule with the generated conformers

    Source: https://github.com/openmm/spice-dataset/blob/main/pubchem/createPubchem.py#L25-L83
    """

    if n_conformers % n_conformers_after_rmsd_filter != 0:
        raise ValueError("`n_conformers` should be a multiple of `n_conformers_after_rmsd_filter`")

    # Parameterize the system
    logger.debug("Parameterizing the simulation system")
    mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
    ret = mol_to_openmm_topology_and_system(mol, ff_engine)
    # print(mol, ret)

    if ret is None:
        return

    topology, system = ret

    # Generate diverse starting points for MD simulation.
    # Run MD from each one to generate the total set of higher energy (due to the high tempature) conformations.
    # NOTE (Cas): We need the rms_cutoff set to 0, as the None would result in the non-zero default
    logger.debug(f"Generating {n_starting_points} starting points")
    mol.generate_conformers(n_conformers=n_starting_points, rms_cutoff=None)

    fn = partial(
        simulate_low_energy_conformers,
        topology=topology,
        system=system,
        max_conformers=n_conformers_per_starting_point,
    )
    logger.debug(
        f"Running a MD simulation to find {n_starting_points * n_conformers_per_starting_point} conformers"
    )
    all_conformers = dm.utils.parallelized(
        fn,
        # NOTE: Conversion is needed for now as openff and openmm are in the
        #  middle of a refactoring to a unified system.
        [m for m in mol.conformers],
        progress=True,
        tqdm_kwargs={"desc": "MD simulations", "leave": False},
    )
    all_conformers = flatten(all_conformers)
    logger.info(f"{len(all_conformers)}, {n_conformers_after_rmsd_filter}")
    # Select a subset that is maximally different from one another
    if len(all_conformers) < n_conformers_after_rmsd_filter:
        msg = f"Filtered to generate at least {n_conformers_after_rmsd_filter} low-energy conformers"
        logger.warning(msg)
        return

    logger.debug(
        f"Filtering conformers by their RMSD to end up with {n_conformers_after_rmsd_filter} conformers"
    )
    all_conformers = filter_by_rmsd(all_conformers, topology)

    n = (n_conformers // n_conformers_after_rmsd_filter) - 1
    logger.info(f"{len(all_conformers)}, {n}")

    if n >= 1:
        # Create a nearby, lower energy conformation from each conformer we have so far.
        fn = partial(
            simulate_low_energy_conformers,
            topology=topology,
            system=system,
            max_conformers=n,
            temperature=100 * unit.kelvin,
            max_energy_minimization_iters=5,
            n_simulation_steps=1000,
            energy_threshold=None,
        )
        logger.debug(
            f"Generating an additional {n} conformer(s) for each of the {n_conformers_after_rmsd_filter} "
            f"conformers we have selected so far."
        )
        new_conformers = dm.utils.parallelized(
            fn, all_conformers, progress=True, tqdm_kwargs={"desc": "MD simulations", "leave": False}
        )
        new_conformers = flatten(new_conformers)

        all_conformers += new_conformers

    assert len(all_conformers) == n_conformers, f"We don't have {n_conformers} conformers as expected"

    # Set all the generated conformers for the mol
    mol._conformers = None
    for state in all_conformers:
        mol.add_conformer(state.getPositions(asNumpy=True))

    return mol


def filter_by_rmsd(conformer_states, topology, n_conformers: int = 25):
    """
    Reduce the total set of conformers.

    Parameters
    ----------
    conformer_states: list of simtk.openmm.app.Simulation
        The set of conformers to filter
    topology: simtk.openmm.app.Topology
        The topology of the molecule
    n_conformers: int
        The number of conformers to keep
    Source: https://github.com/openmm/spice-dataset/blob/main/pubchem/createPubchem.py#L11-L23
    NOTE (Cas): Could be replaced by similar functionality in Datamol to remove the mdtraj dependency
    """

    xyz = np.array([s.getPositions().value_in_unit(unit.nanometer) for s in conformer_states])

    traj = mdtraj.Trajectory(xyz, mdtraj.Topology.from_openmm(topology))
    traj.center_coordinates()

    final_states = {0}
    min_rmsd = mdtraj.rmsd(traj, traj, 0, precentered=True)

    for i in range(n_conformers - 1):
        best = np.argmax(min_rmsd)
        min_rmsd = np.minimum(min_rmsd, mdtraj.rmsd(traj, traj, best, precentered=True))
        final_states.add(best)

    return [conformer_states[i] for i in final_states]


def mol_to_openmm_topology_and_system(mol: Union[str, Molecule], ff_engine: ForceField):
    """
    Initiates a MD simulation system given a molecule and a ForceField.

    Parameters
    ----------
    mol: Union[str, Molecule]
        The molecule to simulate. If a string is provided, it is assumed to be a SMILES string and the user is for defining the stereochemistry.
    ff_engine: openff.toolkit.typing.engines.smirnoff.ForceField
        The force field to use for the simulation.

    Returns
    -------
    openmm_topology: openmm.app.Topology
        The OpenMM topology for the molecule.
    system: openmm.System
        The OpenMM system for the molecule.

    Source: https://github.com/openmm/spice-dataset/blob/main/pubchem/createPubchem.py#L28-L36
    """
    if isinstance(mol, str):
        mol = Molecule.from_smiles(mol, allow_undefined_stereo=True)

    # Parameterize the simulation
    openff_topology = Topology()
    openff_topology.add_molecule(mol)
    openmm_topology = openff_topology.to_openmm()

    try:
        system = ff_engine.create_openmm_system(openff_topology)
    except (
        UnassignedValenceParameterException,
        UnassignedProperTorsionParameterException,
        ValueError,
    ) as error:
        logger.warning(f"Failed to setup the simulation for {mol.to_smiles()} due to a {type(error)}.")
        return
    return openmm_topology, system


def simulate_low_energy_conformers(
    starting_conformation: Union[Quantity, State],
    topology,
    system,
    max_conformers: int = 10,
    n_simulation_steps: int = 1e4,
    energy_threshold: Optional[Quantity] = 1e4 * unit.kilojoules_per_mole,
    temperature: Quantity = 500 * unit.kelvin,
    friction_coeff: Quantity = 1 / unit.picosecond,
    step_size: Quantity = 0.001 * unit.picosecond,
    max_energy_minimization_iters: int = 0,
):
    """
    Find conformers through a MD simulation from a given start confirmation

    Parameters
    ----------
    starting_conformation: Union[Quantity, State]
        The initial conformer used as starting point for the simulation
    topology: openmm.app.Topology
        The topology for the molecule
    system: openmm.System
        The system that is being simulated
    max_conformers: int
        The max number of conformers to return
    n_simulation_steps: int
        The number of simulation steps
    energy_threshold: Optional[Quantity]
        The energy threshold for a conformer to be considered low energy.
    temperature: Quantity
        The temperature for the LangevinMiddleIntegrator
    friction_coeff: Quantity
        The friction_coeff for the LangevinMiddleIntegrator
    step_size: Quantity
        The step_size for the LangevinMiddleIntegrator
    max_energy_minimization_iters: int
        The max number of energy minimization iterations

    Source: https://github.com/openmm/spice-dataset/blob/main/pubchem/createPubchem.py#L44-L56
    """

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction_coeff, step_size)

    # Setup the simulation
    simulation = app.Simulation(topology, system, integrator, 
                                openmm.Platform.getPlatformByName("Reference"))

    if isinstance(starting_conformation, State):
        simulation.context.setState(starting_conformation)
    else:
        simulation.context.setPositions(starting_conformation)

    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.minimizeEnergy(max_energy_minimization_iters)

    # Run the simulation
    states = []
    for i in range(max_conformers):
        simulation.step(n_simulation_steps)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        if energy_threshold is None or state.getPotentialEnergy() < energy_threshold:
            states.append(state)
    return states
