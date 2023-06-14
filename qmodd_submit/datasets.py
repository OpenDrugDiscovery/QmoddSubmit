import os
import sys
import time
import json
import h5py
import tqdm
import typer
import fsspec
import logging
import datetime
import subprocess
import numpy as np
import pandas as pd
import datamol as dm

from collections import Counter, OrderedDict
from typing import Optional, Union
from loguru import logger
from openff.units import unit
from openff.toolkit.topology import Molecule
from openff.qcsubmit.factories import BasicDatasetFactory, OptimizationDatasetFactory
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.qcsubmit.datasets import BasicDataset
from qcfractal.interface.models import TaskStatusEnum
import qcfractal.interface as interface
from openff.qcsubmit import workflow_components

from qmodd_submit.conformers import generate_conformers
from qmodd_submit.spec import (single_basis_dft_specification, 
                        geometry_optimization_specification, 
                        semi_empirical_specification)

# disable openff-toolkit warnings
logging.getLogger("openff.toolkit").setLevel(logging.ERROR)


def parse_server_info(server_info_file=None):
    if server_info_file is not None:
        with open(server_info_file, "r") as f:
            server_info = json.load(f)
    else:
        server_info = {"address": "localhost:7777", "verify": False}

    address = server_info["address"]
    verify = server_info.get("verify", False)
    username = server_info.get("username", None)
    password = server_info.get("password", None)

    return dict(
        address=address, 
        verify=verify, 
        username=username, 
        password=password
    )

class CollectionEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (datetime.date, datetime.datetime)):
            return str(obj)
        return json.JSONEncoder.default(self, obj)


def check_computation_status(server_info_file=None):
    server_info = parse_server_info(server_info_file)
    client = interface.FractalClient(**server_info)

    # res = client.get_collection("OptimizationDataset", 
    #                             "qmodd_semiempirical_optimization_geometry_dataset")
    # print(res.)
    completed = client.query_procedures()
    results = list(map(lambda x: x.dict(), completed))
    with open("qmodd_semi_geo_dataset.json", "w") as f:
        json.dump(results, f, indent=2, sort_keys=True, cls=CollectionEncoder)
    # print(res)
    # print(res.)

    # client.


class BaseQCDataset:
    def __init__(
        self,
        server_info_file: Union[str, None] = None,
        loguru_level: str = "INFO",
    ):
        logger.remove()
        logger.add(sys.stderr, level=loguru_level.upper())

        if loguru_level != "DEBUG":
            dm.disable_rdkit_log()

        res = subprocess.run(["git", "config", "user.name"], stdout=subprocess.PIPE)
        self._git_username = res.stdout.strip().decode()

        self._parse_server_info(server_info_file)
        self.connect_to_server()

    def connect_to_server(self):
        try:
            client = interface.FractalClient(**self.server_info)
        except ConnectionRefusedError as error:
            if self.server_info["verify"]:
                raise typer.BadParameter("Connection refused. Try passing the --no-verify flag") from error
            raise error
        self.computation_server = client

    def _parse_server_info(self, server_info_file=None):
        self.server_info = parse_server_info(server_info_file)


class ConformersDataset(BaseQCDataset):
    """
    This script reads a file of molecules, generates a number of diverse, low energy conformations for each one
    using simulations and writes the results to a HDF5 file.
    
    Sources:"""

    NAME = "qmodd_semiempirical_optimization_geometry_dataset"
    TAGLINE="Optimize the genometry before QM modeling"
    DESCRIPTION="ODD Semi-empirical geometry optimization dataset"

    def __init__(
        self,
        input_file: str,
        server_info_file: Union[str, None] = None,
        ff_spec: str = "openff_unconstrained-2.0.0.offxml",
        n_conformers_sampled_per_molecule: int = 200,
        n_jobs: int = 1,
        loguru_level: str = "INFO",
    ):
        super().__init__(server_info_file=server_info_file, loguru_level=loguru_level)
        self.input_file = input_file
        self.ff_spec = ff_spec
        self.ff_engine = ForceField(ff_spec)
        self.n_jobs = n_jobs
        self.n_conformers_sampled_per_molecule = n_conformers_sampled_per_molecule
        
        with fsspec.open(input_file, "r") as fd:
            tmp = pd.read_csv(fd).to_dict(orient="list")
            self.molecules = dict(zip(tmp["ID"], tmp["SMILES"]))
        logger.info(f"Loaded {len(self.molecules)} molecules from {input_file}.")
        

    def _generate_conformers_(self):
        """
        Reads a file of molecules, generates a number of diverse, low energy conformations for each one
        using MD simulations.
        
        Sources:
            -
        """
        n_successful = 0
        
        self.mol_with_conformers = []
        for mol_name in tqdm.tqdm(self.molecules, desc="Molecules"):
            smi = self.molecules[mol_name]
            mol = generate_conformers(smi, self.ff_engine)
            if mol is not None:
                self.mol_with_conformers.append(mol)
                n_successful += 1

        logger.info(f"Successfully embedded {n_successful}/{len(self.molecules)} compounds")

    def submit(self):
        """
        Create a dataset and submit to the server
        This does not yet run any QM computations. It just preprocesses the molecular structures (e.g. deduplication).
        """  

        self._generate_conformers_()

        # Create the factory
        spec = geometry_optimization_specification()
        factory = OptimizationDatasetFactory(
            qc_specifications={spec.qc_spec.spec_name: spec.qc_spec},
        )
        print(factory.dict())

        # Create the dataset
        self._dataset = factory.create_dataset(
            dataset_name=self.NAME,
            molecules=self.mol_with_conformers,
            tagline=self.TAGLINE,
            description=self.DESCRIPTION,
            verbose=True,
        )

        self._dataset.metadata.submitter = self._git_username
        
        responses = self._dataset.submit(self.computation_server, verbose=True)
        logger.info(f"Submitted {len(responses)} tasks to {self.server_info['address']}")


    # def prepare_dft_dataset(self):
    #     """
    #     A dataset factory preprocesses a set of SMILES and contains meta-data about the dataset.
    #     This does not yet run any QM computations. It just preprocesses the molecular structures (e.g. deduplication).
    #     """  

    #     self.opt_dataset

    #     # Create the factory
    #     spec = single_basis_dft_specification()
    #     factory = BasicDatasetFactory(
    #         driver=spec.driver,
    #         qc_specifications={spec.qc_spec.spec_name: spec.qc_spec},
    #         store_wavefunction=spec.qc_spec.store_wavefunction,
    #     )

    #     # Create the dataset
    #     dataset = factory.create_dataset(
    #         dataset_name="qmodd_dft_optimization",
    #         molecules=molecules,
    #         tagline="Optimize the genometry before QM modeling",
    #         description="ODD DFT geometry optimization",
    #         verbose=True,
    #     )

    #     self.dft_dataset = dataset


class SemiEmpiricalDataset(BaseQCDataset):
    def __init__(
        self,
        conformers_dataset_name: str,
        server_info_file: Union[str, None] = None,
        loguru_level: str = "INFO",
    ):
        super().__init__(server_info_file=server_info_file, loguru_level=loguru_level)
        self.conformers_dataset_name = conformers_dataset_name


    def submit(self):
        """
        A dataset factory preprocesses a set of SMILES and contains meta-data about the dataset.
        This does not yet run any QM computations. It just preprocesses the molecular structures (e.g. deduplication).
        """  

        # Create the factory
        spec = semi_empirical_specification()
        factory = BasicDatasetFactory(
            driver=spec.driver,
            qc_specifications={spec.qc_spec.spec_name: spec.qc_spec},
        )

        # Create the dataset
        self._dataset = factory.create_dataset(
            dataset_name=self.NAME,
            molecules=self.mol_with_conformers,
            tagline=self.TAGLINE,
            description=self.DESCRIPTION,
            verbose=True,
        )

        responses = self._dataset.submit(self.computation_server, verbose=True)
        logger.info(f"Submitted {len(responses)} tasks to {self.computation_server['address']}")


class DFTDataset(BaseQCDataset):
    def __init__(
        self,
        conformers_dataset_name: str,
        server_info_file: Union[str, None] = None,
        loguru_level: str = "INFO",
    ):
        super().__init__(server_info_file=server_info_file, loguru_level=loguru_level)
        self.conformers_dataset_name = conformers_dataset_name


    def submit(self):
        """
        A dataset factory preprocesses a set of SMILES and contains meta-data about the dataset.
        This does not yet run any QM computations. It just preprocesses the molecular structures (e.g. deduplication).
        """  

        # Create the factory
        spec = single_basis_dft_specification()
        factory = BasicDatasetFactory(
            driver=spec.driver,
            qc_specifications={spec.qc_spec.spec_name: spec.qc_spec},
            store_wavefunction=spec.qc_spec.store_wavefunction,
        )

        # Create the dataset
        self._dataset = factory.create_dataset(
            dataset_name=self.NAME,
            molecules=self.mol_with_conformers,
            tagline=self.TAGLINE,
            description=self.DESCRIPTION,
            verbose=True,
        )

        responses = self._dataset.submit(self.computation_server, verbose=True)
        logger.info(f"Submitted {len(responses)} tasks to {self.computation_server['address']}")

