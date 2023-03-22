import dataclasses

from qcelemental.models import DriverEnum
from qcelemental.models.results import WavefunctionProtocolEnum
from openff.qcsubmit.common_structures import QCSpec, SCFProperties


@dataclasses.dataclass
class Specification:
    driver: DriverEnum
    qc_spec: QCSpec


def semi_empirical_specification() -> Specification:
    """ Specification for semi-empirical computations """

    # Create the spec
    qc_specifications = QCSpec(
        program="xtb",
        method="GFN2-xTB",
        spec_name="odd_se_default",
        spec_description="ODD SE quantum chemistry specification",
        maxiter=200,
        keywords={
            "wcombine": True,
            "scf_type": "df",
            "accuracy": 1.0,
            "electronic_temperature": 300.0,
            "max_iterations": 50,
            "solvent": "none", # "water",
        }
    )

    return Specification(DriverEnum.gradient, qc_specifications)


def single_basis_dft_specification() -> Specification:
    """ Specification for single basis DFT computations """

    # Create the spec
    qc_specifications = QCSpec(
        program="psi4",
        method="wb97m-d3bj",
        basis="def2-tzvppd",
        spec_name="odd_dft_default",
        spec_description="ODD DFT quantum chemistry specification",
        implicit_solvent=None, # "water",
        maxiter=200,
        store_wavefunction=WavefunctionProtocolEnum.all,
        keywords={
            "wcombine": False,
            "scf_type": "df",
            "scf_properties": [
                SCFProperties.Dipole,
                SCFProperties.Quadrupole,
                SCFProperties.MullikenCharges,
                SCFProperties.LowdinCharges,
                SCFProperties.WibergLowdinIndices,
                SCFProperties.MayerIndices,
                SCFProperties.MBISCharges,
            ],
        },
    )

    return Specification(DriverEnum.gradient, qc_specifications)


def geometry_optimization_specification() -> Specification:
    """ Specification for semi-empirical geometry optimization """
    # Create the spec
    qc_specifications = QCSpec(
        method="GFN2-xTB",
        program="xtb",
        basis=None,
        spec_name="odd_se_geometry",
        spec_description="ODD SE geometry optimization",
        maxiter=200,
        keywords={
            "wcombine": True,
            "scf_type": "df",
            "accuracy": 1.0,
            "electronic_temperature": 300.0,
            "max_iterations": 50,
            "solvent": "none", # "water",
        }
    )

    return Specification(DriverEnum.gradient, qc_specifications)
