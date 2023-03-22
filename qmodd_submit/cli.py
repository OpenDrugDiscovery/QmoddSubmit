import numpy as np
import typer
import fsspec
import logging
import datamol as dm
from qmodd_submit.datasets import ConformersDataset, SemiEmpiricalDataset, DFTDataset, check_computation_status

# disable openff-toolkit warnings
logging.getLogger("openff.toolkit").setLevel(logging.ERROR)

# CLI
app = typer.Typer(
    help="ODD cli for QM data generation", add_completion=False, pretty_exceptions_enable=False
)

@app.command("geo", help="Compute the geometry optimization for all molecules in the given input_file")
def geo(
    input_file: str = typer.Argument(..., help="Path to the input file"),
    server_info_file: str = typer.Option(default=None, help="Path to the server info file"),
):
    cd = ConformersDataset(input_file=input_file, server_info_file=server_info_file)
    cd.submit()


@app.command("semi", help="Compute semi-empirical properties for all molecules in the optimization dataset")
def semi(
    optimization_dataset: str = typer.Argument(..., help="The dataset name for starting geometries"),
    server_info_file: str = typer.Option(default=None, help="Path to the server info file"),
):
    cd = SemiEmpiricalDataset(input_file=input_file, server_info_file=server_info_file)
    cd.submit()


@app.command("dft", help="Compute DFT properties for the molecules in the optimization dataset")
def dft(
    optimization_dataset: str = typer.Argument(..., help="The dataset name for starting geometries"),
    server_info_file: str = typer.Option(default=None, help="Path to the server info file"),
):

    cd = DFTDataset(input_file=input_file, server_info_file=server_info_file)
    cd.submit()


@app.command("status", help="Check the computation status on the server")
def status(
    server_info_file: str = typer.Option(default=None, help="Path to the server info file"),
):
    """
    Submit a geometry optimization for all molecules in a hdf5 file.
    """
    check_computation_status(server_info_file=server_info_file)





if __name__ == "__main__":
    app()
    
