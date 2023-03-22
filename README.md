# QMODD

## Inspiration for this repository
- [openmm/spice-dataset](https://github.com/openmm/spice-dataset): The Spice dataset that is the small scale of ODD.
  - [/pubchem](https://github.com/openmm/spice-dataset/tree/main/pubchem): Describes the procedure for how the PubChem subset was created.
  - [/downloader](https://github.com/openmm/spice-dataset/tree/main/downloader): Describes how to load the original dataset from [QCArchive](https://qcarchive.molssi.org/).
  - [/issues/54](https://github.com/openmm/spice-dataset/issues/54): Relevant issue on reproducing the SPICE dataset.
- [openforcefield/qca-dataset-submission](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-11-09-QMDataset-pubchem-set6-single-points): Repo that contains the submission to QCArchive.

## Installation

| :information_source: INFO                                                                                                                                                                                                                                                                                          |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Installing the environment will take time. If you check the `env.yml` file, you'll see it includes some very specific versions.  There's several issues on Github on improving the overall process. See e.g. [here](https://github.com/psi4/psi4/issues/2300) or [here](https://github.com/psi4/psi4/issues/2621). |

```bash
mamba create -n qmodd -f env.yml
mamba activate qmodd
pip install -e /path/to/qmodd
```

## Usage
This project makes four commands available through the `qmodd_submit` CLI: 
- `geo` to submit computing and saving conformers for a list of molecules from a `.txt` file.
- `semi` to submit computing semi-empirical properties for a conformer dataset.
- `dft` to submit computing dft properties for a conformer dataset.
- `status` to query the status of computations on  the server.

Get started by running: 
```shell
qmodd --help
```

### Start a QCEngine server locally for testing
To use the `submit` argument, you're required to have a QCEngine. There's several possibilities documented [here](http://docs.qcarchive.molssi.org/projects/qcfractal/en/latest/setup_quickstart.html#setup-overview).

For debugging, consider using the local setup (as documented [here](http://docs.qcarchive.molssi.org/projects/qcfractal/en/latest/setup_quickstart.html#single-workstation)): 
```shell
qcfractal-server init # you have to run this only once to create the local database
nohup qcfractal-server start --local-manager 16 &  # 16 is the number of workers to use
qmodd_submit geo /path/to/dataset.csv
```

If you're using a Notebook, you can find a way to debug the client-server setup [here](http://docs.qcarchive.molssi.org/projects/QCFractal/en/stable/quickstart.html?highlight=FractalClient#Importing-QCFractal).