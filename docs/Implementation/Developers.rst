Information for developers
==========================


Github actions
**************


The project uses Github actions to run tests and build the documentation. The actions are defined in the `.github/workflows` directory.

During new releases, a conda package is built and uploaded to the `ComplexFluidsUAM` channel on Anaconda Cloud. The package is built using the `devtools/conda-build` directory, which contains the necessary files for building the conda package. The workflow `.github/workflows/conda.yml` is responsible for building the conda package and uploading it to the Anaconda Cloud. Note that it requires the `ANACONDA_TOKEN` secret to be set in the repository settings, which is used to authenticate with Anaconda Cloud.
