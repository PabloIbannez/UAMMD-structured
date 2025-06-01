Installation
============

pyUAMMD is part of the UAMMD-structured repository, but can be installed as a standalone package. It is recommended to install it using conda, as this will ensure that all dependencies are correctly resolved and installed.

Installing from conda (recommended)
-----------------------------------

To install pyUAMMD using conda, you can use the following command:

.. code-block:: bash

    conda install -c ComplexFluidsUAM pyuammd-structured

This command will install pyUAMMD along with all its dependencies from the `ComplexFluidsUAM` channel. Note that `uammd-structured` is one of these dependencies.


Installing from source
----------------------

pyUAMMD consists of two parts: the Python package that deals with preparing inputs for UAMMD-structured and a C++ wrapper that allows to run simulations directly from Python.

Start by ensuring you have the necessary dependencies for building installed. The UAMMD-structured repository provides an `environment.yml` file that lists all required packages. You can create a conda environment with these dependencies by running the following at the root of the repository:

.. code-block:: bash

    conda env create -f environment.yml
    

Building UAMMD-structured and pyUAMMD
*************************************

pyUAMMD can be compiled and installed alongside UAMMD-structured by calling CMake with the right options:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/desired/install/path -DBUILD_PYTHON_WRAPPER=ON -DINSTALL_PYTHON_PACKAGE=ON -DBUILD_LIBRARY=ON ..
    make -j4 install


.. hint:: The install path in a conda environment is usally set to `$CONDA_PREFIX`, so you can use `-DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX` to install it directly into your conda environment.

The command above, if successful, will build UAMMD-structured, the C++ wrapper and finally install it all to the provided path. The Python wrapper (pyUAMMD) is installed via `pip`, which CMake is configured to do automatically.


Building pyUAMMD when UAMMD-structured is already installed
***********************************************************


You might already have UAMMD-structured installed (either via conda or by building it yourself), and want to build pyUAMMD only. The instructions are similar to the previous case. In this case you can instruct CMake to not build the UAMMD-structured library, in which case CMake will try to find the UAMMD-structured library in the system. You can do this by running:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/desired/install/path -DBUILD_PYTHON_WRAPPER=ON -DINSTALL_PYTHON_PACKAGE=ON -DBUILD_LIBRARY=OFF ..
    make -j4 install



