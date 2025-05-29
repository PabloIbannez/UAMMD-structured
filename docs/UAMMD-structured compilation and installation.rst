UAMMD-structured compilation and installation
=============================================

UAMMD-structured is composed of several components:

- A CUDA/C++ library containing the core components and utilities called **libuammd_structured**.
- An command-line utility called **UAMMDlauncher** that serves as the main interface for running simulations.
- A Python wrapper called **pyUAMMD** that allows users to interact with the C++ library and the command-line utility from Python scripts.



Installing from conda (recommended)
-----------------------------------

You can install UAMMD-structured (library plus CLI utility) and pyUAMMD via conda:

.. code-block:: bash

   conda install -c ComplexFluidsUAM uammd-structured pyuammd-structured

.. hint:: Conda is a package manager that simplifies the installation of software and its dependencies, making it easier to manage different environments and versions. You can read more about conda at https://docs.conda.io/en/latest/. We recommend installing conda via `Miniforge <https://github.com/conda-forge/miniforge>`_. 
   

Installing from source
----------------------

If you need to customize the build or are looking to contribute to the project, you can compile UAMMD-structured from source.

This process involves several steps, including cloning the repository, installing dependencies, and configuring the build system.

Getting the source code
***********************

UAMMD-structured is hosted on GitHub, and you can clone the repository to your local machine using the following command:

.. code-block:: bash

   git clone https://github.com/PabloIbannez/UAMMD-structured.git
   cd UAMMD-structured

This will download the latest development version of UAMMD-structured and place you in the project directory.

Getting the build dependencies
******************************

The `environment.yml` file lists all dependencies, which can be installed with conda, by running the following at the root of the project:

.. code-block:: bash

   conda env create -n uammd-structured

   
.. hint:: conda is a package manager that simplifies the installation of dependencies. If you don't have conda installed, we recommend using the `Miniforge <https://github.com/conda-forge/miniforge>`_ distribution.

.. warning:: Make sure the CUDA version installed in the environment is compatible with your current NVIDIA driver. You can check with `nvidia-smi`.

	     You can choose a different cuda version (e.g. 12.3) when creating the environment by running:

	     .. code:: bash

		       conda env create -n uammd-structured cuda-version==12.3

An environment called `uammd-structured` will be created, which you can activate with:

.. code-block:: bash

   conda activate uammd-structured



Compiling UAMMD-structured
**************************

Once you have the source code and the dependencies installed, you can compile UAMMD-structured. The compilation process uses CMake, a cross-platform build system that generates makefiles or project files for various platforms.

.. code-block:: bash

   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
   make -j4 install

.. hint:: You can adjust the number of cores used for compilation by changing the `-j4` option. For example, `-j10` will use 10 cores for compilation. Beware, though, compiling CUDA code requires vasts amounts of RAM memory.

The above will compile the UAMMD-structured library and the command-line utility. Including the Python wrapper (pyUAMMD) requires adding the following options to the CMake command:

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DBUILD_PYTHON_WRAPPER=ON -DINSTALL_PYTHON_PACKAGE=ON ..


.. hint:: The `CMAKE_INSTALL_PREFIX` variable specifies where the compiled binaries and libraries will be installed. Setting it to `$CONDA_PREFIX` is standard when using conda, and installs everything into the conda environment, making it easier to manage dependencies and versions.

Additional Compilation Options
******************************

UAMMD-structured offers several configuration options that can be set during the CMake configuration step. These options allow you to customize the build to your specific needs.

Floating precision
%%%%%%%%%%%%%%%%%%

Precision can be set to either single (default) or double. While single precision is faster, especially on GPUs, double precision may be necessary for certain simulations requiring higher accuracy. To enable double precision, use:

.. code-block:: bash

   cmake -DUAMMD_PRECISION=DOUBLE ..

CUDA architectures
%%%%%%%%%%%%%%%%%%

You can specify CUDA architectures to compile for, which can significantly reduce compilation time if you know the specific GPU architectures you'll be using. For example:

.. code-block:: bash

   cmake -DCUDA_ARCHITECTURES="70;75" ..
   #Or
   cmake -DCUDA_ARCHITECTURES=OFF ..
   #Or
   cmake -DCUDA_ARCHITECTURES=all-major ..
   

.. warning:: Make sure the CUDA version installed in the environment is compatible with your current NVIDIA driver. You can check with `nvidia-smi`


Not building UAMMDlauncher and the library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The `BUILD_LIBRARY` option (ON by default) allows you to skip building the UAMMD-structured library and the command-line utility (UAMMDlauncher) when set to OFF. This is useful if you only want to build the Python wrapper. To do this, use:

.. code-block:: bash

   cmake -DBUILD_LIBRARY=OFF ..

Unity builds
%%%%%%%%%%%%%%

Unity builds can significantly speed up compilation by combining multiple source files into a single compilation unit. This reduces the overhead of compiling each file separately. To enable unity builds, use:

.. code-block:: bash

   cmake -DCMAKE_UNITY_BUILD=ON -DCMAKE_UNITY_BUILD_BATCH_SIZE=256 ..

The above will compile the entire UAMMD-structured library in a single compilation unit. This will greatly reduce the compilation time and the size of the final binary, but may increase the memory usage during compilation. You can adjust the `CMAKE_UNITY_BUILD_BATCH_SIZE` to control how many files are combined in each unit.

Debian Package
%%%%%%%%%%%%%%

Creating a Debian package allows for easy installation on Debian-based systems. To build the .deb package, configure CMake with:

.. code-block:: bash

   cmake -DBUILD_DEB_PACKAGE=ON ..

Then build the package:

.. code-block:: bash

   make package

This will create a .deb file in the build directory, which can be installed using the `dpkg` command.


Building the conda package
%%%%%%%%%%%%%%%%%%%%%%%%%%

A conda package can be created for easier distribution and installation across different systems. To build the conda package, first create the build environment, from the root of the repository, using:

.. code-block:: bash

   conda env create -f devtools/conda-envs/build_env.yml

Then, build the package using the conda build command:

.. code-block:: bash
		
   conda build -c conda-forge  --output-folder output devtools/conda-build/uammd-structured


Then, the package for pyUAMMD can be built using:

.. code-block:: bash

   conda build -c conda-forge  --output-folder output devtools/conda-build/pyuammd-structured --use-local


.. hint:: The `--use-local` option allows you to use the locally built UAMMD-structured package when building the pyUAMMD package, ensuring that the two are compatible.


The packages will be created in the `output` directory, and can be installed using the `conda install` command:

.. code-block:: bash

   conda create -n structured --use-local uammd-structured pyuammd-structured
   conda activate structured



----

Testing
-------

After compilation, you can test your installation. Navigate to a test directory within the UAMMD-structured project, run the test generation script, move to the results directory, and execute the simulation. Finally, analyze the results using the provided Python script. This process helps ensure that all components are working correctly after compilation.

