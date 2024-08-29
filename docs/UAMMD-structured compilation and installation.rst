UAMMD-structured compilation and installation
=============================================

**Introduction**

UAMMD-structured offers multiple compilation options to suit various user needs. This guide will walk you through the process of obtaining, configuring, and compiling UAMMD-structured. We'll cover the general compilation process and then we'll explore the three main compilation options: executable binary, Python wrapper, and Debian package.

Before we begin, ensure you have the necessary prerequisites installed on your system. These include Git for version control, CMake (version 3.8 or higher) for build configuration, CUDA Toolkit (nvcc and CUDA libraries), a C++ compiler compatible with C++14, and Python if you plan to build the Python wrapper.

**Obtaining the Source Code**

To get started, you'll need to clone the UAMMD-structured repository. Open a terminal and run the following commands:

.. code-block:: bash

   git clone https://github.com/PabloIbannez/UAMMD-structured.git
   cd UAMMD-structured

This will download the latest version of UAMMD-structured and place you in the project directory.

**General Compilation Process**

The compilation process for UAMMD-structured follows a standard CMake workflow. First, create a build directory to keep your source tree clean:

.. code-block:: bash

   mkdir build
   cd build

Next, configure the build using CMake:

.. code-block:: bash

   cmake ..

This command generates the necessary makefiles based on your system configuration and the options you've specified. After configuration, you can compile the project:

.. code-block:: bash

   make

The compilation time may vary depending on your system specifications and the options you've chosen.

**Compilation Options**

UAMMD-structured offers several configuration options that can be set during the CMake configuration step. These options allow you to customize the build to your specific needs.

Precision can be set to either single (default) or double. While single precision is faster, especially on GPUs, double precision may be necessary for certain simulations requiring higher accuracy. To enable double precision, use:

.. code-block:: bash

   cmake -DUAMMD_PRECISION=DOUBLE ..

You can specify CUDA architectures to compile for, which can significantly reduce compilation time if you know the specific GPU architectures you'll be using. For example:

.. code-block:: bash

   cmake -DCUDA_ARCHITECTURES="70;75" ..

If you need to use a specific host compiler for CUDA, you can specify it with:

.. code-block:: bash

   cmake -DCUDA_HOST_COMPILER=/path/to/compiler ..

For an interactive configuration interface, you can use `ccmake ..` instead of `cmake ..` (if available on your system). This will open a CMake configuration screen where you can set the different options.

**Specific Compilation Options**

----

**Executable Binary**

The executable binary is the default compilation option. It produces a standalone program that can be run from the command line. To compile the executable binary, simply follow the general compilation process outlined above. The resulting executable, named `UAMMDlauncher`, will be created in the build directory.

If you wish to install the executable, you can run:

.. code-block:: bash

   make install

By default, this installs to `~/bin`. You can change the installation directory during the configuration step:

.. code-block:: bash

   cmake -DINSTALL_DIR=/path/to/install ..

----

**Python Wrapper**

The Python wrapper allows UAMMD-structured to be used within Python scripts, integrating it with other Python-based tools. To build the Python wrapper, configure CMake with the following options:

.. code-block:: bash

   cmake -DINSTALL_PYTHON_PACKAGE=ON -DBUILD_PYTHON_WRAPPER=ON ..

Then compile and install as usual:

.. code-block:: bash

   make
   make install

This will build the Python wrapper and install the pyUAMMD package.

----

**Debian Package**

Creating a Debian package allows for easy installation on Debian-based systems. To build the .deb package, configure CMake with:

.. code-block:: bash

   cmake -DBUILD_DEB_PACKAGE=ON ..

Then build the package:

.. code-block:: bash

   make package

This will create a .deb file in the build directory, which can be installed using the `dpkg` command.

----

Troubleshooting
---------------

One common issue you might encounter is BLAS/LAPACK linking errors, especially on Ubuntu systems. These errors can occur due to incompatibilities between the installed libraries and UAMMD-structured's requirements. If you've installed BLAS/LAPACK using `apt install libopenblas-dev`, you might see `undefined reference to ''` errors. This happens because the `libopenblas.so` provided by the package manager may lack some symbols required by UAMMD-structured.

To resolve this, you can compile OpenBLAS from source. Download the OpenBLAS source code from https://www.openblas.net/, then compile and install it:

.. code-block:: bash

   mkdir build && cd build
   cmake -DBUILD_SHARED_LIBS=ON ..
   make
   sudo make install

The `-DBUILD_SHARED_LIBS=ON` flag is crucial as it ensures the creation of shared libraries. You may need sudo permissions to install the library system-wide.

Testing
-------

After compilation, you can test your installation. Navigate to a test directory within the UAMMD-structured project, run the test generation script, move to the results directory, and execute the simulation. Finally, analyze the results using the provided Python script. This process helps ensure that all components are working correctly after compilation.

