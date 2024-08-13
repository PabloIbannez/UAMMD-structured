UAMMD-structured installation
=============================

Known issues
------------

- Sometimes BLAS/LAPACK produces different linking errors:

  - Problems in Ubuntu: If you have installed BLAS/LAPACK using the
    `apt install libopenblas-dev` some linking errors of the type:
    `undefined reference to ''` might appear. This is because the
    `libopenblas.so` provided by the package manager lacks some
    symbols. The best solution to solve this problem is to download the
    OpenBLAS source code and compile it yourself and install all the needed
    libraries. You can download the needed files from the https://www.openblas.net/ website.
    This library can be installed by simply doing mkdir build, cd build, cmake -DBUILD_SHARED_LIBS=ON .., make, make install.
    Note that we have to add the -DBUILD_SHARED_LIBS=ON flag to the cmake command to build the shared libraries.
    Probably you will need sudo permissions to install the library in the system.


