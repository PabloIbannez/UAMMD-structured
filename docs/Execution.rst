Execution
=========

Initiating a simulation with UAMMD-structured begins with preparing the input, the details of which are covered in following sections. There are primarily two methods to conduct simulations:

1. **Command Line Execution**:
   UAMMD-structured can be compiled into a single executable and launched via the command line:

   .. code-block:: bash

      UAMMDlauncher simulation.json

   Here, ``UAMMDlauncher`` represents the compiled executable, and "simulation.json" is the designated input file.

2. **Python Integration**:
   UAMMD-structured is equipped for Python integration. In this approach, the simulation input is formatted similarly to a Python dictionary. The simulation is then executed by invoking the ``run()`` method:

   .. code-block:: python

      import pyUAMMD
      sim = pyUAMMD.simulation()

      sim[...] = {...}  # Set up simulation parameters
      ...

      sim.run()

Both methods are functionally equivalent, giving users the flexibility to choose the one that best suits their workflow. The key here lies in setting up the input correctly to realize the intended simulation, contrasting with UAMMD where new code needs to be written and compiled for each new simulation type.

