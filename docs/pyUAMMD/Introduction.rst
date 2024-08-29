Introduction
============

Manually writing input files can be tedious. To facilitate this task, pyUAMMD has been developed. It is a Python library that adds tools for easily working with UAMMD-structured input files. pyUAMMD includes a single class, ``simulation``, which behaves like a dictionary where different Sections and Data Entries of UAMMD-structured can be specified:

.. code-block:: python

   import pyUAMMD
   sim = pyUAMMD.simulation()

   sim["system"] = {}
   ...

   sim["global"] = {}
   ...

   sim["integrator"] = {}
   ...

   sim["state"] = {}
   sim["state"]["labels"] = ["id", "position", "direction"]
   sim["state"]["data"] = []
   for i in range(N):
       sim["state"]["data"].append([i, [0,0,0], [0,0,0,1.0]])

   sim["topology"] = {}
   sim["topology"]["structure"] = {}
   sim["topology"]["structure"]["labels"] = ["id", "type"]
   sim["topology"]["structure"]["data"] = []
   for i in range(N):
       sim["topology"]["structure"]["data"].append([i, "A"])

   sim["topology"]["forceField"] = {}
   ...
   sim["simulationStep"] = {}
   ...

   sim.write("simulation.json")

The ``simulation`` class has the ``write(...)`` method to write the dictionary content to a JSON file. Writing a dictionary to a JSON file is trivial in Python using the ``json`` library, but this is extremely slow. To address these restrictions, we have developed the **JFIO** (JSON Fast Input Output) library. It uses ``pybind11`` to create a wrapper with C++ functions of ``nlohmann::json``, enabling rapid writing and reading of JSON files where simulation parameters are set.

``simulation`` also includes the ``append`` method, which allows different simulations to be merged into one. This method is central to the implementation of batching. With ``append``, we add the content of another simulation (the added simulation) to the current one (reference simulation). The ``append`` method has two modes:

- "modelId": In this mode, the added simulation is considered part of the same simulation. The result is configured so the number of batches does not increase. Contents of the added simulation (like particle IDs, group definitions, interaction parameters) are adjusted to make the resultant simulation viable. This method is useful for building complex simulations from simpler ones. For example, suppose we have two simulation files, sim1.json and sim2.json, each containing a different protein using the same model. We could create two ``simulation`` objects with pyUAMMD and then add them together. The resultant simulation would contain both proteins, and by adding an interaction potential between the them, it could be used to study their interaction.

- "batchId": Using this mode, the two simulations are joined into one, but each is associated with a different batchId (in addition to configuring the rest of the parameters to make the simulation viable). This is the simplest way to use UAMMD-structured batching and allows for incorporating a large number of simulations into a single large simulation to execute the sub-simulations concurrently.

So far, we have shown how to use pyUAMMD to create simulations with the intention of finalizing this procedure by executing the ``write`` method, which writes a .json file containing the simulation input. This file is then used to run the ``UAMMDlauncher`` binary, initiating the simulation. But pyUAMMD also allows direct execution of the simulation from Python with the ``run`` method:

.. code-block:: python

   import pyUAMMD
   sim = pyUAMMD.simulation()

   sim["system"] = {}
   ...

   sim.run()

This method starts the simulation, yielding a result analogous to using the binary. This approach is intended to make UAMMD-structured more accessible for less advanced users who may not be comfortable using a **command line interface** CLI. This implementation is straightforward with ``pybind11``. Launching the simulation is fully functional and can be executed. It also serves to showcase the potential for integrating UAMMD-structured with Python. Future work aims to extend interoperability with Python, with the first objective being the ability to program simulation steps entirely in Python. This would allow users with basic Python knowledge to access simulation information in real-time to extract or process variables they consider to be relevant.
