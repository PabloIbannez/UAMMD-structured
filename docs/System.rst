System
======

To carry out a simulation, we must not only determine the purely physical variables of the system, but often it is also necessary to specify certain technical aspects. These options are indicated in the *System* Section. The *System* Section is the first segment of the input to be processed in a UAMMD-structured simulation. This section establishes some fundamental parameters of the simulation. Currently, within *System*, we can specify two types of Data Entries: one is mandatory ("Information"), and the other is optional ("Backup"). Both entries are categorized under the primary type "Simulation".

.. code-block:: yaml

   system:
     info:
       type: ["Simulation", "Information"]
       parameters:
         name: "example"
         seed: 123456
     backup:
       type: ["Simulation", "Backup"]
       parameters:
         backupFilePath: "backup"
         backupIntervalStep: 10000

The "Information" Data Entry is essential and must always be included in the simulation input. This entry is employed to assign a name to the simulation, which is specified using the "name" parameter. Optionally, the "seed" parameter can be set to fix the seed for random number generation. Fixing a specific seed value is done for ensuring that simulations are reproducible, providing consistent outputs which are needed for validating results or for debugging purposes.

Handling Errors
---------------

The optional "Backup" Data Entry activates backup and restart mechanism of UAMMD-structured. This mechanism is used for maintaining simulation integrity, especially if an error occurs. If activated, the simulation can revert to a previously saved state and attempt to continue from there. While it might initially appear counterintuitive, presuming a repeated failure at the same simulation point, the utility of this feature depends on the nature of the failure. For instance, a simulation might fail due to a particle being subjected to an excessively large force, leading to instability and senseless numerical values. In molecular dynamics, this situation often arises from an inappropriate setting of the time step integration parameter, "dt". Adjusting "dt" can mitigate the occurrence of these extreme events, although it cannot entirely prevent them. This problem intensifies in coarse-grained simulations where a larger "dt" is usually employed and particularly on GPUs, where work is often done with single precision. To address this, UAMMD-structured modifies the random number seed upon a simulation failure, with the intention that this alteration will sufficiently diverge the simulation path to avoid the previously encountered low-probability failure event. The "Backup" entry allows specifying the file path for saving the simulation state ("backupFilePath") and the interval at which these backups occur ("backupIntervalStep").

In the C++/CUDA environment of UAMMD-structured, error handling is sophisticated. Standard CPU exceptions are captured and handled by regular C++ procedures, but CUDA presents a challenge with its "sticky" errors: error flags set by the CUDA API that remain unless the process is restarted. To counter this, UAMMD-structured employs the following strategy: at the initiation of a simulation with restart enabled, a system call creates two processes: a parent and a child. The child runs the simulation, while the parent supervises. If the child process ends due to an error, the parent detects this and can restart the child process using a backup file, thus bypassing the "sticky" error issue in CUDA and continuing the simulation from a previous save point. This procedure is shown in the following flowchart:

.. figure:: /img/backup.png
   :alt: UAMMD-structured error handling

   Flowchart depicting the subprocess management in UAMMD-structured for handling errors in CUDA.

This approach ensures that simulations can be robustly maintained, allowing for continuity even if unexpected computational errors take place.
