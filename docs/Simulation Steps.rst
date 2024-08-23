Simulation Steps
================

So far, we have shown several ways to configure the elements that make up the simulation. 
However, running such simulations in their current form would cause our calculations to produce no output. 
*Arguably, the most important aspect of conducting a simulation is the extraction and analysis of data from it*. 
To access this data, we utilize **simulation steps**. In essence, a *simulation step* 
is a function that is executed periodically, granting access to the complete simulation information. 
This capability allows us to process and output specific data, either by saving it to a file or displaying it on screen. 
Importantly, simulation steps are not just passive data handlers. They can actively modify the variables of the simulation, 
enabling the implementation of various scenarios and dynamic adjustments during the simulation.

The following example illustrates the addition of three distinct simulation steps. 
The first step, labeled "info", serves to display key simulation data like velocity (steps per second) or remaining time. 
The second, "saveState", records the state of the simulation, particularly the position of the particles, 
into a binary DCD-format file. This format is widely recognized in simulation analysis software, 
facilitating further trajectory analysis. The third step, "thermodynamicMeasurement", 
tracks certain thermodynamic variables, such as energy or temperature. 
Each of these steps includes an "intervalStep" parameter, a standard feature across all *simulation steps* 
that dictates the frequency of which each step is executed:

.. code-block:: yaml

   simulationStep:
     info:
       type: ["UtilsStep", "InfoStep"]
       parameters:
         intervalStep: 10000
     saveState:
       type: ["WriteStep", "WriteStep"]
       parameters:
         intervalStep: 100000
         outputFilePath: "output"
         outputFormat: "dcd"
     thermodynamicMeasurement:
       type: ["ThermodynamicMeasure", "ThermodynamicQuantityMeasure"]
       parameters:
         intervalStep: 10000
         outputFilePath: "thermo.dat"

There are numerous types of simulation steps available, catering to diverse needs. 
For instance, the "[WriteStep, WriteStep]" type, as exemplified above, is used for saving simulation states 
in various formats such as DCD, XYZ, PDB, LAMMPSTRJ, etc.

Using *simulation steps* we are able to not only to measure specific observables, 
like energy distributions or unique experiments such as AFM indentation curves, 
but also to implement a wide array of algorithms. These steps go beyond mere data collection, 
enabling modifications to the simulation state itself. 
This functionality is particularly vital for executing complex algorithms, such as **Thermodynamic Integration** (TI) 
which is a fundamental tool for free energy evaluation. 
We will illustrate how a *simulation step* can be used to provide a protocol for the TI scheme. 
In TI, the objective is to compute the average derivative of a potential energy function U(:math:`\lambda`)
with respect to a parameter :math:`\lambda`, as represented by:

.. math::

   \left\langle\frac{\partial U(\lambda)}{\partial\lambda}\right\rangle_{\lambda}

This calculation extends across various :math:`\lambda` values within the [0,1] interval. 
The aim is to determine the free energy change, ΔF(A → B), by integrating the above quantity over :math:`\lambda`:

.. math::

   \Delta F(A \rightarrow B) = \int_0^1 \left\langle\frac{\partial U(\lambda)}{\partial\lambda}\right\rangle_{\lambda} d\lambda

For successful TI, it is essential to start with selecting suitable :math:`\lambda`-dependent potentials. 
UAMMD-structured provides these potentials for bonds and pairwise interactions. 
The simulation uses the NVT:math:`\lambda` ensemble, defining a global lambda parameter. 
By selecting the right simulation step, ["ThermodynamicMeasure","ThermodynamicIntegration"], the process is set in motion. 
This step requires parameters like the varying values of :math:`\lambda` within the range [1,...,0]. 
The algorithm then computes the derivative in the first equation and records the results for subsequent integration 
as shown in the second equation.

In addition, UAMMD-structured includes *simulation steps* for deriving correlations or essential system properties, 
such as the moment of inertia or the position of the center of mass. 
Calculating pair-defined properties, like the pairwise force matrix (:math:`f_{ij}` for potential-induced force), 
is particularly challenging on GPUs due to their parallel nature. 
Standard GPU functions tend to focus on computing the total force on a particle (:math:`f_i`), 
as each GPU thread is usually linked to a specific particle. UAMMD-structured addresses this by implementing 
a feature that allows the computation of such pairwise matrices, crucial for many coarse-grained parameterization algorithms, 
and the Hessian matrix, relevant for analyses like normal mode analysis. 
Notably, a *simulation step* is available for calculating stress on each particle.
