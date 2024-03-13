Input format
============

UAMMD-Structured is designed with a focus on a flexible and robust input mechanism. 
This design allows users to adjust the simulation's characteristics without the need for recompilation, 
thus optimizing efficiency and adaptability.

The fundamental components of this input mechanism are known as 'Data Entries.' 
Each 'Data Entry' is designed to initialize a specific entity within the simulation, 
facilitating precise control over various simulation aspects.

The input files are organized into sections, each comprising a set of 'Data Entries' 
grouped by their functional role and purpose in the simulation. 
This organization enhances the ease of managing and identifying simulation parameters.

A key aspect of UAMMD-Structured's input design is its hierarchical nature, 
allowing sections to encompass other sections or multiple 'Data Entries.' 
This structure provides users with substantial flexibility and a coherent setup for simulations.

For clarity, within the 'Integrators' section, users would discover 'Data
Entries' purposed to initialize entities related to integration tasks, such as
the integrator itself or other entities influencing the integration process.

UAMMD-Structured's input approach emphasizes consistency and modularity. 
Consistency aids users in quickly learning and adapting to various simulation setups, 
while modularity ensures a unified initialization structure for each simulation component.

Details of the Data Entry:

- **Type**: This is comprised of two strings: the class and its subclass that
  denote the entity to be initialized. For instance, an entry to initialize a
  Langevin integrator would have the class labeled as "integrator" and the
  subclass as "Langevin."

- **Parameters Entry**: (Optional) This section houses a dictionary detailing
  the entity's required parameters. Using a standard key-value structure, it
  accepts a diverse range of data types from floats and strings to lists.

- **Matrix Format for Data Intensive Entities**: (Optional) Entities with
  extensive data requirements utilize a matrix format. Described using "layers"
  (labels) and "data" (values), this matrix is structured to support data in
  multiple formats, be it numbers, strings, or lists.

Regarding the sections, as mentioned above, they are composed of other sections or data entries.
They are identified by their name (e.g., integrators, topology), although there is also a type entry 
(similar to the one used for data entries) available, although its use is optional and less frequent.

This flexible input mechanism caters to a broad spectrum of simulation configurations, 
effectively meeting diverse requirements without necessitating further modifications.

For a better understanding of the structured input format, refer to the accompanying diagram below:

.. figure:: img/input.svg
   :width: 512
   :align: center

   **Figure 1:** An illustrative diagram showing the structured input format in UAMMD-Structured,
   highlighting the hierarchical organization of 'Data Entry' within different sections.
