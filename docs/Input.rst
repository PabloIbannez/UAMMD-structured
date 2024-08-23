Input
=====

UAMMD-structured is designed around an input system that is flexible enough to adjust the characteristics of the simulation without the need to recompile the entire program.

Data Entry
----------

The central element of the input system is the **Data Entry**. A Data Entry is a data structure designed to initialize a specific entity of the simulation. The different Data Entries are organized hierarchically. They can be grouped in Sections. These Sections can contain several Data Entries or also other Sections. Both Sections and Data Entries are identified by a **name**, which has to be *unique*. For clarity, in the *Integrators* Section, users will find different Data Entries, each uniquely identified by names like "langevin", "brownian" ... These entries are designed to initialize entities associated with the integration task. This includes not only the integrator itself but also other entities that have an impact on the integration process.

The UAMMD-structured input system is designed with two key principles in mind: consistency and modularity. Consistency ensures ease of learning, while modularity facilitates easy configuration and extension. Each Data Entry contains a series of fields to store information in different formats:

1. **Type**: Defined by two distinct strings, the class and its subclass, this field identifies the specific entity to be initialized. For instance, a Langevin integrator initialization would label its class as "integrator" and the subclass as "Langevin".

2. **Parameters Entry**: This field contains a dictionary that permits to specify some parameters of the Data Entry in detail, using a standard key-value format. It is versatile enough to support a wide variety of data types, such as floats, strings, lists, and more. As an example, consider a Data Entry for a harmonic spring potential, where a typical parameter like "K", representing the spring constant, might be specified.

3. **Matrix Format for Data-Intensive Entities**: Employed for entities requiring extensive data, this format utilizes a matrix structure. Outlined through "labels" (keys) and "data" (values), it is designed to be flexible, supporting varied data formats such as numerical values, strings, or lists. For instance, in setting parameters for a Lennard-Jones interaction, pairs of particle types are related with their respective ε and σ parameters. Assuming that there are two types of particles, "A" and "B", the parameters of the interaction would be specified in the following format:

   +--------+--------+--------+----------+---------+
   | Labels | Name i | Name j | ε        | σ       |
   +========+========+========+==========+=========+
   |        | A      | A      | ε\_AA    | σ\_AA   |
   +        +--------+--------+----------+---------+
   |  Data  | A      | B      | ε\_AB    | σ\_AB   |
   +        +--------+--------+----------+---------+
   |        | B      | B      | ε\_BB    | σ\_BB   |
   +--------+--------+--------+----------+---------+

   Then to access a particular parameter we must specify its corresponding row and label.

Sections, on the other hand, can have a Type field and contain a set of Data Entries. The presence of parameters, data and labels fields in a Data Entry is optional, depending on the needs of the type of entity to be initialized. A Data Entry is still considered valid whether these components are included or not. By contrast, the Type field is essential as it defines the nature of the entry. While the name serves an identification purpose and can vary widely, the Type field is crucial for the program to understand what the Data Entry represents.

.. figure:: /img/input2.png
    :alt: UAMMD-structured input format

    Diagram showing the input format in UAMMD-structured, highlighting the hierarchical organization of Data Entry within different Sections.

There are notable exceptions to this structure. The six principal elements—**system**, **global**, **state**, **integrators**, **topology**, and **simulation steps**, along with the two main elements within **topology** (**structure** and **forcefield**)—are identified uniquely by their names. These elements are so central to the functioning of UAMMD-structured that they do not require a Type field; their names alone suffice for identification purposes.

Input Formats
-------------

UAMMD-structured employs both JSON (JavaScript Object Notation) and YAML (YAML Ain't Markup Language) for its input mechanism. These formats were selected due to their suitability for the aforementioned requirements. JSON and YAML both offer hierarchical data structuring, flexibility in handling diverse data types, and user-friendliness. They provide a balance between human-readability and adaptability. Key benefits of utilizing JSON and YAML include:

- **Human-Readable and Write-Friendly**: Both formats are clear and accessible, simplifying user interaction and understanding.
- **Hierarchical Structure**: They support nested structures, aligning well with the hierarchical input of UAMMD-structured.
- **Flexibility in Data Types**: Efficient in managing various data types, crucial for different entries like "Type", "parameters", and matrix like format data.
- **Broad Support**: Due to their popularity, numerous tools across languages are available.

Although JSON and YAML were chosen, the design of UAMMD-structured permits integration with other formats, requiring careful implementation to maintain compatibility.

In this documentation, YAML is mainly used for its better readability. Yet, it is important to note that for the specific data structures used in UAMMD-structured input, both YAML and JSON formats are equivalent, allowing for easy conversion between the two.

Example
-------

To illustrate the implementation of the input system of UAMMD-structured using YAML, the following examples demonstrate how the previously described input format can be represented using YAML. These examples are provided to give a clear idea of how the structured input format, as discussed earlier, is translated into practical YAML templates.

The first example shows a basic YAML structure for a simple Section of the input file:

.. code-block:: yaml

   Section1:
     DataEntry1-1:
       type: ["class", "subclass"]
       parameters:
         param1: "value1"
         param2: "value2"
         "...": "..."
       labels: ["label1", "label2", "..."]
       data:
         - ["label1_data1", "label1_data2", "..."]
         - ["label2_data1", "label2_data2", "..."]
         - ["..."]
     DataEntry1-2:
       type: ["class", "subclass"]
       parameters:
         param1: "value1"
         param2: "value2"
         "...": "..."
       labels: ["label1", "label2", "..."]
       data:
         - ["label1_data1", "label1_data2", "..."]
         - ["label2_data1", "label2_data2", "..."]
         - ["..."]

The second example demonstrates how nested Sections are represented in YAML:

.. code-block:: yaml

   Section2:
     DataEntry2-1:
       type: ["class", "subclass"]
       parameters:
         param1: "value1"
         param2: "value2"
         "...": "..."
       labels: ["label1", "label2", "..."]
       data:
         - ["label1_data1", "label1_data2", "..."]
         - ["label2_data1", "label2_data2", "..."]
         - ["..."]
     # ...
     SectionA:
       DataEntryA-1:
         type: ["class", "subclass"]
         parameters:
           param1: "value1"
           param2: "value2"
           "...": "..."
         labels: ["label1", "label2", "..."]
         data:
           - ["label1_data1", "label1_data2", "..."]
           - ["label2_data1", "label2_data2", "..."]
           - ["..."]
     # ...

These examples provide a clear illustration of how the structured input format of UAMMD-structured is implemented using YAML, demonstrating both simple and nested structures.
