Input
=====

To implement the UAMMD-structured input, the primary task involves
creating a series of functions that parse the JSON/YAML file and
establish a C++ interface with its contents. For efficient parsing of
the JSON file, we utilize the popular nlohmann::json library, which
allows for easy navigation through such files. For YAML files, they are
first internally converted to JSON, and then processed by
nlohmann::json, treating them as standard JSON files. Concurrently, we
create a C++ structure that represents the concept of a Data Entry,
equipping it with a series of methods to access its content. This
structure is named ‘DataEntry‘, and many functions and classes of
UAMMD-structured take one of these structures as parameters.

In addition to ‘DataEntry‘, UAMMD-structured also features reduced
versions. These can be used if one wishes to consult certain parts of
the ‘DataEntry‘ without the need to process it in its entirety (which
could lead to inefficient resource use). To access only the type, we can
use ‘TypeEntry‘, and for the parameters, ‘ParameterEntry‘.

The input as a whole is managed by the ‘InputJSON‘ class. When this
class is created, the file is processed by nlohmann::json, and the
result is stored within the ‘InputJSON‘ class. Therefore, to access the
different data entries, we can do so through this class. For example,
let us imagine we want to access the information located in the entry
"DataEntryA-2":

.. code-block:: yaml

    Section2 :
        DataEntry2 -1:
            # ...
        DataEntry2 -2:
            # ...
        SectionA :
            DataEntryA -1:
                # ...
            DataEntryA -2:
                type : [ " type " , " subtype " ]
                parameters :
                    param1 : " value1 "
                    param2 : " value2 "
                     ...   : " ... "
                labels : [ " label1 " , " label2 " , " ... " ]
                data :
                - [ " label1_data1 " , " label1_data2 " , " ... " ]
                - [ " label2_data1 " , " label2_data2 " , " ... " ]
                - [ " ... " ]
        # ...

.. line-block::
    To achieve this, we would request the ‘DataEntry‘ from the ‘input‘ object (of type ‘InputJSON‘) in the path: "Section2/SectionA/DataEntryA-2", which is done as follows:
    ``input->getDataEntry("Section2","SectionA","DataEntryA-2")``. This method would return a ‘DataEntry‘ object with processed information and a series of methods to access it, which return structures with easily manageable data types. The available methods of ‘DataEntry‘ are:

-  **getType()** ``std::string getType()`` and ``std::string getSubType()``: These
   return a string with the type and subtype labels of the entry, respectively.

-  **getParameter()** ``template<typename T> T getParameter(std::string parameterName)``
   and
   ``template<typename T> T getParameter(std::string parameterName, T defaultValue)``:
   To obtain the value of parameters, the second version returns the
   value ‘defaultValue‘ if such parameter is not found.

-  **getData()** ``template<typename T> std::vector<T> getData(std::string label)``
   and ``int getDataSize()``: These methods return the data
   corresponding to the selected label as a vector, facilitating its
   management. With ``getDataSize()``, we can know the size of the data
   before constructing the vector, which can be useful for memory
   management.

There are other methods for accessing data that return different data
structures, which may be more useful depending on the intended use.

Through ‘DataEntry‘, we can also update data using a set of methods.
These allow for updating the input and subsequently exporting it to a
file, which will be key for the backup process.

Surprisingly, one of the most common issues among users is entering
parameters that are ignored due to minor errors in their writing. To
solve these problems, the input keeps a record of all entries that have
been consulted. Before starting the simulation, it checks that all these
entries have been used; otherwise, an error is issued indicating the
entries that have not been employed.

The use of ‘DataEntry‘ is common in UAMMD-structured. In particular,
many classes or functions take a ‘DataEntry‘ as parameters for
initialization or execution. Generally, they tend to take instances of
other objects and subsequently a ‘DataEntry‘, which they process to
obtain the necessary parameters:

.. code-block:: cpp

    void foo(std::shared_ptr<System> sys,
             std::shared_ptr<ParticleData> pd,
             ...,
             DataEntry& data,
             ...){
        float param1 = data.getParameter<float>("param1");
        ...
        std::vector<int> vec1 = data.getData<int>("label1");
        ...
    }


