Input implementation
====================

After thoroughly analyzing the structured input design, JSON (JavaScript Object Notation)
was chosen to implement the input mechanism in UAMMD-Structured. This decision was made
considering several potential file formats that could fit the requirements.

Several formats, including XML, YAML, and TOML, possess the necessary structural capabilities.
Each of these formats supports hierarchical data representation, flexibility in data types, and
ease of use. However, the JSON format offers a balance of human-readability, widespread
support, and adaptability that makes it a preferred choice for UAMMD-Structured.

Some advantages of using JSON in UAMMD-Structured's input mechanism include:

1. **Human-Readable**: JSON is clear and compact, making it easy for users to understand and use.
2. **Hierarchy and Structure**: Naturally, JSON supports nested structures, allowing a seamless
   representation of UAMMD-Structured's input.
3. **Flexibility**: JSON can handle a variety of data types, making it efficient for different
   entries, such as 'Type', 'Parameters Entry', or the 'Matrix Format'.
4. **Widely Supported**: JSON's popularity means there are many tools available for it across
   languages, facilitating integration with UAMMD-Structured.
5. **Extensibility**: JSON's adaptability ensures the system can scale and adapt to future
   advancements in UAMMD-Structured.

While JSON was the chosen format, it's worth noting that UAMMD-Structured's design allows users
to integrate other file formats if desired. While integrating a new format wouldn't be inherently
difficult, it does require meticulous work to ensure complete compatibility and maintain the
integrity of the input design.

By selecting JSON and ensuring the adaptability of UAMMD-Structured, we aim to provide a system
that remains efficient and user-friendly, both now and in the future.

.. code-block:: json

   {
      "Section1":{
        "DataEntry1-1":{
            "type":["class","subclass"],
            "parameters":{
                "param1":"value1",
                "param2":"value2",
                "...":"..."
            }
            "labels":["label1","label2","..."]
            "data":[
                ["label1_data1","label1_data2","..."],
                ["label2_data1","label2_data2","..."],
                ["..."]
            ]
        }
        "DataEntry1-2":{
            "type":["class","subclass"],
            "parameters":{
                "param1":"value1",
                "param2":"value2",
                "...":"..."
            }
            "labels":["label1","label2","..."]
            "data":[
                ["label1_data1","label1_data2","..."],
                ["label2_data1","label2_data2","..."],
                ["..."]
            ]
        }
      },
      "Section2":{
        "DataEntry2-1":{
            "type":["class","subclass"],
            "parameters":{
                "param1":"value1",
                "param2":"value2",
                "...":"..."
            }
            "labels":["label1","label2","..."]
            "data":[
                ["label1_data1","label1_data2","..."],
                ["label2_data1","label2_data2","..."],
                ["..."]
            ]
        }
        "DataEntry2-2":{
            "type":["class","subclass"],
            "parameters":{
                "param1":"value1",
                "param2":"value2",
                "...":"..."
            }
            "labels":["label1","label2","..."]
            "data":[
                ["label1_data1","label1_data2","..."],
                ["label2_data1","label2_data2","..."],
                ["..."]
            ]
        },
        "SectionA":{
          "DataEntryA-1":{
              "type":["class","subclass"],
              "parameters":{
                  "param1":"value1",
                  "param2":"value2",
                  "...":"..."
              }
              "labels":["label1","label2","..."]
              "data":[
                  ["label1_data1","label1_data2","..."],
                  ["label2_data1","label2_data2","..."],
                  ["..."]
              ]
          }
          "DataEntryA-2":{
              "type":["class","subclass"],
              "parameters":{
                  "param1":"value1",
                  "param2":"value2",
                  "...":"..."
              }
              "labels":["label1","label2","..."]
              "data":[
                  ["label1_data1","label1_data2","..."],
                  ["label2_data1","label2_data2","..."],
                  ["..."]
              ]
          }
        }
      }
   }
