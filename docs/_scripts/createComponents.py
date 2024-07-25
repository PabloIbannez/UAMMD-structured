import json
import os

# Load the JSON data
with open("../structured/Components.json", "r") as file:
    components = json.load(file)

subSubType = ["Integrator","SimulationStep","Interactor","DataStructures"]

def generate_rst_content(category,type,subtype,file):
    content = f"""
{subtype}
{'-' * (len(subtype) + 1)}

Category: {category}
Type: {type}
Subtype: {subtype}
File: {file}

[Describe the component]

[Describe the parameters]

[Provide examples]

[Provide links to related components]
"""
    return content.strip()

grouped_categories = {}

createSubSubType = False
for category in components:
    if category in subSubType:
        createSubSubType = True
    else:
        createSubSubType = False
    grouped_categories[category] = {}
    for subcategory in components[category]:
        if createSubSubType:
            grouped_categories[category][subcategory] = {}
        else:
            grouped_categories[category][subcategory] = []
        for component in components[category][subcategory]:
            type, subtype, file = component
            if createSubSubType:
                if type not in grouped_categories[category][subcategory]:
                    grouped_categories[category][subcategory][type] = []
                grouped_categories[category][subcategory][type].append((subtype, file))
            else:
                grouped_categories[category][subcategory].append((type, subtype, file))

for category in grouped_categories:
    if category in subSubType:
        createSubSubType = True
    else:
        createSubSubType = False
    print(f"{category}:")
    # If file category does not exist, create it
    if not os.path.exists(category):
        os.makedirs(category)
    for subcategory in grouped_categories[category]:
        print(f"  {subcategory}:")
        # category/subcategory does not exist, create it
        if not os.path.exists(os.path.join(category, subcategory)):
            os.makedirs(os.path.join(category, subcategory))

        # Create index.rst file
        if createSubSubType:

            with open(os.path.join(category, subcategory, "index.rst"), "w") as findex:
                findex.write(f"{subcategory}\n")
                findex.write("=" * len(subcategory) + "\n\n")
                findex.write(".. toctree::\n")
                findex.write("   :maxdepth: 1\n\n")

                for type in grouped_categories[category][subcategory]:
                    findex.write(f"   {type}/index\n")


            for type in grouped_categories[category][subcategory]:
                print(f"    {type}:")

                # category/subcategory/type does not exist, create it
                if not os.path.exists(os.path.join(category, subcategory, type)):
                    os.makedirs(os.path.join(category, subcategory, type))

                with open(os.path.join(category, subcategory, type, "index.rst"), "w") as findex:
                    findex.write(f"{type}\n")
                    findex.write("=" * len(type) + "\n\n")
                    findex.write(".. toctree::\n")
                    findex.write("   :maxdepth: 1\n\n")

                    for component in grouped_categories[category][subcategory][type]:
                        subtypeName, fileName = component
                        findex.write(f"   {subtypeName}\n")



                for component in grouped_categories[category][subcategory][type]:
                    subtypeName, fileName = component
                    print(f"      {category} {type} {subtypeName} {fileName}")
                    docFilename = os.path.join(category, subcategory, type, subtypeName+".rst")
                    if not os.path.exists(docFilename):
                        with open(os.path.join(docFilename), "w") as file:
                            file.write(generate_rst_content(category,type,subtypeName, fileName))
                    else:
                        print(f"        File {docFilename} already exists")
        else:

            with open(os.path.join(category, subcategory, "index.rst"), "w") as findex:
                findex.write(f"{subcategory}\n")
                findex.write("=" * len(subcategory) + "\n\n")
                findex.write(".. toctree::\n")
                findex.write("   :maxdepth: 1\n\n")

                for type in grouped_categories[category][subcategory]:
                    findex.write(f"   {type[1]}\n")

            for component in grouped_categories[category][subcategory]:
                type, subtypeName, fileName = component
                print(f"      {category} {type} {subtypeName} {fileName}")
                docFilename = os.path.join(category, subcategory, subtypeName+".rst")
                if not os.path.exists(docFilename):
                    with open(os.path.join(docFilename), "w") as file:
                        file.write(generate_rst_content(category,type,subtypeName, fileName))
                else:
                    print(f"        File {docFilename} already exists")
