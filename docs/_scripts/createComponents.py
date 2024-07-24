import json

componentsFile = "../structured/Components.json"

with open(componentsFile, "r") as file:
    components = json.load(file)

for category in components:
    type_subType_file = []
    for subCategory in components[category]:
        for comp in components[category][subCategory]:
            t,s,file = comp
            type_subType_file.append((t,s,file))

    types = set([t for t,s,f in type_subType_file])
    subTypes = {t: set([s for t1,s,f in type_subType_file if t1 == t]) for t in types}
    print(f"Category: {category}")
    for t in types:
        print(f"  Type: {t}")
        print(f"    SubTypes: {subTypes[t]}")



