import os
import sys

inputFile = "index.rst"

#Iterate over the lines in the file
inToc = False
for line in open(inputFile):
    if "toctree::" in line:
        inToc = True
    if inToc and ":" not in line:
        #Check if line contains something different than whitespace
        if line.strip():
            #Remove whitespace at beginning and end of line, and remove newline character
            line = line.strip()
            line = line.rstrip()
            line = line.rstrip("\n")

            #Add .rst extension to line
            file = line + ".rst"

            #Check if file exists
            if os.path.isfile(file):
                #Do nothing
                pass
            else:
                #Create file
                with open(file, "w") as f:
                    f.write(line + "\n")
                    f.write("=" * len(line) + "\n\n")







