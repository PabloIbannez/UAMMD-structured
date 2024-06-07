#!/bin/bash

# Define the source and destination directories
SOURCE_DIR="$(pwd)"
DEST_DIR_MAIN="$HOME/UAMMD"
DEST_DIR="$DEST_DIR_MAIN/extensions"

# Define CUDA target architecture, if not set, it will be compile for all architectures
# If more than one architecture is desired, separate them with a semicolon, e.g. "60;61;70;75;80;86"
export CUDA_ARCH="86"

# Check if DEST_DIR_MAIN exists. If it exists, stop the installation
rm -rf $DEST_DIR_MAIN
#if [ -d "$DEST_DIR_MAIN" ]; then
#                echo "UAMMD is already installed in $DEST_DIR_MAIN"
#                exit 1
#fi

# Create the destination directories if they do not exist
mkdir -p "$DEST_DIR"
mkdir -p "$DEST_DIR_MAIN"

# Copy the files to the destination directories
cp -r "$SOURCE_DIR/structured" "$DEST_DIR/"
cp -r "$SOURCE_DIR/USCM" "$DEST_DIR_MAIN/"
cp "$SOURCE_DIR/preamble.h" "$DEST_DIR_MAIN/extensions"

# Download the UAMMD source code

#mkdir TMP
cd TMP
# git clone https://github.com/PabloIbannez/UAMMD.git

# Copy UAMMD/src to the destination directory
cp -r UAMMD/src "$DEST_DIR_MAIN/"
cd ..
#rm -rf TMP

# Create the folder bin at the destination directory
mkdir -p "$DEST_DIR_MAIN/bin"
mkdir -p "$DEST_DIR_MAIN/launcher"

# Write the launcher code to $DEST_DIR_MAIN/launcher/UAMMDlauncher.cu

echo "
#include \"UAMMDstructured.cuh\"

using namespace uammd::structured;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        uammd::System::log<uammd::System::CRITICAL>(\"No input file provided!\");
        return EXIT_FAILURE;
    }

    std::string inputFilePath = argv[1];
    startSelfStartingSimulation(inputFilePath);

    return EXIT_SUCCESS;
}
" > "$DEST_DIR_MAIN/launcher/UAMMDlauncher.cu"

# Create a symbolic link between DEST_DIR_MAIN/USCM/USCM and DEST_DIR_MAIN/bin/USCM
ln -s "$DEST_DIR_MAIN/USCM/USCM" "$DEST_DIR_MAIN/bin/USCM"

# Check the environment variables UAMMD_PATH and UAMMD_STRUCTURED_PATH exist
# If they do not exist, add them to the .bashrc file
# If they do exist, check they are correct
# If they are not correct, print a warning message
# If they are correct, do nothing
#
# UAMMD_PATH should be DEST_DIR_MAIN
# UAMMD_STRUCTURED_PATH should be DEST_DIR
		
if [ -z "$UAMMD_PATH" ]; then
		echo "##### UAMMD-structured #####" >> ~/.bashrc
		echo "Adding UAMMD_PATH to .bashrc"
		echo "export UAMMD_PATH=$DEST_DIR_MAIN" >> ~/.bashrc
		export UAMMD_PATH=$DEST_DIR_MAIN
else
		if [ "$UAMMD_PATH" != "$DEST_DIR_MAIN" ]; then
				echo "WARNING: UAMMD_PATH is not set to $DEST_DIR_MAIN"
		else
				echo "UAMMD_PATH found and set to $DEST_DIR_MAIN"
		fi
fi

if [ -z "$UAMMD_STRUCTURED_PATH" ]; then
		echo "Adding UAMMD_STRUCTURED_PATH to .bashrc"
		echo "export UAMMD_STRUCTURED_PATH=$DEST_DIR/structured/src" >> ~/.bashrc
		export UAMMD_STRUCTURED_PATH=$DEST_DIR/structured
else
		if [ "$UAMMD_STRUCTURED_PATH" != "$DEST_DIR/structured/src" ]; then
				echo "WARNING: UAMMD_STRUCTURED_PATH is not set to $DEST_DIR"
		else
				echo "UAMMD_STRUCTURED_PATH found and set to $DEST_DIR"
		fi
fi

# Check if DEST_DIR_MAIN/bin/ is in the PATH environment variable
# If it is not, add it to the .bashrc file
# If it is, do nothing

if [[ ":$PATH:" != *":$DEST_DIR_MAIN/bin:"* ]]; then
		echo "Adding $DEST_DIR_MAIN/bin to PATH"
		echo "export PATH=\$PATH:$DEST_DIR_MAIN/bin" >> ~/.bashrc
		export PATH=$PATH:$DEST_DIR_MAIN/bin
		echo "############################" >> ~/.bashrc
else
		echo "$DEST_DIR_MAIN/bin found in PATH"
fi

#python -m pip install --no-cache-dir colorama -U
#python -m pip install --no-cache-dir setuptools_cuda_cpp -U
#
## Install pyUAMMD python lib
#python -m pip install --no-cache-dir pyUAMMD -U

# If folder ./build exists, remove it
if [ -d "./build" ]; then
                rm -rf build
fi

# Create the folder build
mkdir build
cd build
if [ -z "$CUDA_ARCH" ]; then
    cmake ..
else
    cmake -DCUDA_ARCHITECTURES=$CUDA_ARCH ..
fi
make -j VERBOSE=1
make install
cd ..
