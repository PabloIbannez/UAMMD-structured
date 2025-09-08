set -x
SRC_DIR=$(pwd)
TMPDIR=$(mktemp -d)
cd $TMPDIR

if [ -z ${CUDAARCHS+x} ]; then
  CUDAARCHS=$(nvcc --list-gpu-code | tr ' ' '\n' \
      | grep -E '^sm_[0-9]+$' \
      | sed 's/sm_//;s/$/-real/' \
      | paste -sd';' -)
  # Add the virtual architecture of the latest
  LATEST_ARCH=$(echo $CUDAARCHS | tr ';' '\n' | sort -V | tail -n1 | sed 's/-real//')
  CUDAARCHS="${CUDAARCHS};${LATEST_ARCH}"
fi
cmake $SRC_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release -DINSTALL_PYTHON_PACKAGE=ON -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_CUDA_ARCHITECTURES=${CUDAARCHS} -DPython_EXECUTABLE=${PYTHON} -DUAMMD_PRECISION=$UAMMD_PRECISION -DCMAKE_UNITY_BUILD=ON -DCMAKE_UNITY_BUILD_BATCH_SIZE=256
make -j4
make install
rm -rf $TMPDIR
