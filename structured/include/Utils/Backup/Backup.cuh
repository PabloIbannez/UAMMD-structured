#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <stdexcept>
#include <cmath>

#include "System/ExtendedSystem.cuh"

namespace uammd{
namespace structured{
namespace Backup{

    bool openFile(std::shared_ptr<ExtendedSystem> sys,
                  std::string outputFilePath,
                  std::ofstream& outputFile,
                  bool binary = false);
}}}

