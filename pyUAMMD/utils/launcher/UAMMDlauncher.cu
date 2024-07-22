#include "../../ThirdParty/pybind11_json.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "UAMMDstructured.cuh"

using namespace uammd::structured;

int UAMMDlauncher(const nlohmann::json& json) {
    startSelfStartingSimulationFromInput(json);
    return EXIT_SUCCESS;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

PYBIND11_MODULE(UAMMDlauncher,m)
{
  m.doc() = "UAMMDlauncher";

  m.def("UAMMDlauncher", &UAMMDlauncher, "UAMMDlauncher");
}
