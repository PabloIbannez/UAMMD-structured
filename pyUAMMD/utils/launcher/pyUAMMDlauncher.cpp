#include "../../ThirdParty/pybind11_json.hpp"
#include "UAMMDlauncher.hpp"

namespace py = pybind11;
PYBIND11_MODULE(pyUAMMDlauncher, m)
{
    m.doc() = "pyUAMMDlauncher: Python wrapper for UAMMDlauncher";
    m.def("UAMMDlauncher", &UAMMDlauncher, "Run UAMMDlauncher simulation from JSON input");
}
