#include <csignal>
#include <atomic>

#include "../../ThirdParty/pybind11_json.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Simulation/Simulation.cuh"

using namespace uammd::structured;

std::atomic<bool> g_stop_simulation(false);
std::shared_ptr<ExtendedSystem> g_sys;

void sigint_handler(int signal) {
    g_stop_simulation.store(true);
    if (g_sys) {
        g_sys->setState(ExtendedSystem::SIMULATION_STATE::STOPPED);
    }
}

int UAMMDlauncher(const nlohmann::json& json) {

    uammd::System::log<uammd::System::WARNING>("UAMMD-structured Python wrapper"
            " is not compatible with UAMMD-structured self restarting mechanism");

    // Set up the SIGINT handler
    std::signal(SIGINT, sigint_handler);

    auto input = std::make_shared<Input::Input>(json);
    g_sys = std::make_shared<ExtendedSystem>(input);

    std::shared_ptr<Simulation> sim = std::make_shared<Simulation>(g_sys);

    int exitCode = 0;
    try {
        exitCode = sim->run();

        if (g_stop_simulation.load()) {
            uammd::System::log<uammd::System::MESSAGE>("Simulation stopped by user (SIGINT)");
            uammd::System::log<uammd::System::MESSAGE>("Finishing simulation...");
        }
    }
    catch (const std::exception& e) {
        uammd::System::log<uammd::System::ERROR>("Exception caught: %s", e.what());
        exitCode = -1;
    }

    g_sys->finish();
    g_sys.reset();

    return exitCode;
}

// ----------------
// Python interface
// ----------------
namespace py = pybind11;
PYBIND11_MODULE(pyUAMMDlauncher, m)
{
    m.doc() = "pyUAMMDlauncher: Python wrapper for UAMMDlauncher";
    m.def("UAMMDlauncher", &UAMMDlauncher, "Run UAMMDlauncher simulation from JSON input");
}
