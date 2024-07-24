#include "Simulation/Simulation.cuh"

using namespace uammd::structured;

int main(int argc, char *argv[]) {


    if (argc < 2) {
        uammd::System::log<uammd::System::CRITICAL>("No input file provided!");
        return EXIT_FAILURE;
    }

    std::string inputFilePath = argv[1];
    startSelfStartingSimulation(inputFilePath);

    return EXIT_SUCCESS;
}
