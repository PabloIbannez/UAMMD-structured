#pragma once

#include "SimulationStep/SimulationStep.cuh"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace uammd {
namespace structured {
namespace SimulationStep {

class SimulationStepFactory {
public:
    using Creator = std::function<
        std::shared_ptr<SimulationStepBase>(
            std::shared_ptr<ParticleGroup>,
            std::shared_ptr<IntegratorManager>,
            std::shared_ptr<ForceField>,
            DataEntry&,
            std::string)>;

    static SimulationStepFactory& getInstance() {
        static SimulationStepFactory instance;
        return instance;
    }

    void registerSimulationStep(const std::string& simulationStepType,
                                const std::string& simulationStepSubType,
                                Creator creator) {
        System::log<System::DEBUG>("[SimulationStepFactory] Registering SimulationStep in factory: %s, %s",
                                     simulationStepType.c_str(), simulationStepSubType.c_str());

        std::pair<std::string, std::string> key(simulationStepType, simulationStepSubType);
        if (isSimulationStepRegistered(simulationStepType, simulationStepSubType)) {
            System::log<System::CRITICAL>("[SimulationStepFactory] SimulationStep already registered: %s, %s",
                                          simulationStepType.c_str(), simulationStepSubType.c_str());
            throw std::runtime_error("SimulationStep already registered");
        }
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<SimulationStepBase> createSimulationStep(const std::string& simulationStepType,
                                                             const std::string& simulationStepSubType,
                                                             std::shared_ptr<ParticleGroup>     pg,
                                                             std::shared_ptr<IntegratorManager> integrator,
                                                             std::shared_ptr<ForceField>        ff,
                                                             DataEntry&  data,
                                                             std::string name) {
        System::log<System::DEBUG>("[SimulationStepFactory] Creating SimulationStep: %s (type: %s, subType: %s)",
                                     name.c_str(), simulationStepType.c_str(), simulationStepSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(simulationStepType, simulationStepSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(pg, integrator, ff, data, name);
        }
        System::log<System::CRITICAL>("[SimulationStepFactory] Unknown SimulationStep type: %s, subType: %s",
                                      simulationStepType.c_str(), simulationStepSubType.c_str());
        throw std::runtime_error("Unknown SimulationStep type");
    }

    bool isSimulationStepRegistered(const std::string& simulationStepType,
                                    const std::string& simulationStepSubType) {
        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(simulationStepType, simulationStepSubType);
        return creators.find(key) != creators.end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    SimulationStepFactory() = default;

    static std::unordered_map<std::pair<std::string, std::string>,
                              Creator,
                              PairHash>& getCreatorsRef() {
        static std::unordered_map<std::pair<std::string, std::string>,
                                  Creator,
                                  PairHash> creators;
        return creators;
    }
};

}}}

#define REGISTER_SIMULATION_STEP(type, subType, ...) \
    namespace { \
        struct registerSimulationStep##type##subType { \
            registerSimulationStep##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::SimulationStep::SimulationStepFactory::getInstance().registerSimulationStep( \
                        #type, #subType,\
                        [](std::shared_ptr<uammd::ParticleGroup>                 pg, \
                           std::shared_ptr<uammd::structured::IntegratorManager> integrator, \
                           std::shared_ptr<uammd::structured::ForceField>        ff, \
                           uammd::structured::DataEntry&  data, \
                           std::string name ) -> std::shared_ptr<uammd::structured::SimulationStep::SimulationStepBase> { \
                        return std::make_shared<__VA_ARGS__>(pg, integrator, ff, data, name); \
                    }); \
                } \
            } \
        }; \
        registerSimulationStep##type##subType registerSimulationStep##type##subType##Instance; \
    }
