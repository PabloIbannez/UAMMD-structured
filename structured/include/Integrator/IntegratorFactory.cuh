#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace uammd {
namespace structured {
namespace Integrator {

class IntegratorFactory {
public:
    using Creator = std::function<
        std::shared_ptr<uammd::Integrator>(
            std::shared_ptr<GlobalData>,
            std::shared_ptr<ParticleGroup>,
            DataEntry&,
            std::string)>;

    static IntegratorFactory& getInstance() {
        static IntegratorFactory instance;
        return instance;
    }

    void registerIntegrator(const std::string& integratorType,
                            const std::string& integratorSubType,
                            Creator creator) {
        System::log<System::MESSAGE>("[IntegratorFactory] Registering Integrator in factory: %s, %s",
                                     integratorType.c_str(), integratorSubType.c_str());

        std::pair<std::string, std::string> key(integratorType, integratorSubType);
        if (isIntegratorRegistered(integratorType, integratorSubType)) {
            System::log<System::CRITICAL>("[IntegratorFactory] Integrator already registered: %s, %s",
                                          integratorType.c_str(), integratorSubType.c_str());
            throw std::runtime_error("Integrator already registered");
        }
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<uammd::Integrator> createIntegrator(const std::string& integratorType,
                                                        const std::string& integratorSubType,
                                                        std::shared_ptr<GlobalData>    gd,
                                                        std::shared_ptr<ParticleGroup> pg,
                                                        DataEntry&  data,
                                                        std::string name) {
        System::log<System::MESSAGE>("[IntegratorFactory] Creating Integrator: %s (type: %s, subType: %s)",
                                     name.c_str(), integratorType.c_str(), integratorSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(integratorType, integratorSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(gd,pg,data,name);
        }
        System::log<System::CRITICAL>("[IntegratorFactory] Unknown Integrator type: %s, subType: %s",
                                      integratorType.c_str(), integratorSubType.c_str());
        throw std::runtime_error("Unknown Integrator type");
    }

    bool isIntegratorRegistered(const std::string& integratorType,
                                const std::string& integratorSubType) const {
        std::pair<std::string, std::string> key(integratorType, integratorSubType);
        return getCreatorsRef().find(key) != getCreatorsRef().end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    IntegratorFactory() = default;

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

#define REGISTER_INTEGRATOR(type, subType, ...) \
    namespace { \
        struct registerIntegrator##type##subType { \
            registerIntegrator##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Integrator::IntegratorFactory::getInstance().registerIntegrator( \
                        #type,#subType,\
                            [](std::shared_ptr<uammd::structured::GlobalData>    gd, \
                               std::shared_ptr<uammd::ParticleGroup> pg, \
                               uammd::structured::DataEntry&  data, \
                               std::string name \
                               ) -> std::shared_ptr<uammd::Integrator> { \
                        return std::make_shared<__VA_ARGS__>(gd,pg,data,name); \
                    }); \
                } \
            } \
        }; \
        registerIntegrator##type##subType registerIntegrator##type##subType##Instance; \
    }
