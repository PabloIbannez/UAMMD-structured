#pragma once
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace uammd {
namespace structured {
namespace Integrator {

// Custom hash function for std::pair<std::string, std::string>
struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

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
        struct registerVCLS##type##subType { \
            registerVCLS##type##subType() { \
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
        }; \
        registerVCLS##type##subType registerVCLS##type##subType##Instance; \
    }
