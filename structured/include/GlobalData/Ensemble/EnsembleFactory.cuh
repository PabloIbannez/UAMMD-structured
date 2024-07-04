#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "Definitions/Hashes.cuh"

namespace uammd {
namespace structured {
namespace Ensemble {

class EnsembleFactory {
public:

    using Creator = std::function<std::shared_ptr<EnsembleHandler>(DataEntry&)>;

    static EnsembleFactory& getInstance() {
        static EnsembleFactory instance;
        return instance;
    }

    void registerEnsemble(const std::string& ensembleType,
                          const std::string& ensembleSubType,
                          Creator creator) {
        System::log<System::MESSAGE>("[EnsembleFactory] Registering Ensemble in factory: %s, %s",
                                     ensembleType.c_str(), ensembleSubType.c_str());

        std::pair<std::string, std::string> key(ensembleType, ensembleSubType);
        if (isEnsembleRegistered(ensembleType, ensembleSubType)) {
            System::log<System::CRITICAL>("[EnsembleFactory] Ensemble already registered: %s, %s",
                                         ensembleType.c_str(), ensembleSubType.c_str());
            throw std::runtime_error("Ensemble already registered");
        }
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<EnsembleHandler> createEnsemble(const std::string& ensembleType,
                                                    const std::string& ensembleSubType,
                                                    DataEntry&  data){

        System::log<System::MESSAGE>("[EnsembleFactory] Creating Ensemble. Type: %s, SubType: %s)",
                                      ensembleType.c_str(), ensembleSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(ensembleType, ensembleSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(data);
        }
        System::log<System::CRITICAL>("[EnsembleFactory] Unknown Ensemble type: %s, subType: %s",
                                      ensembleType.c_str(), ensembleSubType.c_str());
        throw std::runtime_error("Unknown Ensemble type");
    }

    bool isEnsembleRegistered(const std::string& ensembleType,
                              const std::string& ensembleSubType) const {
        std::pair<std::string, std::string> key(ensembleType, ensembleSubType);
        return getCreatorsRef().find(key) != getCreatorsRef().end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    EnsembleFactory() = default;

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

#define REGISTER_ENSEMBLE(type, subType, ...) \
    namespace { \
        struct registerEnsemble##type##subType { \
            registerEnsemble##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Ensemble::EnsembleFactory::getInstance().registerEnsemble( \
                        #type,#subType,\
                            [](uammd::structured::DataEntry&  data \
                               ) -> std::shared_ptr<uammd::structured::Ensemble::EnsembleHandler> { \
                        return std::make_shared<__VA_ARGS__>(data); \
                    }); \
                } \
            } \
        }; \
        registerEnsemble##type##subType registerEnsemble##type##subType##Instance; \
    }
