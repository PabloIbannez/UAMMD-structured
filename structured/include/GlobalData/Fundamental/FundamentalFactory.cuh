#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "Definitions/Hashes.cuh"
#include "Utils/Plugins/Plugins.cuh"

namespace uammd {
namespace structured {
namespace Fundamental {

class FundamentalFactory {
public:

    using Creator = std::function<std::shared_ptr<FundamentalHandler>(DataEntry&)>;

    static FundamentalFactory& getInstance() {
        static FundamentalFactory instance;
        return instance;
    }

    void registerFundamental(const std::string& fundamentalType,
                             const std::string& fundamentalSubType,
                             Creator creator) {
        System::log<System::DEBUG>("[FundamentalFactory] Registering Fundamental in factory: %s, %s",
                                     fundamentalType.c_str(), fundamentalSubType.c_str());

        std::pair<std::string, std::string> key(fundamentalType, fundamentalSubType);
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<FundamentalHandler> createFundamental(const std::string& fundamentalType,
                                                          const std::string& fundamentalSubType,
                                                          DataEntry&  data){

        System::log<System::DEBUG>("[FundamentalFactory] Creating Fundamental. Type: %s, SubType: %s)",
                                      fundamentalType.c_str(), fundamentalSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(fundamentalType, fundamentalSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(data);
        }
        System::log<System::CRITICAL>("[FundamentalFactory] Unknown Fundamental type: %s, subType: %s",
                                      fundamentalType.c_str(), fundamentalSubType.c_str());
        throw std::runtime_error("Unknown Fundamental type");
    }

    bool isFundamentalRegistered(const std::string& fundamentalType,
                                 const std::string& fundamentalSubType) const {
        std::pair<std::string, std::string> key(fundamentalType, fundamentalSubType);
        return getCreatorsRef().find(key) != getCreatorsRef().end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    FundamentalFactory() = default;

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

#define REGISTER_FUNDAMENTAL(type, subType, ...)                                                                       \
    namespace {                                                                                                        \
    struct registerFundamental##type##subType                                                                          \
    {                                                                                                                  \
        registerFundamental##type##subType()                                                                           \
        {                                                                                                              \
            PLUGIN_REGISTRATION_GUARD("Fundamental" + std::string(#type) + std::string(#subType));                     \
            uammd::structured::Fundamental::FundamentalFactory::getInstance().registerFundamental(                     \
                #type,                                                                                                 \
                #subType,                                                                                              \
                [](uammd::structured::DataEntry& data)                                                                 \
                    -> std::shared_ptr<uammd::structured::Fundamental::FundamentalHandler> {                           \
                    return std::make_shared<__VA_ARGS__>(data);                                                        \
                });                                                                                                    \
        }                                                                                                              \
    };                                                                                                                 \
    registerFundamental##type##subType registerFundamental##type##subType##Instance;                                   \
    }
