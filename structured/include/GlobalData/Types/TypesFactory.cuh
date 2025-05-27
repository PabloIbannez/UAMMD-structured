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
namespace Types {

class TypesFactory {
public:

    using Creator = std::function<std::shared_ptr<TypesHandler>(DataEntry&)>;

    static TypesFactory& getInstance() {
        static TypesFactory instance;
        return instance;
    }

    void registerTypes(const std::string& typesType,
                       const std::string& typesSubType,
                       Creator creator) {
        System::log<System::DEBUG>("[TypesFactory] Registering Types in factory: %s, %s",
                                     typesType.c_str(), typesSubType.c_str());

        std::pair<std::string, std::string> key(typesType, typesSubType);
        if (isTypesRegistered(typesType, typesSubType)) {
            System::log<System::CRITICAL>("[TypesFactory] Types already registered: %s, %s",
                                         typesType.c_str(), typesSubType.c_str());
            throw std::runtime_error("Types already registered");
        }
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<TypesHandler> createTypes(const std::string& typesType,
                                              const std::string& typesSubType,
                                              DataEntry&  data){

        System::log<System::DEBUG>("[TypesFactory] Creating Types. Type: %s, SubType: %s)",
                                      typesType.c_str(), typesSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(typesType, typesSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(data);
        }
        System::log<System::CRITICAL>("[TypesFactory] Unknown Types type: %s, subType: %s",
                                      typesType.c_str(), typesSubType.c_str());
        throw std::runtime_error("Unknown Types type");
    }

    bool isTypesRegistered(const std::string& typesType,
                           const std::string& typesSubType) const {
        std::pair<std::string, std::string> key(typesType, typesSubType);
        return getCreatorsRef().find(key) != getCreatorsRef().end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    TypesFactory() = default;

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

#include <source_location>

#define REGISTER_TYPES(type, subType, ...)                                                                             \
    namespace {                                                                                                        \
    struct registerTypes##type##subType                                                                                \
    {                                                                                                                  \
        registerTypes##type##subType()                                                                                 \
        {                                                                                                              \
            uammd::structured::PluginUtils::registrationGuard("TypesFactory" + std::string(#type) +                    \
                                                              std::string(#subType) +                                  \
                                                              std::source_location::current().file_name());            \
            uammd::structured::Types::TypesFactory::getInstance().registerTypes(                                       \
                #type,                                                                                                 \
                #subType,                                                                                              \
                [](uammd::structured::DataEntry& data) -> std::shared_ptr<uammd::structured::Types::TypesHandler> {    \
                    return std::make_shared<__VA_ARGS__>(data);                                                        \
                });                                                                                                    \
        }                                                                                                              \
    };                                                                                                                 \
    registerTypes##type##subType registerTypes##type##subType##Instance;                                               \
    }
