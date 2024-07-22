#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "Definitions/Hashes.cuh"

namespace uammd {
namespace structured {
namespace Units {

class UnitsFactory {
public:

    using Creator = std::function<std::shared_ptr<UnitsHandler>(DataEntry&)>;

    static UnitsFactory& getInstance() {
        static UnitsFactory instance;
        return instance;
    }

    void registerUnits(const std::string& unitsType,
                       const std::string& unitsSubType,
                       Creator creator) {
        System::log<System::DEBUG>("[UnitsFactory] Registering Units in factory: %s, %s",
                                     unitsType.c_str(), unitsSubType.c_str());

        std::pair<std::string, std::string> key(unitsType, unitsSubType);
        if (isUnitsRegistered(unitsType, unitsSubType)) {
            System::log<System::CRITICAL>("[UnitsFactory] Units type already registered: %s, %s",
                                         unitsType.c_str(), unitsSubType.c_str());
            throw std::runtime_error("Units type already registered");
        }
        getCreatorsRef()[key] = creator;
    }

    std::shared_ptr<UnitsHandler> createUnits(const std::string& unitsType,
                                              const std::string& unitsSubType,
                                              DataEntry&  data){

        System::log<System::DEBUG>("[UnitsFactory] Creating Units. Type: %s, SubType: %s)",
                                      unitsType.c_str(), unitsSubType.c_str());

        auto& creators = getCreatorsRef();
        std::pair<std::string, std::string> key(unitsType, unitsSubType);
        auto it = creators.find(key);

        if (it != creators.end()) {
            return it->second(data);
        }
        System::log<System::CRITICAL>("[UnitsFactory] Unknown Units type: %s, subType: %s",
                                      unitsType.c_str(), unitsSubType.c_str());
        throw std::runtime_error("Unknown Units type");
    }

    bool isUnitsRegistered(const std::string& unitsType,
                           const std::string& unitsSubType) const {
        std::pair<std::string, std::string> key(unitsType, unitsSubType);
        return getCreatorsRef().find(key) != getCreatorsRef().end();
    }

    const std::unordered_map<std::pair<std::string, std::string>,
                             Creator,
                             PairHash>& getCreators() const {
        return getCreatorsRef();
    }

private:
    UnitsFactory() = default;

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

#define REGISTER_UNITS(type, subType, ...) \
    namespace { \
        struct registerUnits##type##subType { \
            registerUnits##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Units::UnitsFactory::getInstance().registerUnits( \
                        #type,#subType,\
                            [](uammd::structured::DataEntry&  data \
                               ) -> std::shared_ptr<uammd::structured::Units::UnitsHandler> { \
                        return std::make_shared<__VA_ARGS__>(data); \
                    }); \
                } \
            } \
        }; \
        registerUnits##type##subType registerUnits##type##subType##Instance; \
    }
