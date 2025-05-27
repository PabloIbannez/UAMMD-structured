#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "Utils/Plugins/Plugins.cuh"

namespace uammd {
namespace structured {
namespace Interactor {

    class PatchesFactory {

        public:

            using InteractorPtr    = std::shared_ptr<uammd::Interactor>;
            using GlobalDataPtr    = std::shared_ptr<GlobalData>;
            using ParticleGroupPtr = std::shared_ptr<ParticleGroup>;
            using VCLSPtr          = std::shared_ptr<VerletConditionalListSetBase>;

            using PatchesCreator          = std::function<InteractorPtr(GlobalDataPtr, ParticleGroupPtr,
                                                                        GlobalDataPtr, ParticleGroupPtr,
                                                                        DataEntry&, std::string)>;

            using NonBondedPatchesCreator = std::function<InteractorPtr(GlobalDataPtr, ParticleGroupPtr,
                                                                        GlobalDataPtr, ParticleGroupPtr,
                                                                        VCLSPtr,
                                                                        DataEntry&, std::string)>;

            static PatchesFactory& getInstance() {
                static PatchesFactory instance;
                return instance;
            }

            template<typename CreatorType>
            void registerInteractor(const std::string& integratorType, const std::string& integratorSubType, CreatorType creator) {
                System::log<System::DEBUG>("[PatchesFactory] Registering Interactor in factory: %s, %s",
                                             integratorType.c_str(), integratorSubType.c_str());

                std::pair<std::string, std::string> key(integratorType, integratorSubType);
                getCreatorsRef<CreatorType>()[key] = creator;
            }

            template<typename CreatorType, typename... Args>
            InteractorPtr createInteractor(const std::string& integratorType, const std::string& integratorSubType, Args&&... args) {
                System::log<System::DEBUG>("[PatchesFactory] Creating Interactor: %s (type: %s, subType: %s)",
                                             std::get<sizeof...(Args) - 1>(std::forward_as_tuple(args...)).c_str(), // name is always the last argument
                                             integratorType.c_str(), integratorSubType.c_str());

                std::pair<std::string, std::string> key(integratorType, integratorSubType);
                auto& creators = getCreatorsRef<CreatorType>();
                auto it = creators.find(key);

                if (it != creators.end()) {
                    return it->second(std::forward<Args>(args)...);
                }

                System::log<System::CRITICAL>("[PatchesFactory] Unknown Interactor type: %s, subType: %s",
                                              integratorType.c_str(), integratorSubType.c_str());
                throw std::runtime_error("Unknown Interactor type");
            }

            bool isInteractorRegistered(const std::string& integratorType, const std::string& integratorSubType) const {
                std::pair<std::string, std::string> key(integratorType, integratorSubType);
                return isRegistered<PatchesCreator>(key) || isRegistered<NonBondedPatchesCreator>(key);
            }

            bool isInteractorNonBonded(const std::string& integratorType, const std::string& integratorSubType) const {
                std::pair<std::string, std::string> key(integratorType, integratorSubType);
                return isRegistered<NonBondedPatchesCreator>(key);
            }

            template<typename CreatorType>
            const std::unordered_map<std::pair<std::string, std::string>, CreatorType, PairHash>& getCreators() const {
                return getCreatorsRef<CreatorType>();
            }

    private:
        PatchesFactory() = default;

        template<typename CreatorType>
        static std::unordered_map<std::pair<std::string, std::string>, CreatorType, PairHash>& getCreatorsRef() {
            static std::unordered_map<std::pair<std::string, std::string>, CreatorType, PairHash> creators;
            return creators;
        }

        template<typename CreatorType>
        bool isRegistered(const std::pair<std::string, std::string>& key) const {
            return getCreatorsRef<CreatorType>().find(key) != getCreatorsRef<CreatorType>().end();
        }
    };
}}}

// Macro for registering a patches interactor
#define REGISTER_PATCHES_INTERACTOR(type, subType, ...)                                                                \
    namespace {                                                                                                        \
    struct registerInteractor##type##subType                                                                           \
    {                                                                                                                  \
        registerInteractor##type##subType()                                                                            \
        {                                                                                                              \
            PLUGIN_REGISTRATION_GUARD("PatchesInteractor" + std::string(#type) + std::string(#subType));               \
            uammd::structured::Interactor::PatchesFactory::getInstance()                                               \
                .registerInteractor<uammd::structured::Interactor::PatchesFactory::PatchesCreator>(                    \
                    #type,                                                                                             \
                    #subType,                                                                                          \
                    [](std::shared_ptr<uammd::structured::GlobalData> gd,                                              \
                       std::shared_ptr<uammd::ParticleGroup> pg,                                                       \
                       std::shared_ptr<uammd::structured::GlobalData> patchesGd,                                       \
                       std::shared_ptr<uammd::ParticleGroup> patchesPg,                                                \
                       uammd::structured::DataEntry& data,                                                             \
                       std::string name) -> std::shared_ptr<uammd::Interactor> {                                       \
                        return std::make_shared<__VA_ARGS__>(gd, pg, patchesGd, patchesPg, data, name);                \
                    });                                                                                                \
        }                                                                                                              \
    };                                                                                                                 \
    registerInteractor##type##subType registerInteractor##type##subType##Instance;                                     \
    }

// Macro for registering a non-bonded patches interactor
#define REGISTER_NONBONDED_PATCHES_INTERACTOR(type, subType, ...)                                                      \
    namespace {                                                                                                        \
    struct registerInteractor##type##subType                                                                           \
    {                                                                                                                  \
        registerInteractor##type##subType()                                                                            \
        {                                                                                                              \
            PLUGIN_REGISTRATION_GUARD("NonBondedPatchesInteractor" + std::string(#type) + std::string(#subType));      \
            uammd::structured::Interactor::PatchesFactory::getInstance()                                               \
                .registerInteractor<uammd::structured::Interactor::PatchesFactory::NonBondedPatchesCreator>(           \
                    #type,                                                                                             \
                    #subType,                                                                                          \
                    [](std::shared_ptr<uammd::structured::GlobalData> gd,                                              \
                       std::shared_ptr<uammd::ParticleGroup> pg,                                                       \
                       std::shared_ptr<uammd::structured::GlobalData> patchesGd,                                       \
                       std::shared_ptr<uammd::ParticleGroup> patchesPg,                                                \
                       std::shared_ptr<uammd::structured::VerletConditionalListSetBase> vcls,                          \
                       uammd::structured::DataEntry& data,                                                             \
                       std::string name) -> std::shared_ptr<uammd::Interactor> {                                       \
                        return std::make_shared<__VA_ARGS__>(gd, pg, patchesGd, patchesPg, vcls, data, name);          \
                    });                                                                                                \
        }                                                                                                              \
    };                                                                                                                 \
    registerInteractor##type##subType registerInteractor##type##subType##Instance;                                     \
    }
