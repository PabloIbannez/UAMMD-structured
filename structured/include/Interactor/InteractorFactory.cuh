#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "Definitions/Hashes.cuh"

#include "Interactor/Interactor.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"

namespace uammd {
namespace structured {
namespace Interactor {

    class InteractorFactory {

        public:

            using InteractorPtr    = std::shared_ptr<uammd::Interactor>;
            using GlobalDataPtr    = std::shared_ptr<GlobalData>;
            using ParticleGroupPtr = std::shared_ptr<ParticleGroup>;
            using VCLSPtr          = std::shared_ptr<VerletConditionalListSetBase>;

            using Creator                 = std::function<InteractorPtr(GlobalDataPtr, ParticleGroupPtr,
                                                                        DataEntry&, std::string)>;

            using NonBondedCreator        = std::function<InteractorPtr(GlobalDataPtr, ParticleGroupPtr,
                                                                        VCLSPtr,
                                                                        DataEntry&, std::string)>;

            using PatchyParticleCreator        = std::function<InteractorPtr(GlobalDataPtr, ParticleGroupPtr,
                                                                             std::vector<std::string>,
                                                                             std::string)>;

            static InteractorFactory& getInstance() {
                static InteractorFactory instance;
                return instance;
            }

            template<typename CreatorType>
            void registerInteractor(const std::string& interactorType, const std::string& interactorSubType, CreatorType creator) {
                System::log<System::MESSAGE>("[InteractorFactory] Registering Interactor in factory: %s, %s",
                                             interactorType.c_str(), interactorSubType.c_str());

                std::pair<std::string, std::string> key(interactorType, interactorSubType);
                if (getCreatorsRef<CreatorType>().find(key) != getCreatorsRef<CreatorType>().end()) {
                    System::log<System::CRITICAL>("[InteractorFactory] Interactor already registered: %s, %s",
                                                  interactorType.c_str(), interactorSubType.c_str());
                    throw std::runtime_error("Interactor already registered");
                }
                getCreatorsRef<CreatorType>()[key] = creator;
            }

            template<typename CreatorType, typename... Args>
            InteractorPtr createInteractor(const std::string& interactorType, const std::string& interactorSubType, Args&&... args) {
                System::log<System::MESSAGE>("[InteractorFactory] Creating Interactor: %s (type: %s, subType: %s)",
                                             std::get<sizeof...(Args) - 1>(std::forward_as_tuple(args...)).c_str(), // name is always the last argument
                                             interactorType.c_str(), interactorSubType.c_str());

                std::pair<std::string, std::string> key(interactorType, interactorSubType);
                auto& creators = getCreatorsRef<CreatorType>();
                auto it = creators.find(key);

                if (it != creators.end()) {
                    return it->second(std::forward<Args>(args)...);
                }

                System::log<System::CRITICAL>("[InteractorFactory] Unknown Interactor type: %s, subType: %s",
                                              interactorType.c_str(), interactorSubType.c_str());
                throw std::runtime_error("Unknown Interactor type");
            }

            bool isInteractorRegistered(const std::string& interactorType, const std::string& interactorSubType) const {
                std::pair<std::string, std::string> key(interactorType, interactorSubType);
                return isRegistered<Creator>(key) || isRegistered<NonBondedCreator>(key) || isRegistered<PatchyParticleCreator>(key);
            }

            bool isInteractorNonBonded(const std::string& interactorType, const std::string& interactorSubType) const {
                std::pair<std::string, std::string> key(interactorType, interactorSubType);
                return isRegistered<NonBondedCreator>(key);
            }

            bool isInteractorPatchyParticle(const std::string& interactorType, const std::string& interactorSubType) const {
                std::pair<std::string, std::string> key(interactorType, interactorSubType);
                return isRegistered<PatchyParticleCreator>(key);
            }

            template<typename CreatorType>
            const std::unordered_map<std::pair<std::string, std::string>, CreatorType, PairHash>& getCreators() const {
                return getCreatorsRef<CreatorType>();
            }

    private:
        InteractorFactory() = default;

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

// Macro for registering a standard interactor
#define __REGISTER_INTERACTOR__(type, subType, ...) \
    namespace { \
        struct registerInteractor##type##subType { \
            registerInteractor##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Interactor::InteractorFactory::getInstance().registerInteractor<uammd::structured::Interactor::InteractorFactory::Creator>( \
                        #type, #subType, \
                        [](std::shared_ptr<uammd::structured::GlobalData> gd, \
                           std::shared_ptr<uammd::ParticleGroup> pg, \
                           uammd::structured::DataEntry&  data, \
                           std::string name) -> std::shared_ptr<uammd::Interactor> { \
                        return std::make_shared<__VA_ARGS__>(gd, pg, data, name); \
                    }); \
                } \
            } \
        }; \
        registerInteractor##type##subType registerInteractor##type##subType##Instance; \
    }

#define REGISTER_BOND_INTERACTOR(type, subType, ...) \
    __REGISTER_INTERACTOR__(type, subType, __VA_ARGS__)

#define REGISTER_SET_INTERACTOR(type, subType, ...) \
    __REGISTER_INTERACTOR__(type, subType, __VA_ARGS__)

#define REGISTER_SINGLE_INTERACTOR(type, subType, ...) \
    __REGISTER_INTERACTOR__(type, subType, __VA_ARGS__)

#define REGISTER_AFM_INTERACTOR(type, subType, ...) \
    __REGISTER_INTERACTOR__(type, subType, __VA_ARGS__)

// Macro for registering a non-bonded interactor
#define REGISTER_NONBONDED_INTERACTOR(type, subType, ...) \
    namespace { \
        struct registerInteractor##type##subType { \
            registerInteractor##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Interactor::InteractorFactory::getInstance().registerInteractor<uammd::structured::Interactor::InteractorFactory::NonBondedCreator>( \
                        #type, #subType, \
                        [](std::shared_ptr<uammd::structured::GlobalData> gd, \
                           std::shared_ptr<uammd::ParticleGroup> pg, \
                           std::shared_ptr<uammd::structured::VerletConditionalListSetBase> vcls, \
                           uammd::structured::DataEntry&  data, \
                           std::string name) -> std::shared_ptr<uammd::Interactor> { \
                        return std::make_shared<__VA_ARGS__>(gd, pg, vcls, data, name); \
                    }); \
                } \
            } \
        }; \
        registerInteractor##type##subType registerInteractor##type##subType##Instance; \
    }

// Macro for registering a patchy particle interactor
#define REGISTER_PATCHY_PARTICLE_INTERACTOR(type, subType, ...) \
    namespace { \
        struct registerInteractor##type##subType { \
            registerInteractor##type##subType() { \
                if (__INCLUDE_LEVEL__ == 0) { \
                    uammd::structured::Interactor::InteractorFactory::getInstance().registerInteractor<uammd::structured::Interactor::InteractorFactory::PatchyParticleCreator>( \
                        #type, #subType, \
                        [](std::shared_ptr<uammd::structured::GlobalData> gd, \
                           std::shared_ptr<uammd::ParticleGroup> pg, \
                           std::vector<std::string> path, \
                           std::string name) -> std::shared_ptr<uammd::Interactor> { \
                        return std::make_shared<__VA_ARGS__>(gd, pg, path, name); \
                    }); \
                } \
            } \
        }; \
        registerInteractor##type##subType registerInteractor##type##subType##Instance; \
    }
