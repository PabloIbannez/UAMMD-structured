#pragma once

#include <set>
#include <source_location>
#include <string>
namespace uammd {
namespace structured {
namespace PluginUtils {

/* This function is used to ensure that a plugin is registered only once.
 * It checks if the plugin identifier is already in the set of registered plugins.
 * If it is, it throws an exception; otherwise, it adds the identifier to the set.
 */
void registrationGuard(std::string identifier, std::string location)
{
    static std::set<std::string> registeredPlugins;
    auto key = identifier + " at " + location;
    if (registeredPlugins.find(key) != registeredPlugins.end()) {
        System::log<System::EXCEPTION>("[PluginUtils] Plugin already registered: %s", key.c_str());
        throw std::runtime_error("Plugin already registered");
    }
    registeredPlugins.insert(key);
    System::log<System::DEBUG>("[PluginUtils] Registering plugin: %s", key.c_str());
}

} // namespace PluginUtils
} // namespace structured
} // namespace uammd
// namespace {
// struct registerInteractorBond3KratkyPorod
// {
//     registerInteractorBond3KratkyPorod()
//     {
//         uammd::structured::PluginUtils::registrationGuard(
//             "Interactor" + std::string("Bond3") + std::string("KratkyPorod"),
//             std::source_location::current().file_name() + ":" +
//             std::to_string(std::source_location::current().line()) +
//                 " in " + std::source_location::current().function_name());
//         ;
//         uammd::structured::Interactor::InteractorFactory::getInstance()
//             .registerInteractor<uammd::structured::Interactor::InteractorFactory::Creator>(
//                 "Bond3",
//                 "KratkyPorod",
//                 [](std::shared_ptr<uammd::structured::GlobalData> gd,
//                    std::shared_ptr<uammd::ParticleGroup> pg,
//                    uammd::structured::DataEntry& data,
//                    std::string name) -> std::shared_ptr<uammd::Interactor> {
//                     return std::make_shared<uammd::structured::Interactor::BondsInteractor<
//                         uammd::structured::Potentials::Bond3::KratkyPorod>>(gd, pg, data, name);
//                 });
//     }
// };
// registerInteractorBond3KratkyPorod registerInteractorBond3KratkyPorodInstance;
// } // namespace

#define PLUGIN_REGISTRATION_GUARD(identifier)                                                                          \
    {                                                                                                                  \
        auto __location = std::source_location::current();                                                             \
        uammd::structured::PluginUtils::registrationGuard(                                                             \
            identifier,                                                                                                \
            __location.file_name());   \
    }
