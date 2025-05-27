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

#define PLUGIN_REGISTRATION_GUARD(identifier)                                                                          \
    {                                                                                                                  \
        auto __location = std::source_location::current();                                                             \
        uammd::structured::PluginUtils::registrationGuard(                                                             \
            identifier,                                                                                                \
            __location.file_name());   \
    }
