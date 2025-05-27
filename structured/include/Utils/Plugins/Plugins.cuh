#pragma once

#include <set>
#include <string>
namespace uammd {
namespace structured {
namespace PluginUtils {

/* This function is used to ensure that a plugin is registered only once.
 * It checks if the plugin identifier is already in the set of registered plugins.
 * If it is, it throws an exception; otherwise, it adds the identifier to the set.
 */
void registrationGuard(std::string identifier)
{
    static std::set<std::string> registeredPlugins;
    if (registeredPlugins.find(identifier) != registeredPlugins.end()) {
        System::log<System::EXCEPTION>("[PluginUtils] Plugin already registered: %s", identifier.c_str());
        throw std::runtime_error("Plugin already registered");
    }
    registeredPlugins.insert(identifier);
    System::log<System::DEBUG>("[PluginUtils] Registering plugin: %s", identifier.c_str());
}

} // namespace PluginUtils
} // namespace structured
} // namespace uammd
