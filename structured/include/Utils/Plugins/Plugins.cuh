#pragma once

#include <map>
#include <source_location>
#include <string>
namespace uammd {
namespace structured {
namespace PluginUtils {

/* This function is used to ensure that no two different plugins can have the same identifier.
 * If a plugin has already been registered with the same name but in a different location,
 * it will throw an exception. This is useful to prevent conflicts between plugins.
 */
inline void registrationGuard(std::string identifier, std::string location)
{
    static std::map<std::string, std::string> registeredPlugins;
    const bool alreadyRegistered = registeredPlugins.find(identifier) != registeredPlugins.end();
    if (alreadyRegistered) {
        if (registeredPlugins[identifier] != location) {
            throw std::runtime_error("Plugin with identifier '" + identifier +
                                     "' already registered at a different location: " + registeredPlugins[identifier]);
        }
    }
    registeredPlugins[identifier] = location;
}

} // namespace PluginUtils
} // namespace structured
} // namespace uammd

#define PLUGIN_REGISTRATION_GUARD(identifier)                                                                          \
    {                                                                                                                  \
        auto __location = std::source_location::current();                                                             \
        uammd::structured::PluginUtils::registrationGuard(                                                             \
            identifier, std::string(__location.file_name()) + " at " + std::to_string(__location.line()));             \
    }
