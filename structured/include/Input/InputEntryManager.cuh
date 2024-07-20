#pragma once

#include "System/ExtendedSystem.cuh"

namespace uammd{
namespace structured{

    class InputEntryManager{

        public:

            struct entryInfo{

                std::string name;

                std::vector<std::string> path;

                std::string entryType;
                std::string entrySubType;

                bool used = false;

            };

        private:

            std::shared_ptr<ExtendedSystem>      sys;

            std::vector<std::string> path;

            //////////////////////////////

            std::map<std::string,entryInfo> entriesInfo;

            //////////////////////////////

            void loadEntriesInfo();

        public:

            InputEntryManager(std::shared_ptr<ExtendedSystem> sys,
                              std::vector<std::string> path);

            bool isEntryPresent(std::string name);

            bool isEntrySubTypePresent(std::string entryType,std::string entrySubType);

            bool isEntryClassPresent(std::string entryType);

            entryInfo& getEntryInfo(std::string name);

            std::map<std::string,entryInfo>& getEntriesInfo();

            std::vector<std::string> getEntriesBySubType(std::string entryType,
                                                         std::string entrySubType);
            std::vector<std::string> getEntriesByClass(std::string entryType);

            void checkEntriesUsed();
    };

}}
