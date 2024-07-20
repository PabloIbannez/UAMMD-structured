#pragma once

#include "ThirdParty/json.hpp"

#include "Input/Formats/json.cuh"

#include "System/System.h"
#include "DataStructures/GroupsList/GroupsList.cuh"

#include "Utils/String/StringUtils.cuh"

#include <string>
#include <fstream>

#include <vector>
#include <set>
#include <map>

namespace uammd{
namespace structured{
namespace Input{

    //Forward declaration of Input
    class Input;

    struct entryBase{

        std::string name;
        JSON_TYPE* entry_ptr;

        std::vector<std::string> ref_path;
        Input* ref_ptr;

        std::vector<std::string>              keys;
        std::vector<std::vector<std::string>> usedPaths;

        entryBase(std::string name,JSON_TYPE* entry_ptr,
                  std::vector<std::string> ref_path,Input* ref_ptr);
        ~entryBase();

        bool isInfoAdded(std::string info);

    };

    struct parametersEntry : public entryBase {

        parametersEntry(std::string name,JSON_TYPE* entry_ptr,
                        std::vector<std::string> ref_path,
                        Input* ref_ptr):entryBase(name,entry_ptr,ref_path,ref_ptr){}

        template<typename T>
        T getParameter(std::string parameterName){
            if (!isInfoAdded(parameterName)){
                System::log<System::CRITICAL>("[ParametersEntry] (%s) Parameter %s"
                                              " is not present",
                                              name.c_str(),parameterName.c_str());
            }

            this->usedPaths.push_back({parameterName});
            return entry_ptr->at(parameterName).get<T>();
        }

        template<typename T>
        T getParameter(std::string parameterName, T defaultValue){
            if (!isInfoAdded(parameterName)){
                std::stringstream ss;
                ss << defaultValue;
                System::log<System::DEBUG>("[ParametersEntry] (%s) Parameter \"%s\""
                                           " is not present, using default value: %s",
                                           name.c_str(),parameterName.c_str(),ss.str().c_str());
                return defaultValue;
            } else{
                this->usedPaths.push_back({parameterName});
            }

            return entry_ptr->at(parameterName).get<T>();
        }

        template<typename T>
        void setParameter(std::string parameterName, T value){
            if (!isInfoAdded(parameterName)){
                System::log<System::DEBUG>("[ParametersEntry] (%s) Parameter \"%s\""
                                             " is not present, adding it",
                                             name.c_str(),parameterName.c_str());
                keys.push_back(parameterName);
                (*entry_ptr)[parameterName] = value;
            } else {
                entry_ptr->at(parameterName) = value;
            }

        }

        bool isParameterAdded(std::string parameterName);

        std::map<std::string,JSON_TYPE> getParametersMap();
    };

    struct typeEntry : public entryBase {

        std::vector<std::string> cls = {"Unknow","Unknow"};

        bool isType = false;

        typeEntry(std::string name,
                  JSON_TYPE* entry_ptr,
                  std::vector<std::string> ref_path,
                  Input* ref_ptr);

        std::string getType();
        std::string getSubType();
    };

    struct dataEntry : public typeEntry {

        bool isLabels     = false;
        bool isData       = false;
        bool isParameters = false;

        std::vector<std::string>  labels;
        std::map<std::string,int> labelIndices;
        JSON_TYPE* data;
        std::shared_ptr<parametersEntry> parameters;

        dataEntry(std::string name,
                  JSON_TYPE* entry_ptr,
                  std::vector<std::string> ref_path,
                  Input* ref_ptr);

        std::vector<std::string> getLabels();
        void setLabels(const std::vector<std::string>& newLabels);

        bool isDataAdded(std::string label);
        bool isParameterAdded(std::string parameterName);

        //Data access
        int getDataSize();

        template<typename T>
        std::vector<T> getData(std::string label){
            if (!isData){
                System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                              name.c_str());
            }

            if (!isDataAdded(label)){
                System::log<System::CRITICAL>("[DataEntry] (%s) Label \"%s\" is not present",
                                              name.c_str(),label.c_str());
            }

            std::vector<T> ret;
            for(auto e : data->items()){
                ret.push_back(e.value().at(labelIndices[label]).get<T>());
            }

            this->usedPaths.push_back({"data"});
            return ret;
        }

        std::vector<std::map<std::string,JSON_TYPE>>
        getDataMap(std::vector<std::string>& selLabels);

        std::vector<std::map<std::string,JSON_TYPE>>
        getDataMapIgnoringLabels(std::vector<std::string>& ignoredLabels);

        std::vector<std::map<std::string,JSON_TYPE>>
        getDataMap();

        ///////////////////////////////////////////////////////

        std::shared_ptr<groupsList>
        getDataAsGroupsList(std::string idLabel, std::string groupsLabel,
                            int idMin, int idMax, bool sym = false);

        std::shared_ptr<groupsList>
        getDataAsGroupsList(std::string idLabel, std::string groupsLabel, bool sym = false);

        std::shared_ptr<groupsList>
        getDataAsGroupsList(std::string groupsLabel, bool sym = false);

        ///////////////////////////////////////////////////////

        void setData(uint index,nlohmann::json::initializer_list_t inputData);

        template<typename T>
        void setData(uint index1,uint index2,T inputData){
            if (!isData){
                System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                              name.c_str());
            }

            if (index1 >= data->size()){
                System::log<System::CRITICAL>("[DataEntry] (%s) Index out of range",
                                              name.c_str());
            }

            if (index2 >= data->at(index1).size()){
                //Append
                for(uint i=data->at(index1).size();i<index2;i++){
                    data->at(index1).push_back(nlohmann::json());
                }
                data->at(index1).push_back(inputData);
            } else {
                data->at(index1).at(index2) = inputData;
            }
        }

        ///////////////////////////////////////////////////////

        //Parameter access

        template<typename T>
        T getParameter(std::string parameterName){
            if (!isParameters){
                System::log<System::CRITICAL>("[DataEntry] (%s) Parameters are not present",
                                              name.c_str());
            }

            return parameters->getParameter<T>(parameterName);
        }

        template<typename T>
        T getParameter(std::string parameterName, T defaultValue){
            if (!isParameters){
                return defaultValue;
            }
            return parameters->getParameter<T>(parameterName,defaultValue);
        }

        template<typename T>
        void setParameter(std::string parameterName,T value){
            if (!isParameters){
                System::log<System::DEBUG>("[DataEntry] (%s) Parameters are not present. Adding them.",
                                            name.c_str());

                nlohmann::json param = {};
                (*entry_ptr)[path2jsonPointer({"parameters"})] = param;

                std::vector<std::string> param_path = ref_path;
                param_path.push_back("parameters");

                parameters = std::make_shared<parametersEntry>(name,&entry_ptr->at("parameters"),
                                                               param_path,this->ref_ptr);
                isParameters = true;
                this->usedPaths.push_back({"parameters"});
            }

            parameters->setParameter<T>(parameterName,value);
        }

        std::map<std::string,JSON_TYPE> getParametersMap();
    };

    class Input{

        friend struct entryBase;
        friend struct parametersEntry;
        friend struct typeEntry;
        friend struct dataEntry;

        public:

            using DataType = typename JSON_TYPE;

        private:

            JSON_TYPE* getJSONFromPath(const std::vector<std::string>& path);

            void extractPaths(const JSON_TYPE& j,
                              std::vector<std::string>& currentPath,
                              std::vector<std::vector<std::string>>& paths);

            std::vector<std::vector<std::string>> getAllPaths();

            void addUsedPath(const std::vector<std::string>& path);

            std::shared_ptr<JSON_TYPE>  js;
            std::set<std::vector<std::string>> usedPaths;

        public:

            Input(std::string inputFilePath);

            Input(const JSON_TYPE& jsBuffer);

            std::string getFileExtension();
            std::vector<std::string> getEntriesList(std::vector<std::string> path);

            bool isEntryAdded(const std::vector<std::string>& path);

            void addEntry(const std::vector<std::string>& path,
                          std::string cls,std::string subCls);
            void removeEntry(const std::vector<std::string>& path);

            // Data structures getters

            parametersEntry getParametersEntry(const std::vector<std::string>& path);
            typeEntry       getTypeEntry(const std::vector<std::string>& path);
            dataEntry       getDataEntry(const std::vector<std::string>& path);

            // Check all paths
            bool checkAllPathsUsed();

            // Data structure update

            void write(std::string fileName);
    };

}


using ParametersEntry = Input::parametersEntry;
using TypeEntry       = Input::typeEntry;
using DataEntry       = Input::dataEntry;

}}
