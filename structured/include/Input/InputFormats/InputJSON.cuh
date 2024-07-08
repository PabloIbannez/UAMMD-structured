#pragma once

#include "System/System.h"

#include "ThirdParty/json.hpp"

#include "DataStructures/GroupsList/GroupsList.cuh"

#include<string>
#include<fstream>

#include<vector>
#include<set>
#include<map>

#define JSON_TYPE nlohmann::json

inline void to_json(JSON_TYPE& j, const float2& p) {
    j = JSON_TYPE{p.x,p.y};
}

inline void from_json(const JSON_TYPE& j, float2& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
}


inline void to_json(JSON_TYPE& j, const float3& p) {
    j = JSON_TYPE{p.x,p.y,p.z};
}

inline void from_json(const JSON_TYPE& j, float3& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
    p.z=float(j[2]);
}

inline void to_json(JSON_TYPE& j, const float4& p) {
    j = JSON_TYPE{p.x,p.y,p.z,p.w};
}

inline void from_json(const JSON_TYPE& j, float4& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
    p.z=float(j[2]);
    p.w=float(j[3]);
}

inline void to_json(JSON_TYPE& j, const double2& p) {
    j = JSON_TYPE{p.x,p.y};
}

inline void from_json(const JSON_TYPE& j, double2& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
}

inline void to_json(JSON_TYPE& j, const double3& p) {
    j = JSON_TYPE{p.x,p.y,p.z};
}

inline void from_json(const JSON_TYPE& j, double3& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
    p.z=double(j[2]);
}

inline void to_json(JSON_TYPE& j, const double4& p) {
    j = JSON_TYPE{p.x,p.y,p.z,p.w};
}

inline void from_json(const JSON_TYPE& j, double4& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
    p.z=double(j[2]);
    p.w=double(j[3]);
}

namespace uammd{
namespace structured{
namespace InputJSON{

    namespace InputJSON_ns{
        inline std::string strip(const std::string& str_in)
        {
            std::string str = str_in;

            if (str.length() == 0) {
                return str;
            }

            auto start_it = str.begin();
            auto end_it = str.rbegin();
            while (std::isspace(*start_it)) {
                ++start_it;
                if (start_it == str.end()) break;
            }
            while (std::isspace(*end_it)) {
                ++end_it;
                if (end_it == str.rend()) break;
            }
            int start_pos = start_it - str.begin();
            int end_pos = end_it.base() - str.begin();
            str = start_pos <= end_pos ? std::string(start_it, end_it.base()) : "";

            return str;
        }

        inline nlohmann::json_pointer<nlohmann::json::basic_json::string_t> path2jsonPointer(const std::vector<std::string>& path){
            std::string path_ptr = "";

            for(const std::string& p : path){
                path_ptr+="/"+p;
            }

            return nlohmann::json_pointer<std::string>(path_ptr);
        }
    }


    class InputJSON{

        public:

            using DataType = typename JSON_TYPE;

        private:

            JSON_TYPE* getJSONFromPath(const std::vector<std::string>& path){

                auto jsonPointer = InputJSON_ns::path2jsonPointer(path);

                JSON_TYPE* jsonBuffer;

                try{
                    jsonBuffer = &js->at(jsonPointer);
                } catch (nlohmann::json::out_of_range& e) {
                    std::string err = path[0];
                    for(uint i=1;i<path.size();i++){err+=std::string(", ") + path[i];}

                    System::log<System::CRITICAL>("[InputJSON] Entry \"%s\""
                                                  " is not present",
                                                  err.c_str());
                }

                return jsonBuffer;
            }

            void extractPaths(const JSON_TYPE& j,
                              std::vector<std::string>& currentPath,
                              std::vector<std::vector<std::string>>& paths) {

                if (j.is_object()) {
                    for (auto it = j.begin(); it != j.end(); ++it) {
                        currentPath.push_back(it.key());
                        paths.push_back(currentPath);
                        extractPaths(it.value(), currentPath, paths);
                        currentPath.pop_back();
                    }
                }
            }

            std::vector<std::vector<std::string>> getAllPaths(){
                std::vector<std::vector<std::string>> paths;
                std::vector<std::string> currentPath;

                extractPaths(*js,currentPath,paths);

                return paths;
            }

            void addUsedPath(const std::vector<std::string>& path){
                std::vector<std::string> tmpPath;
                for (std::string p : path){
                    tmpPath.push_back(p);
                    usedPaths.insert(tmpPath);
                }
            }

            struct entryBase{
                std::string name;
                JSON_TYPE* entry_ptr;

                std::vector<std::string> ref_path;
                InputJSON* ref_ptr;

                std::vector<std::string>              keys;
                std::vector<std::vector<std::string>> usedPaths;

                entryBase(std::string name,JSON_TYPE* entry_ptr,
                          std::vector<std::string> ref_path,InputJSON* ref_ptr):
                          name(name),entry_ptr(entry_ptr),
                          ref_path(ref_path),ref_ptr(ref_ptr){

                    for(auto itm : entry_ptr->items()){keys.push_back(itm.key());}
                }

                ~entryBase(){
                    for(auto p : usedPaths){
                        p.insert(p.begin(),ref_path.begin(),ref_path.end());
                        ref_ptr->addUsedPath(p);
                    }
                }

                bool isInfoAdded(std::string info){
                    return std::find(keys.begin(),keys.end(),info) != keys.end();
                }
            };

            std::shared_ptr<JSON_TYPE>  js;
            std::set<std::vector<std::string>> usedPaths;

        public:

            InputJSON(std::string inputFilePath){
                std::ifstream ifs(inputFilePath);
                JSON_TYPE jsBuffer = JSON_TYPE::parse(ifs,nullptr,true,true); //Ignore comments

                js = std::make_shared<JSON_TYPE>(jsBuffer);

                System::log<System::MESSAGE>("[InputJSON] File %s loaded",inputFilePath.c_str());
            }

            InputJSON(const JSON_TYPE& jsBuffer){
                js = std::make_shared<JSON_TYPE>(jsBuffer);
            }

            std::string getFileExtension(){
                return ".json";
            }

            std::vector<std::string> getEntriesList(std::vector<std::string> path){

                JSON_TYPE* jsonBuffer = getJSONFromPath(path);

                std::vector<std::string> entries;

                for(auto e : jsonBuffer->items()){
                    entries.push_back(e.key());
                }

                return entries;
            }

            bool isEntryAdded(const std::vector<std::string>& path){
                auto jsonPointer = InputJSON_ns::path2jsonPointer(path);
                return js->contains(jsonPointer);
            }

            void addEntry(const std::vector<std::string>& path,std::string cls,std::string subCls){
                {
                    std::string msg = path[0];
                    for(uint i=1;i<path.size();i++){msg+=std::string(", ") + path[i];}

                    System::log<System::DEBUG>("[InputJSON] Adding entry %s",msg.c_str());
                }

                if(!isEntryAdded(path)){
                    nlohmann::json typeInfo = {{"type",{cls,subCls}}};
                    (*js.get())[InputJSON_ns::path2jsonPointer(path)] = typeInfo;
                } else {
                    std::string err = path[0];
                    for(uint i=1;i<path.size();i++){err+=std::string(", ") + path[i];}

                    System::log<System::CRITICAL>("[InputJSON] Entry \"%s\""
                                                  " is already present",
                                                  err.c_str());
                }
            }

            void removeEntry(const std::vector<std::string>& path){
                {
                    std::string msg = path[0];
                    for(uint i=1;i<path.size();i++){msg+=std::string(", ") + path[i];}

                    System::log<System::DEBUG>("[InputJSON] Removing entry %s",msg.c_str());
                }

                if(isEntryAdded(path)){
                    auto pathBuffer = path;

                    std::string toRemove = pathBuffer.back();
                    pathBuffer.pop_back();

                    auto jsonPointer = InputJSON_ns::path2jsonPointer(pathBuffer);

                    auto& entry = js->at(jsonPointer);
                    entry.erase(toRemove);
                } else {
                    std::string err = path[0];
                    for(uint i=1;i<path.size();i++){err+=std::string(", ") + path[i];}

                    System::log<System::CRITICAL>("[InputJSON] Entry \"%s\""
                                                  " is not present",
                                                  err.c_str());
                }
            }

            struct parametersEntry : public entryBase {

                parametersEntry(std::string name,JSON_TYPE* entry_ptr,
                                std::vector<std::string> ref_path,
                                InputJSON* ref_ptr):entryBase(name,entry_ptr,ref_path,ref_ptr){}

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

                bool isParameterAdded(std::string parameterName){
                    return isInfoAdded(parameterName);
                }

                std::map<std::string,nlohmann::json> getParametersMap(){

                    std::map<std::string,nlohmann::json> tmp;
                    for(std::string p : keys){
                        tmp[p] = entry_ptr->at(p);
                        this->usedPaths.push_back({p});
                    }

                    return tmp;

                }

            };

            struct typeEntry : public entryBase {

                std::vector<std::string> cls = {"Unknow","Unknow"};

                bool isType = false;

                typeEntry(std::string name,
                          JSON_TYPE* entry_ptr,
                          std::vector<std::string> ref_path,
                          InputJSON* ref_ptr):entryBase(name,entry_ptr,ref_path,ref_ptr){
                    if(std::find(keys.begin(),keys.end(),"type") != keys.end()){
                        isType = true;

                        cls = entry_ptr->at("type").get<std::vector<std::string>>();
                    }
                }

                std::string getType(){
                    if (!isType){
                        System::log<System::CRITICAL>("[TypeEntry] (%s) Type is not present",
                                                      name.c_str());
                    }

                    this->usedPaths.push_back({"type"});
                    return InputJSON_ns::strip(cls[0]);
                }

                std::string getSubType(){
                    if (!isType){
                        System::log<System::CRITICAL>("[TypeEntry] (%s) Type is not present",
                                                      name.c_str());
                    }

                    this->usedPaths.push_back({"type"});
                    return InputJSON_ns::strip(cls[1]);
                }

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
                          InputJSON* ref_ptr):typeEntry(name,entry_ptr,ref_path,ref_ptr){

                    if(std::find(keys.begin(),keys.end(),"labels") != keys.end())    {isLabels = true;}
                    if(std::find(keys.begin(),keys.end(),"data") != keys.end())      {isData = true;}
                    if(std::find(keys.begin(),keys.end(),"parameters") != keys.end()){isParameters = true;}

                    if(isLabels)    {
                        labels = entry_ptr->at("labels").get<std::vector<std::string>>();
                        for(uint i = 0; i < labels.size(); i++){
                            labelIndices[labels[i]] = i;
                        }
                        this->usedPaths.push_back({"labels"});
                    }
                    if(isData)      {data = &entry_ptr->at("data");}
                    if(isParameters){std::vector<std::string> param_path = ref_path;
                                     param_path.push_back("parameters");

                                     parameters = std::make_shared<parametersEntry>(name,&entry_ptr->at("parameters"),
                                                                                    param_path,this->ref_ptr);

                                     this->usedPaths.push_back({"parameters"});
                    }
                }

                std::vector<std::string> getLabels(){
                    if (!isLabels){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                                      name.c_str());
                    }

                    return labels;
                }

                void setLabels(const std::vector<std::string>& newLabels){
                    if (!isLabels){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                                      name.c_str());
                    }

                    labels = newLabels;

                    this->entry_ptr->at("labels") = newLabels;
                }

                bool isDataAdded(std::string label){
                    if (!isLabels){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                                      name.c_str());
                    }

                    return std::find(labels.begin(),labels.end(),label) != labels.end();
                }

                bool isParameterAdded(std::string parameterName){
                    if (!isParameters){
                        return false;
                    }

                    return parameters->isInfoAdded(parameterName);
                }

                //Data access

                int getDataSize(){
                    if (!isData){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                                      name.c_str());
                    }

                    this->usedPaths.push_back({"data"});
                    return data->size();
                }

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

                std::vector<std::map<std::string,nlohmann::json>> getDataMap(std::vector<std::string>& selLabels){
                    if (!isData){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                                      name.c_str());
                    }

                    for(auto l : selLabels){
                        if (!isDataAdded(l)){
                            System::log<System::CRITICAL>("[DataEntry] (%s) Label \"%s\" is not present",
                                                          name.c_str(),l.c_str());
                        }
                    }

                    std::vector<std::map<std::string,nlohmann::json>> ret;
                    for(auto e : data->items()){
                        std::map<std::string,nlohmann::json> tmp;
                        for(auto l : selLabels){
                            tmp[l] = e.value().at(labelIndices[l]);
                        }
                        ret.push_back(tmp);
                    }

                    this->usedPaths.push_back({"data"});
                    return ret;
                }

                std::vector<std::map<std::string,nlohmann::json>> getDataMapIgnoringLabels(std::vector<std::string>& ignoredLabels){

                    std::vector<std::string> selLabels;

                    for(auto l : labels){
                        if (std::find(ignoredLabels.begin(),ignoredLabels.end(),l) == ignoredLabels.end()){
                            selLabels.push_back(l);
                        }
                    }

                    this->usedPaths.push_back({"data"});
                    return getDataMap(selLabels);
                }

                std::vector<std::map<std::string,nlohmann::json>> getDataMap(){
                    if (!isLabels){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                                      name.c_str());
                    }

                    this->usedPaths.push_back({"data"});
                    return getDataMap(labels);
                }

                ///////////////////////////////////////////////////////

                std::shared_ptr<groupsList> getDataAsGroupsList(std::string idLabel, std::string groupsLabel, int idMin, int idMax, bool sym = false){

                    std::vector<int>              idVector     = getData<int>(idLabel);
                    std::vector<std::vector<int>> groupsVector = getData<std::vector<int>>(groupsLabel);

                    //Complete the id and groups vectors
                    for(int i=idMin;i<=idMax;i++){
                        if(std::find(idVector.begin(),idVector.end(),i) == idVector.end()){
                            //Not found
                            idVector.push_back(i);
                            groupsVector.push_back(std::vector<int>());
                            System::log<System::DEBUG>("[DataEntry] (%s) Added id %d to the data when reading data as groups list",
                                                          name.c_str(),i);
                        }
                    }

                    std::shared_ptr<groupsList> ret = std::make_shared<groupsList>(idVector,groupsVector,sym);

                    this->usedPaths.push_back({"data"});
                    return ret;
                }

                std::shared_ptr<groupsList> getDataAsGroupsList(std::string idLabel, std::string groupsLabel, bool sym = false){

                    std::vector<int>              idVector     = getData<int>(idLabel);
                    std::vector<std::vector<int>> groupsVector = getData<std::vector<int>>(groupsLabel);

                    std::shared_ptr<groupsList> ret = std::make_shared<groupsList>(idVector,groupsVector,sym);

                    this->usedPaths.push_back({"data"});
                    return ret;
                }

                std::shared_ptr<groupsList> getDataAsGroupsList(std::string groupsLabel, bool sym = false){

                    std::vector<std::vector<int>> groupsVector = getData<std::vector<int>>(groupsLabel);

                    std::shared_ptr<groupsList> ret = std::make_shared<groupsList>(groupsVector,sym);

                    this->usedPaths.push_back({"data"});
                    return ret;
                }

                ///////////////////////////////////////////////////////

                void setData(uint index,nlohmann::json::initializer_list_t inputData){
                    if (!isData){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                                      name.c_str());
                    }

                    if (index >= data->size()){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Index out of range",
                                                      name.c_str());
                    }

                    if (inputData.size() != labels.size()){
                        System::log<System::CRITICAL>("[DataEntry] (%s) Input data size does not match labels size",
                                                      name.c_str());
                    }

                    data->at(index) = nlohmann::json(inputData);
                }

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
                        (*entry_ptr)[InputJSON_ns::path2jsonPointer({"parameters"})] = param;

                        std::vector<std::string> param_path = ref_path;
                        param_path.push_back("parameters");

                        parameters = std::make_shared<parametersEntry>(name,&entry_ptr->at("parameters"),
                                                                       param_path,this->ref_ptr);
                        isParameters = true;
                        this->usedPaths.push_back({"parameters"});
                    }

                    parameters->setParameter<T>(parameterName,value);
                }

                std::map<std::string,nlohmann::json> getParametersMap(){
                    return parameters->getParametersMap();
                }
            };

            // Data structures getters

            parametersEntry getParametersEntry(const std::vector<std::string>& path){

                this->addUsedPath(path);

                JSON_TYPE* jsonBuffer = getJSONFromPath(path);

                return parametersEntry(path.back(),jsonBuffer,path,this);
            }

            typeEntry getTypeEntry(const std::vector<std::string>& path){

                this->addUsedPath(path);

                JSON_TYPE* jsonBuffer = getJSONFromPath(path);

                return typeEntry(path.back(),jsonBuffer,path,this);
            }

            dataEntry getDataEntry(const std::vector<std::string>& path){

                this->addUsedPath(path);

                JSON_TYPE* jsonBuffer = getJSONFromPath(path);

                return dataEntry(path.back(),jsonBuffer,path,this);
            }

            // Check all paths
            bool checkAllPathsUsed(){
                bool allPathsUsed = true;

                std::vector<std::vector<std::string>> allPaths = getAllPaths();

                for(auto k : allPaths){
                    if (std::find(usedPaths.begin(),usedPaths.end(),k) == usedPaths.end()){
                        allPathsUsed = false;

                        std::string path = "{";
                        for(auto s : k){
                            path += s + ",";
                        }
                        //Remove last comma
                        path.pop_back();
                        path += "}";

                        System::log<System::WARNING>("[InputJSON] Path not used: %s",path.c_str());
                    }
                }

                return allPathsUsed;
            }

            // Data structure update

            void write(std::string fileName){
                //Check if fileName has extension, else add it
                if (fileName.find(".") == std::string::npos){
                    System::log<System::DEBUG>("[InputJSON] Adding extension \".json\" to output file name");
                    fileName += this->getFileExtension();
                }

                std::string outputBuffer = js->dump();
                std::ofstream out(fileName);
                //Write outputBuffer to file, char by char
                uint dictLevel = 0;
                uint listLevel = 0;

                std::string indent = "  ";
                for (uint i = 0 ; i < outputBuffer.size(); i++){

                    char& c = outputBuffer[i];

                    char& prevC = std::string(" ")[0];
                    char& nextC = std::string(" ")[0];

                    if (i > 0){
                        prevC = outputBuffer[i-1];
                    }
                    if (i < outputBuffer.size()-1){
                        nextC = outputBuffer[i+1];
                    }

                    if (c == '{'){dictLevel++;}
                    if (c == '}'){dictLevel--;}

                    if (c == '['){listLevel++;}
                    if (c == ']'){listLevel--;}

                    std::string currentIndent = std::string(dictLevel*indent.size(),' ');

                    bool written = false;
                    if (c == ',' and (listLevel == 0 or (listLevel == 1 and prevC == ']'))){
                        out << c << std::endl << currentIndent ;
                        written = true;
                    }

                    if (c == '[' and prevC == ':' and nextC == '['){
                        out << c << std::endl << currentIndent ;
                        written = true;
                    }

                    if (c == ']' and nextC == '}'){
                        out << c ;
                        written = true;
                    }

                    if (c == '{'){
                        out << c << std::endl << currentIndent ;
                        written = true;
                    }

                    if (c == '}'){
                        out << std::endl << currentIndent << c ;
                        written = true;
                    }

                    if (!written){
                        out << c;
                    }
                }
            }
    };


}}}
