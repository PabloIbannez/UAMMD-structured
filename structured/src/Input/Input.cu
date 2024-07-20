#include "Input/Input.cuh"

namespace uammd{
namespace structured{
namespace Input{

    //EntryBase
    entryBase::entryBase(std::string name,JSON_TYPE* entry_ptr,
                         std::vector<std::string> ref_path,Input* ref_ptr):
                         name(name),entry_ptr(entry_ptr),
                         ref_path(ref_path),ref_ptr(ref_ptr){

            for(auto itm : entry_ptr->items()){keys.push_back(itm.key());}
    }

    entryBase::~entryBase(){
        for(auto p : usedPaths){
            p.insert(p.begin(),ref_path.begin(),ref_path.end());
            ref_ptr->addUsedPath(p);
        }
    }

    bool entryBase::isInfoAdded(std::string info){
        return std::find(keys.begin(),keys.end(),info) != keys.end();
    }

    //ParametersEntry
    bool parametersEntry::isParameterAdded(std::string parameterName){
        return isInfoAdded(parameterName);
    }

    std::map<std::string,JSON_TYPE> parametersEntry::getParametersMap(){

        std::map<std::string,nlohmann::json> tmp;
        for(std::string p : keys){
            tmp[p] = entry_ptr->at(p);
            this->usedPaths.push_back({p});
        }

        return tmp;
    }

    //TypeEntry
    typeEntry::typeEntry(std::string name,
                         JSON_TYPE* entry_ptr,
                         std::vector<std::string> ref_path,
                         Input* ref_ptr):entryBase(name,entry_ptr,ref_path,ref_ptr){
        if(std::find(keys.begin(),keys.end(),"type") != keys.end()){
            isType = true;

            cls = entry_ptr->at("type").get<std::vector<std::string>>();
        }
    }

    std::string typeEntry::getType(){
        if (!isType){
            System::log<System::CRITICAL>("[TypeEntry] (%s) Type is not present",
                                          name.c_str());
        }

        this->usedPaths.push_back({"type"});
        return StringUtils::strip(cls[0]);
    }

    std::string typeEntry::getSubType(){
        if (!isType){
            System::log<System::CRITICAL>("[TypeEntry] (%s) Type is not present",
                                          name.c_str());
        }

        this->usedPaths.push_back({"type"});
        return StringUtils::strip(cls[1]);
    }

    //DataEntry
    dataEntry::dataEntry(std::string name,
                         JSON_TYPE* entry_ptr,
                         std::vector<std::string> ref_path,
                         Input* ref_ptr):typeEntry(name,entry_ptr,ref_path,ref_ptr){

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

    std::vector<std::string> dataEntry::getLabels(){
        if (!isLabels){
            System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                          name.c_str());
        }

        return labels;
    }

    void dataEntry::setLabels(const std::vector<std::string>& newLabels){
        if (!isLabels){
            System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                          name.c_str());
        }

        labels = newLabels;

        this->entry_ptr->at("labels") = newLabels;
    }

    bool dataEntry::isDataAdded(std::string label){
        if (!isLabels){
            System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                          name.c_str());
        }

        return std::find(labels.begin(),labels.end(),label) != labels.end();
    }

    bool dataEntry::isParameterAdded(std::string parameterName){
        if (!isParameters){
            return false;
        }

        return parameters->isInfoAdded(parameterName);
    }

    int dataEntry::getDataSize(){
        if (!isData){
            System::log<System::CRITICAL>("[DataEntry] (%s) Data is not present",
                                          name.c_str());
        }

        this->usedPaths.push_back({"data"});
        return data->size();
    }

    std::vector<std::map<std::string,JSON_TYPE>>
    dataEntry::getDataMap(std::vector<std::string>& selLabels){
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

        std::vector<std::map<std::string,JSON_TYPE>> ret;
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

    std::vector<std::map<std::string,JSON_TYPE>>
    dataEntry::getDataMapIgnoringLabels(std::vector<std::string>& ignoredLabels){

        std::vector<std::string> selLabels;

        for(auto l : labels){
            if (std::find(ignoredLabels.begin(),ignoredLabels.end(),l) == ignoredLabels.end()){
                selLabels.push_back(l);
            }
        }

        this->usedPaths.push_back({"data"});
        return getDataMap(selLabels);
    }

    std::vector<std::map<std::string,JSON_TYPE>>
    dataEntry::getDataMap(){
        if (!isLabels){
            System::log<System::CRITICAL>("[DataEntry] (%s) Labels are not present",
                                          name.c_str());
        }

        this->usedPaths.push_back({"data"});
        return getDataMap(labels);
    }

    ///////////////////////////////////////////////////////

    std::shared_ptr<groupsList>
    dataEntry::getDataAsGroupsList(std::string idLabel, std::string groupsLabel, int idMin, int idMax, bool sym){

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

    std::shared_ptr<groupsList>
    dataEntry::getDataAsGroupsList(std::string idLabel, std::string groupsLabel, bool sym){

        std::vector<int>              idVector     = getData<int>(idLabel);
        std::vector<std::vector<int>> groupsVector = getData<std::vector<int>>(groupsLabel);

        std::shared_ptr<groupsList> ret = std::make_shared<groupsList>(idVector,groupsVector,sym);

        this->usedPaths.push_back({"data"});
        return ret;
    }

    std::shared_ptr<groupsList>
    dataEntry::getDataAsGroupsList(std::string groupsLabel, bool sym){

        std::vector<std::vector<int>> groupsVector = getData<std::vector<int>>(groupsLabel);

        std::shared_ptr<groupsList> ret = std::make_shared<groupsList>(groupsVector,sym);

        this->usedPaths.push_back({"data"});
        return ret;
    }

    ///////////////////////////////////////////////////////

    void dataEntry::setData(uint index,nlohmann::json::initializer_list_t inputData){
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

    //Parameter access
    std::map<std::string,JSON_TYPE>  dataEntry::getParametersMap(){
        return parameters->getParametersMap();
    }

    JSON_TYPE* Input::getJSONFromPath(const std::vector<std::string>& path){

        auto jsonPointer = path2jsonPointer(path);

        JSON_TYPE* jsonBuffer;

        try{
            jsonBuffer = &js->at(jsonPointer);
        } catch (nlohmann::json::out_of_range& e) {
            std::string err = path[0];
            for(uint i=1;i<path.size();i++){err+=std::string(", ") + path[i];}

            System::log<System::CRITICAL>("[Input] Entry \"%s\""
                                          " is not present",
                                          err.c_str());
        }

        return jsonBuffer;
    }

    void Input::extractPaths(const JSON_TYPE& j,
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

    std::vector<std::vector<std::string>> Input::getAllPaths(){
        std::vector<std::vector<std::string>> paths;
        std::vector<std::string> currentPath;

        extractPaths(*js,currentPath,paths);

        return paths;
    }

    void Input::addUsedPath(const std::vector<std::string>& path){
        std::vector<std::string> tmpPath;
        for (std::string p : path){
            tmpPath.push_back(p);
            usedPaths.insert(tmpPath);
        }
    }

    Input::Input(std::string inputFilePath){

        std::ifstream ifs(inputFilePath);
        JSON_TYPE jsBuffer = JSON_TYPE::parse(ifs,nullptr,true,true); //Ignore comments

        js = std::make_shared<JSON_TYPE>(jsBuffer);

        System::log<System::MESSAGE>("[Input] File %s loaded",inputFilePath.c_str());
    }

    Input::Input(const JSON_TYPE& jsBuffer){
        js = std::make_shared<JSON_TYPE>(jsBuffer);
    }

    std::string Input::getFileExtension(){
        return ".json";
    }

    std::vector<std::string> Input::getEntriesList(std::vector<std::string> path){

        JSON_TYPE* jsonBuffer = getJSONFromPath(path);

        std::vector<std::string> entries;

        for(auto e : jsonBuffer->items()){
            entries.push_back(e.key());
        }

        return entries;
    }

    bool Input::isEntryAdded(const std::vector<std::string>& path){
        auto jsonPointer = path2jsonPointer(path);
        return js->contains(jsonPointer);
    }

    void Input::addEntry(const std::vector<std::string>& path,
                         std::string cls,std::string subCls){
        {
            std::string msg = path[0];
            for(uint i=1;i<path.size();i++){msg+=std::string(", ") + path[i];}

            System::log<System::DEBUG>("[InputJSON] Adding entry %s",msg.c_str());
        }

        if(!isEntryAdded(path)){
            nlohmann::json typeInfo = {{"type",{cls,subCls}}};
            (*js.get())[path2jsonPointer(path)] = typeInfo;
        } else {
            std::string err = path[0];
            for(uint i=1;i<path.size();i++){err+=std::string(", ") + path[i];}

            System::log<System::CRITICAL>("[InputJSON] Entry \"%s\""
                                          " is already present",
                                          err.c_str());
        }
    }

    void Input::removeEntry(const std::vector<std::string>& path){
        {
            std::string msg = path[0];
            for(uint i=1;i<path.size();i++){msg+=std::string(", ") + path[i];}

            System::log<System::DEBUG>("[InputJSON] Removing entry %s",msg.c_str());
        }

        if(isEntryAdded(path)){
            auto pathBuffer = path;

            std::string toRemove = pathBuffer.back();
            pathBuffer.pop_back();

            auto jsonPointer = path2jsonPointer(pathBuffer);

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


    // Data structures getters

    parametersEntry Input::getParametersEntry(const std::vector<std::string>& path){

        this->addUsedPath(path);

        JSON_TYPE* jsonBuffer = getJSONFromPath(path);

        return parametersEntry(path.back(),jsonBuffer,path,this);
    }

    typeEntry Input::getTypeEntry(const std::vector<std::string>& path){

        this->addUsedPath(path);

        JSON_TYPE* jsonBuffer = getJSONFromPath(path);

        return typeEntry(path.back(),jsonBuffer,path,this);
    }

    dataEntry Input::getDataEntry(const std::vector<std::string>& path){

        this->addUsedPath(path);

        JSON_TYPE* jsonBuffer = getJSONFromPath(path);

        return dataEntry(path.back(),jsonBuffer,path,this);
    }

    // Check all paths
    bool Input::checkAllPathsUsed(){
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

    void Input::write(std::string fileName){
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

}}}
