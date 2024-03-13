import json
import logging

def getLoader(condName):

    condType = f"conditions::{condName}"

    whiteSpace = len(f"     loadVerletConditionalListSet_{condName}")*" "

    loaderTemplate=f"""
    std::shared_ptr<uammd::structured::VerletConditionalListSet<{condType}>>
    loadVerletConditionalListSet_{condName}(std::shared_ptr<ExtendedSystem> sys,
{whiteSpace}std::shared_ptr<GlobalData>           gd,
{whiteSpace}std::shared_ptr<ParticleGroup>        pg,
{whiteSpace}std::vector<std::string>      path){{

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<{condType}>  cond =
        std::make_shared<{condType}>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<{condType}>>(gd,pg,data,cond,path.back());
    }}\n"""

    return loaderTemplate

def getGlobalLoader(conditionsNames):

    loaderCondition =""
    for c in conditionsNames:
        loaderTemplate = f"""
        if("{c}" == condition){{
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet {c}",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_{c}(sys,gd,pg,path);
            found = true;
        }}"""
        loaderCondition += loaderTemplate

    whiteSpace = len(f"     loadVerletConditionalListSet")*" "

    availableConditions = ", ".join(conditionsNames)
    globalLoaderTemplate=f"""
    std::shared_ptr<uammd::structured::VerletConditionalListSetBase>
    loadVerletConditionalListSet(std::shared_ptr<ExtendedSystem> sys,
{whiteSpace}std::shared_ptr<GlobalData>    gd,
{whiteSpace}std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
{whiteSpace}std::vector<std::string>       path){{


        DataEntry data = sys->getInput()->getDataEntry(path);
        //Check data type is VerletConditionalListSet
        if(data.getType() != "VerletConditionalListSet"){{
            System::log<System::CRITICAL>("[VerletConditionalListSetLoader] (%s) Data type is not VerletConditionalListSet!",path.back().c_str());
        }}

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        //////////////////////////////////////////////////////////////////////////////

        std::shared_ptr<uammd::structured::VerletConditionalListSetBase> vCondListSet;

        std::string condition = data.getSubType();
        bool found=false;
{loaderCondition}

        if(not found){{
            std::string available = "{availableConditions}";
            System::log<System::CRITICAL>("[VerletConditionalListSetLoader] (%s) Condition '%s' not available. Available conditions are: %s",path.back().c_str(),condition.c_str(),available.c_str());
        }}

        return vCondListSet;

    }}\n"""

    return globalLoaderTemplate

def generateVerletConditionalListSetLoaders(componentsPath,condComponentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")
    output = "/".join(condComponentPath)+"/"+condComponentPath[-1]+"Loaders.cuh"

    with open(componentsPath,"r") as f:
        condComp = json.load(f)
        for path in condComponentPath:
            condComp = condComp[path]

    #Discard first element, should be VerletConditionalListSet for all
    condComp = [c[1] for c in condComp]

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write("#ifndef __VERLET_CONDITIONAL_LIST_SET_LOADERS__\n")
        fout.write("#define __VERLET_CONDITIONAL_LIST_SET_LOADERS__\n\n")

        fout.write("namespace uammd{\n")
        fout.write("namespace structured{\n")
        fout.write("namespace VerletConditionalListSetLoaders{\n")

        for cond in condComp:
            logger.debug("Generating loader for: VerletConditionalListSet %s",cond)

            loader = getLoader(cond)
            fout.write(loader)

        loader = getGlobalLoader(condComp)
        fout.write(loader)

        fout.write("}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""
