import json
import logging

def getIsIntegrator(integrators):

    isIntegrator = ""
    for integrator in integrators:
        _,integratorType,integratorSubType,file = integrator
        isIntegrator += f"""if(\"{integratorType}\" == integratorType and \"{integratorSubType}\" == integratorSubType){{
            return true;
        }}
        """

    isIntegratorTemplate=f"""
    bool isIntegratorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string integratorType    = data.getType();
        std::string integratorSubType = data.getSubType();
        {isIntegrator}
        return false;

    }}\n
    """
    return isIntegratorTemplate


def getIntegratorLoader(integrators):

    intLoader=""
    for integrator in integrators:
        integratorSubClass,integratorType,integratorSubType,file = integrator

        intLoader += f"""
        if(\"{integratorType}\" == integratorType and \"{integratorSubType}\" == integratorSubType){{
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected {integratorSubClass}::{integratorType}::{integratorSubType} integrator",path.back().c_str());
            integrator = std::make_shared<{integratorSubClass}::{integratorType}::{integratorSubType}>(gd,pg,data,path.back());
            found = true;
        }}"""


    loaderTemplate=f"""
    std::shared_ptr<typename uammd::Integrator>
    loadIntegrator(std::shared_ptr<ExtendedSystem> sys,
                   std::shared_ptr<GlobalData>     gd,
                   std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                   std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string integratorType    = data.getType();
        std::string integratorSubType = data.getSubType();

        std::shared_ptr<typename uammd::Integrator> integrator;
        bool found = false;
        {intLoader}

        if(not found){{
            System::log<System::CRITICAL>("[IntegratorLoader] (%s) Could not find integrator %s::%s",
                                            path.back().c_str(),integratorType.c_str(),integratorSubType.c_str());
        }}

        return integrator;

    }}\n
    """
    return loaderTemplate

    loaderTemplate=f"""
    """
    return loaderTemplate

def generateGenericIntegratorLoader(componentsPath,integratorComponentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    output = "/".join(integratorComponentPath)+"/"+integratorComponentPath[-1]+"Loaders.cuh"

    with open(componentsPath,"r") as f:
        inteComp = json.load(f)
        for path in integratorComponentPath:
            inteComp = inteComp[path]

    #Merge all keys in the dictionary
    inteComp = [[key]+decl for key,value in inteComp.items() for decl in value]

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write("#ifndef __INTEGRATOR_LOADER__\n")
        fout.write("#define __INTEGRATOR_LOADER__\n")

        fout.write("namespace uammd{\n")
        fout.write("namespace structured{\n")
        fout.write("namespace IntegratorLoader{\n")

        logger.debug("Generating isInteractor function")
        logger.debug("Generating generic integrator loader")

        isIntegrator = getIsIntegrator(inteComp)
        loader       = getIntegratorLoader(inteComp)

        fout.write(isIntegrator)
        fout.write(loader)
        fout.write("}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""

