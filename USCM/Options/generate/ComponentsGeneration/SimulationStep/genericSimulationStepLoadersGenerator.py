import json
import logging

def getIsSimulationStep(simulationSteps):

    isSimulationStep = ""
    for simulationStep in simulationSteps:
        subClass,simulationStepType,simulationStepSubType,f = simulationStep
        isSimulationStep += f"""if(\"{simulationStepType}\" == simulationStepType and \"{simulationStepSubType}\" == simulationStepSubType){{
            return true;
        }}
        """

    isSimulationStepTemplate=f"""
    bool isSimulationStepAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();
        {isSimulationStep}
        return false;

    }}\n
    """
    return isSimulationStepTemplate


def getSimulationStepLoader(simulationSteps):

    logger = logging.getLogger("USCM")

    simLoader=""
    for simulationStep in simulationSteps:
        subClass,simulationStepType,simulationStepSubType,f = simulationStep

        logger.debug("Processing simulations step: %s %s",simulationStepType,simulationStepSubType)

        simLoader += f"""
        if(\"{simulationStepType}\" == simulationStepType and \"{simulationStepSubType}\" == simulationStepSubType){{
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected {simulationStepType}::{simulationStepSubType} simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::{subClass}::{simulationStepSubType}>(pg,integrator,ff,data,path.back());
            found = true;
        }}"""


    loaderTemplate=f"""
    std::shared_ptr<SimulationStep::SimulationStepBase>
    loadSimulationStep(std::shared_ptr<ExtendedSystem> sys,
                       std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField> ff,
                       std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();

        std::shared_ptr<SimulationStep::SimulationStepBase> simulationStep;
        bool found = false;
        {simLoader}

        if(not found){{
            System::log<System::CRITICAL>("[SimulationStepLoader] (%s) Could not find simulationStep %s::%s",
                                            path.back().c_str(),simulationStepType.c_str(),simulationStepSubType.c_str());
        }}

        return simulationStep;

    }}\n
    """
    return loaderTemplate

    loaderTemplate=f"""
    """
    return loaderTemplate

def generateGenericSimulationStepLoader(componentsPath,simStepComponentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    output = "/".join(simStepComponentPath)+"/GenericSimulationStepLoader.cuh"

    with open(componentsPath,"r") as f:
        simStepComp = json.load(f)
        for path in simStepComponentPath:
            simStepComp = simStepComp[path]

    #Merge all keys in the dictionary
    subTypeSimStepComp = []
    for subClass,simStep in simStepComp.items():
        subTypeSimStepComp += [ [subClass]+s for s in simStep]

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write("#ifndef __SIMULATION_STEP_LOADER__\n")
        fout.write("#define __SIMULATION_STEP_LOADER__\n")

        fout.write("namespace uammd{\n")
        fout.write("namespace structured{\n")
        fout.write("namespace SimulationStepLoader{\n")

        logger.debug("Generating isInteractor function")
        logger.debug("Generating generic simulation step loader")

        isSimulationStep = getIsSimulationStep(subTypeSimStepComp)
        loader       = getSimulationStepLoader(subTypeSimStepComp)

        fout.write(isSimulationStep)
        fout.write(loader)

        fout.write("}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""
