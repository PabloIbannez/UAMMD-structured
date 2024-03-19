#ifndef __SIMULATION_STEP_LOADER__
#define __SIMULATION_STEP_LOADER__
namespace uammd{
namespace structured{
namespace SimulationStepLoader{

    bool isSimulationStepAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();
        if("WriteStep" == simulationStepType and "WriteStep" == simulationStepSubType){
            return true;
        }
        if("WriteStep" == simulationStepType and "WritePatchyParticlesStep" == simulationStepSubType){
            return true;
        }
        if("UtilsStep" == simulationStepType and "InfoStep" == simulationStepSubType){
            return true;
        }
        if("UtilsStep" == simulationStepType and "SortStep" == simulationStepSubType){
            return true;
        }
        if("FlowControl" == simulationStepType and "AFMMaxForce" == simulationStepSubType){
            return true;
        }
        if("FlowControl" == simulationStepType and "LambdaCycle" == simulationStepSubType){
            return true;
        }
        if("FlowControl" == simulationStepType and "LambdaActivation" == simulationStepSubType){
            return true;
        }
        if("ParticlesListMeasure" == simulationStepType and "DistancesMeasure" == simulationStepSubType){
            return true;
        }
        if("ParticlesListMeasure" == simulationStepType and "AnglesMeasure" == simulationStepSubType){
            return true;
        }
        if("ParticlesListMeasure" == simulationStepType and "DihedralsMeasure" == simulationStepSubType){
            return true;
        }
        if("MechanicalMeasure" == simulationStepType and "StressMeasure" == simulationStepSubType){
            return true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "ThermodynamicQuantityMeasure" == simulationStepSubType){
            return true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "ThermodynamicIntegration" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanRadius" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "GyrationRadius" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanSquareDisplacement" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanAngularCorrelation" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "CenterOfMassPosition" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "Height" == simulationStepSubType){
            return true;
        }
        if("TopologicalMeasures" == simulationStepType and "PatchPolymers" == simulationStepSubType){
            return true;
        }
        if("ParticlesListMeasure" == simulationStepType and "PotentialMeasure" == simulationStepSubType){
            return true;
        }
        if("MagneticMeasure" == simulationStepType and "MeasureTotalMagnetization" == simulationStepSubType){
            return true;
        }
        if("MagneticMeasure" == simulationStepType and "MeasureMeanMagnetization" == simulationStepSubType){
            return true;
        }
        if("ExperimentMeasures" == simulationStepType and "AFMMeasure" == simulationStepSubType){
            return true;
        }
        if("GeometricalMeasure" == simulationStepType and "DistanceBetweenCentersOfMass" == simulationStepSubType){
            return true;
        }
        if("ParticlesListMeasure" == simulationStepType and "ContactsMeasure" == simulationStepSubType){
            return true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "InteractorsListEnergyMeasure" == simulationStepSubType){
            return true;
        }
        if("MechanicalMeasure" == simulationStepType and "HessianMeasure" == simulationStepSubType){
            return true;
        }
        if("MechanicalMeasure" == simulationStepType and "PairwiseForceMeasure" == simulationStepSubType){
            return true;
        }
        if("MechanicalMeasure" == simulationStepType and "ForceBetweenSetsMeasure" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasure" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasureFromMobility" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VAFMMeasure" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasure" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasureFromMobility" == simulationStepSubType){
            return true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VAFMMeasure" == simulationStepSubType){
            return true;
        }
        
        return false;

    }

    
    std::shared_ptr<SimulationStep::SimulationStepBase>
    loadSimulationStep(std::shared_ptr<ExtendedSystem> sys,
                       std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField> ff,
                       std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();

        std::shared_ptr<SimulationStep::SimulationStepBase> simulationStep;
        bool found = false;
        
        if("WriteStep" == simulationStepType and "WriteStep" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected WriteStep::WriteStep simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationOutput::WriteStep>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("WriteStep" == simulationStepType and "WritePatchyParticlesStep" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected WriteStep::WritePatchyParticlesStep simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationOutput::WritePatchyParticlesStep>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("UtilsStep" == simulationStepType and "InfoStep" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected UtilsStep::InfoStep simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationUtils::InfoStep>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("UtilsStep" == simulationStepType and "SortStep" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected UtilsStep::SortStep simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationUtils::SortStep>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("FlowControl" == simulationStepType and "AFMMaxForce" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected FlowControl::AFMMaxForce simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationUtils::AFMMaxForce>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("FlowControl" == simulationStepType and "LambdaCycle" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected FlowControl::LambdaCycle simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationUtils::LambdaCycle>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("FlowControl" == simulationStepType and "LambdaActivation" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected FlowControl::LambdaActivation simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationUtils::LambdaActivation>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ParticlesListMeasure" == simulationStepType and "DistancesMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ParticlesListMeasure::DistancesMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::DistancesMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ParticlesListMeasure" == simulationStepType and "AnglesMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ParticlesListMeasure::AnglesMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::AnglesMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ParticlesListMeasure" == simulationStepType and "DihedralsMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ParticlesListMeasure::DihedralsMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::DihedralsMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MechanicalMeasure" == simulationStepType and "StressMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MechanicalMeasure::StressMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::StressMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "ThermodynamicQuantityMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ThermodynamicMeasure::ThermodynamicQuantityMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::ThermodynamicQuantityMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "ThermodynamicIntegration" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ThermodynamicMeasure::ThermodynamicIntegration simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::ThermodynamicIntegration>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanRadius" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::MeanRadius simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::MeanRadius>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "GyrationRadius" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::GyrationRadius simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::GyrationRadius>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanSquareDisplacement" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::MeanSquareDisplacement simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::MeanSquareDisplacement>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "MeanAngularCorrelation" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::MeanAngularCorrelation simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::MeanAngularCorrelation>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "CenterOfMassPosition" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::CenterOfMassPosition simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::CenterOfMassPosition>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "Height" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::Height simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::Height>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("TopologicalMeasures" == simulationStepType and "PatchPolymers" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected TopologicalMeasures::PatchPolymers simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::PatchPolymers>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ParticlesListMeasure" == simulationStepType and "PotentialMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ParticlesListMeasure::PotentialMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::PotentialMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MagneticMeasure" == simulationStepType and "MeasureTotalMagnetization" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MagneticMeasure::MeasureTotalMagnetization simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::MeasureTotalMagnetization>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MagneticMeasure" == simulationStepType and "MeasureMeanMagnetization" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MagneticMeasure::MeasureMeanMagnetization simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::MeasureMeanMagnetization>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ExperimentMeasures" == simulationStepType and "AFMMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ExperimentMeasures::AFMMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::AFMMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("GeometricalMeasure" == simulationStepType and "DistanceBetweenCentersOfMass" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected GeometricalMeasure::DistanceBetweenCentersOfMass simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::DistanceBetweenCentersOfMass>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ParticlesListMeasure" == simulationStepType and "ContactsMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ParticlesListMeasure::ContactsMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::ContactsMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("ThermodynamicMeasure" == simulationStepType and "InteractorsListEnergyMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected ThermodynamicMeasure::InteractorsListEnergyMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::InteractorsListEnergyMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MechanicalMeasure" == simulationStepType and "HessianMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MechanicalMeasure::HessianMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::HessianMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MechanicalMeasure" == simulationStepType and "PairwiseForceMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MechanicalMeasure::PairwiseForceMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::PairwiseForceMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("MechanicalMeasure" == simulationStepType and "ForceBetweenSetsMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected MechanicalMeasure::ForceBetweenSetsMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::ForceBetweenSetsMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VQCMMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VQCMMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasureFromMobility" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VQCMMeasureFromMobility simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VQCMMeasureFromMobility>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VAFMMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VAFMMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VAFMMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VQCMMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VQCMMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VQCMMeasureFromMobility" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VQCMMeasureFromMobility simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VQCMMeasureFromMobility>(pg,integrator,ff,data,path.back());
            found = true;
        }
        if("OscillatingFluidMeasure" == simulationStepType and "VAFMMeasure" == simulationStepSubType){
            System::log<System::MESSAGE>("[SimulationStepLoader] (%s) Detected OscillatingFluidMeasure::VAFMMeasure simulationStep",path.back().c_str());
            simulationStep = std::make_shared<SimulationStep::SimulationMeasures::VAFMMeasure>(pg,integrator,ff,data,path.back());
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[SimulationStepLoader] (%s) Could not find simulationStep %s::%s",
                                            path.back().c_str(),simulationStepType.c_str(),simulationStepSubType.c_str());
        }

        return simulationStep;

    }

    }}}
#endif
