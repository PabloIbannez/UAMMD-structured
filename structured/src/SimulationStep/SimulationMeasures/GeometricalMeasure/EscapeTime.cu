#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{


class insideCheckerBase_{
    public:
        virtual bool isInsideSystem(real3 pos){
            return false;
        }
};

class insideSphereChecker{
    private:
        std::vector<real3> center;
        std::vector<real> radius2;
        std::string logicalSum;
    public:
        insideSphereChecker(DataEntry& data){
            //Read Data from input file.
            auto centerData = data.getData<real3>("center");
            auto radiusData = data.getData<real>("radius");
            logicalSum      = data.getParameter<std::string>("logicalSum","all");
            if(logicalSum!="all" and logicalSum!="any"){
                System::log<System::CRITICAL>("[EscapeSphereTime] logicalSum must be 'any' or 'all' ");
            }

            for (int i = 0; i < centerData.size(); i++){
                center.push_back({centerData[i].x,centerData[i].y,centerData[i].z});
                radius2.push_back({radiusData[i]*radiusData[i]});
            }
        }
    private:
        bool isInsideAll(real3 pos){ //Check if the particle is inside ALL the spheres, not any.
            for (int sphere = 0; sphere < center.size(); sphere++){
                real3 diff = pos - center[sphere];
                real diff2 = dot(diff,diff);

                if(diff2 > radius2[sphere]){
                    return false;
                }
            }
            return true;
        }

        bool isInsideAny(real3 pos){ //Check if the particle is inside ANY sphere, not all.
            for (int sphere = 0; sphere < center.size(); sphere++){
                real3 diff = pos - center[sphere];
                real diff2 = dot(diff,diff);

                if(diff2 < radius2[sphere]){
                    return true;
                }
            }
            return false;
        }

    public:
        bool isInsideSystem(real3 pos){
            if(logicalSum == "all"){
                return isInsideAll(pos);
            }

            if(logicalSum == "any"){
                return isInsideAny(pos);
            }
            System::log<System::CRITICAL>("[EscapeSphereTime] logicalSum must be 'any' or 'all' ");
            return false;
        }

};

class insidePolygonChecker{
    private:
        std::vector<real3> normalVector;
        std::vector<real3> independentVector;
    public:
        insidePolygonChecker(DataEntry& data){
            //Read Data from input file.
            auto normalVectorData      = data.getData<real3>("normalVector");
            auto independentVectorData = data.getData<real3>("independentVector");

            for (int i = 0; i < normalVectorData.size(); i++){
                normalVector.push_back({normalVectorData[i].x,normalVectorData[i].y,normalVectorData[i].z});
                independentVector.push_back({independentVectorData[i].x,independentVectorData[i].y,independentVectorData[i].z});
            }
        }

        bool isInsideSystem(real3 pos){ //Check if the particle is below all the planes contained in normalVector and independentVector. By computing the dot product of normalVector * (pos-independentVector)
            for (int plane = 0; plane < normalVector.size(); plane++){
                real3 particlePositionMinusIndependentVector = pos - independentVector[plane];
                real dotProduct = dot(normalVector[plane],particlePositionMinusIndependentVector);

                if(dotProduct < 0){
                    return false;
                }
            }
            return true;
        }
};

bool stringToBool(const std::string& str) {
    if (str == "true") {
        return true;
    } else if (str == "false") {
        return false;
    } else {
        System::log<System::CRITICAL>("[EscapeTime] stopSimulation must be 'true' or 'false'");
    }
}

template <class InsideSystemChecker>
class EscapeTime: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        bool stopSimulation;
        std::shared_ptr<InsideSystemChecker> checker;
        std::map<int,bool> isParticleInsideSystem;

        int particlesIn = 0;
        int particlesToEscape;

    public:

        EscapeTime(std::shared_ptr<ParticleGroup>         pg,
               std::shared_ptr<IntegratorManager> integrator,
               std::shared_ptr<ForceField>                ff,
               DataEntry& data,
               std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            //Read parameters from input file.
            outputFilePath = data.getParameter<std::string>("outputFilePath");

            checker = std::make_shared<InsideSystemChecker>(data);

            auto pos = pd->getPos(access::location::cpu, access::mode::read);

            stopSimulation  = stringToBool(data.getParameter<std::string>("stopSimulation","false"));
            particlesToEscape = data.getParameter<int>("particlesToEscape",pos.size());

            const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
            for(int i = 0; i < pos.size(); i++){
                int index = sortedIndex[i];
                real3 particlePosition = make_real3(pos[index]);
                isParticleInsideSystem[index] = checker->isInsideSystem(particlePosition);
                if (!isParticleInsideSystem[index]){
                    outputFile << 0 << " " << index << std::endl;
                    particlesIn++;
                }
            }
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);
        }

        void applyStep(ullint step, cudaStream_t st) override{
            auto pos = pd->getPos(access::location::cpu, access::mode::read);
            const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
            for(int i = 0; i < pos.size(); i++){
                int index = sortedIndex[i];
                real3 particlePosition = make_real3(pos[index]);
                bool isInsideSystem = checker->isInsideSystem(particlePosition);

                if(isParticleInsideSystem[index] && !isInsideSystem){
                    outputFile << step << " " << index << std::endl;
                    isParticleInsideSystem[index] = false;
                    particlesIn++;
                    if(stopSimulation and (particlesIn >= particlesToEscape)){
                        //System::log<System::CRITICAL>("[EscapeTime] Simulation interrumpted, all particles have escaped");
                        System::log<System::MESSAGE>("[EscapeTime] All particlesToEscape have escaped!");
                        this->sys->setState(ExtendedSystem::SIMULATION_STATE::STOPPED);
                    }
                }
            }
        }
};

}}}}

using insideSphereChecker = uammd::structured::SimulationStep::SimulationMeasures::insideSphereChecker;
using insidePolygonChecker = uammd::structured::SimulationStep::SimulationMeasures::insidePolygonChecker;

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,EscapeSphereTime,
    uammd::structured::SimulationStep::SimulationMeasures::EscapeTime<insideSphereChecker>
)

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,EscapePolygonTime,
    uammd::structured::SimulationStep::SimulationMeasures::EscapeTime<insidePolygonChecker>
)
