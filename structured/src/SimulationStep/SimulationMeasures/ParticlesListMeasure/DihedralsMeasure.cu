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

    class DihedralsMeasure: public SimulationStepBase{

            std::string   outputFilePath;
            std::ofstream outputFile;

            std::vector<int4> ids4;

        public:

            DihedralsMeasure(std::shared_ptr<ParticleGroup>             pg,
                             std::shared_ptr<IntegratorManager> integrator,
                             std::shared_ptr<ForceField>                ff,
                             DataEntry& data,
                             std::string name):SimulationStepBase(pg,integrator,ff,data,name){

                //Read parameters

                outputFilePath = data.getParameter<std::string>("outputFilePath");

                //////////////////////////////////////////////

                auto id_i = data.getData<int>("id_i");
                auto id_j = data.getData<int>("id_j");
                auto id_k = data.getData<int>("id_k");
                auto id_l = data.getData<int>("id_l");

                for(int index=0;index<data.getDataSize();index++){
                    ids4.push_back({id_i[index],id_j[index],id_k[index],id_l[index]});
                }
            }

            void init(cudaStream_t st) override{
                bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);
            }

            void applyStep(ullint step, cudaStream_t st) override{

                Box box = gd->getEnsemble()->getBox();

                auto pos = pd->getPos(access::location::cpu,
                                      access::mode::read);

                const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                for(auto ids : ids4){
                    int i = sortedIndex[ids.x];
                    int j = sortedIndex[ids.y];
                    int k = sortedIndex[ids.z];
                    int l = sortedIndex[ids.w];

                    real3 posi = make_real3(pos[i]);
                    real3 posj = make_real3(pos[j]);
                    real3 posk = make_real3(pos[k]);
                    real3 posl = make_real3(pos[l]);

                    const real3 dij = box.apply_pbc(posi - posj);
                    const real3 djk = box.apply_pbc(posj - posk);
                    const real3 dlk = box.apply_pbc(posl - posk);

                    const real3 aijk = cross(dij,djk);
                    const real3 ajkl = cross(dlk,djk);

                    const real raijk2=dot(aijk,aijk);
                    const real rajkl2=dot(ajkl,ajkl);

                    const real inv_raijkl = rsqrt(raijk2*rajkl2);

                    real cos_dih = dot(aijk,ajkl)*inv_raijkl;
                    cos_dih=min(real( 1.0),cos_dih);
                    cos_dih=max(real(-1.0),cos_dih);

                    const real rjk     = sqrt(dot(djk,djk));

                    real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
                    sin_dih=min(real( 1.0),sin_dih);
                    sin_dih=max(real(-1.0),sin_dih);

                    const real ang = atan2(sin_dih,cos_dih);

                    outputFile << ang << " ";

                }

                outputFile << std::endl;
            }
    };


}}}}

REGISTER_SIMULATION_STEP(
    ParticlesListMeasure,DihedralsMeasure,
    uammd::structured::SimulationStep::SimulationMeasures::DihedralsMeasure
)
