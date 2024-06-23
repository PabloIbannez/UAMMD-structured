#ifndef __SIMULATION_STEP_PARTICLE_LIST_ANGLES_MEASURES__
#define __SIMULATION_STEP_PARTICLE_LIST_ANGLES_MEASURES__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

    class AnglesMeasure: public SimulationStepBase{

            std::string   outputFilePath;
            std::ofstream outputFile;

            std::vector<int3> ids3;

        public:

            AnglesMeasure(std::shared_ptr<ParticleGroup>             pg,
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

                for(int index=0;index<data.getDataSize();index++){
                    ids3.push_back({id_i[index],id_j[index],id_k[index]});
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

                for(auto ids : ids3){
                    int i = sortedIndex[ids.x];
                    int j = sortedIndex[ids.y];
                    int k = sortedIndex[ids.z];

                    real3 posi = make_real3(pos[i]);
                    real3 posj = make_real3(pos[j]);
                    real3 posk = make_real3(pos[k]);

                    //         i -------- j -------- k
                    //             <- rji     rjk ->
                    //Compute distances and vectors
                    //---rji---
                    const real3 rji = box.apply_pbc(posi - posj);
                    const real rji2 = dot(rji, rji);
                    //---rkj---
                    const real3 rjk = box.apply_pbc(posk - posj);
                    const real rjk2 = dot(rjk, rjk);

                    const real inv_rjirjk = rsqrt(rji2*rjk2);

                    real cijk = dot(rji, rjk)*inv_rjirjk; //cijk = cos (theta) = rji*rkj / mod(rji)*mod(rkj)
                    //Cos must stay in range
                    cijk = min(real( 1.0),cijk);
                    cijk = max(real(-1.0),cijk);

                    const real ang = acos(cijk);

                    outputFile << ang << " ";

                }

                outputFile << std::endl;
            }
    };

}}}}

#endif
