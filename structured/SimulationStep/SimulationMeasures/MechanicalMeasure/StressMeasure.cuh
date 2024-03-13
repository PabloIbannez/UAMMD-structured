#ifndef __STRESS_MEASURE__
#define __STRESS_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

    namespace StressMeasure_ns{

    __global__ void computeVolume(int numberParticles,
                                  uammd::structured::VerletConditionalListSet<conditions::all>::NeighbourListData nlData,
                                  Box    box,
                                  real4 *position,
                                  int*   id_ptr,
                                  real*  volume,
                                  int*   volumeId){

        int id = blockIdx.x*blockDim.x + threadIdx.x;

        if(id>=numberParticles){return;}

        const int nn = nlData.numberNeighbours[id];

        if(nn>int(1)){

            int const * nl = nlData.neighbourList + nlData.neighbourStart[id];
            const int i_global = nl[0];

            const real4 pos_i = position[i_global];

            real invr_sum  = 0.0;
            real invr2_sum = 0.0;
            for(int n = int(1); n<nn ; n += int(1)) {

                const int j_global = nl[n*numberParticles];

                const real4 pos_j = position[j_global];

                const real3 dr = box.apply_pbc(make_real3(pos_j - pos_i));

                const real invr = rsqrt(dot(dr,dr));

                invr_sum  += invr;
                invr2_sum += invr*invr;
            }

            const real radius = real(0.5)*invr_sum/invr2_sum;

            volume[id]   = real(4.0/3.0)*real(M_PI)*radius*radius*radius;
            volumeId[id] = id_ptr[i_global];

        } else {

            int const * nl = nlData.neighbourList + nlData.neighbourStart[id];
            const int i_global = nl[0];

            volume[id]   = real(0.0);
            volumeId[id] = id_ptr[i_global];
        }
    }

    }

class StressMeasure: public SimulationStepBase{

        std::string    outputFilePath;
        std::ofstream  outputFile;

        std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::all>> nl;
        std::shared_ptr<conditions::all>  cond;

        real radiusCutOff;

        thrust::device_vector<real> volume;
        thrust::device_vector<int>  volumeId;

        thrust::host_vector<real> volume_host;
        thrust::host_vector<int>  volumeId_host;

    public:

        StressMeasure(std::shared_ptr<ParticleGroup>             pg,
                      std::shared_ptr<IntegratorManager> integrator,
                      std::shared_ptr<ForceField>                ff,
                      DataEntry& data,
                      std::string name):SimulationStepBase(pg,
                                                           integrator,ff,
                                                           data,name){

            //Read parameters
            outputFilePath = data.getParameter<std::string>("outputFilePath");
            radiusCutOff   = data.getParameter<real>("radiusCutOff");

            cond = std::make_shared<conditions::all>(gd,pd,data);
            nl = std::make_shared<uammd::structured::VerletConditionalListSet<conditions::all>>(gd,pg,data,cond,name);

        }


        void init(cudaStream_t st) override{

            Backup::openFile(this->sys,outputFilePath,outputFile);

            //Init neighbor list
            nl->setCutOff(radiusCutOff);

            //Init volume
            volume.resize(pg->getNumberParticles());
            volumeId.resize(pg->getNumberParticles());
        }

        void applyStep(ullint step, cudaStream_t st) override{

            //Set stress to zero
				    {
              auto stress = this->pd->getStress(access::location::gpu, access::mode::write);
				    	thrust::fill(thrust::cuda::par.on(st),
				    			         stress.begin(),
				    			         stress.end(),
				    			         tensor3(0));
				    }
            //Sum stress
            {
                for(auto& interactor : this->topology->getInteractors()){

                    //Create computable
                    uammd::Interactor::Computables compTmp;
                    compTmp.stress = true;

                    interactor.second->sum(compTmp,st);
                }
            }
            //Compute volume
            {
                nl->update(st);

                auto nlData = nl->getNeighbourList("all");

                int numberParticles = nlData.N; //numberParticles in the neig list

                int Nthreads=512;
                int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                Box box = gd->getEnsemble()->getBox();

                real4* position_ptr = pd->getPos(access::location::gpu,access::mode::read).raw();
                int*   id_ptr       = pd->getId(access::location::gpu,access::mode::read).raw();

                real*  volume_ptr   = thrust::raw_pointer_cast(volume.data());
                int*   volumeId_ptr = thrust::raw_pointer_cast(volumeId.data());

                StressMeasure_ns::computeVolume<<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                                            nlData,
                                                                            box,
                                                                            position_ptr,
                                                                            id_ptr,
                                                                            volume_ptr,
                                                                            volumeId_ptr);

                cudaDeviceSynchronize();
                CudaCheckError();

                //Copy to host
                volume_host   = volume;
                volumeId_host = volumeId;
            }

            {
                auto id   = pd->getId(access::location::cpu,
                                      access::mode::read);

                auto pos   = pd->getPos(access::location::cpu,
                                        access::mode::read);

                auto rad   = pd->getRadius(access::location::cpu,
                                           access::mode::read);

                auto stress = pd->getStress(access::location::cpu,
                                            access::mode::readwrite);


                auto groupIndex  = pg->getIndexIterator(access::location::cpu);
                auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                Box box = gd->getEnsemble()->getBox();

                outputFile <<"#Lx="<<box.boxSize.x*0.5
                           <<";Ly="<<box.boxSize.y*0.5
                           <<";Lz="<<box.boxSize.z*0.5<<";"<< std::endl;
                std::map<int,int> id_index;
                fori(0,pg->getNumberParticles()){
                    int id_   = id[groupIndex[i]];
                    int index = sortedIndex[id_];

                    id_index[id_]=index;
                }

                std::map<int,int> id_index_volume;
                fori(0,pg->getNumberParticles()){
                    int id_   = id[groupIndex[i]];
                    forj(0,volumeId_host.size()){
                        if(id_ == volumeId_host[j]){
                            id_index_volume[id_]=j;
                            break;
                        }
                    }
                }

                for(const auto& ii : id_index){

                    int id_   = ii.first;

                    int index        = ii.second;
                    int index_volume = id_index_volume.at(id_);

                    real4 pc = pos[index];
                    real3 p  = box.apply_pbc(make_real3(pc));
                    int  type   = int(pc.w);
                    real radius = rad[index];

                    real vol = volume_host[index_volume];
                    tensor3 stress_ = stress[index];

                    outputFile << std::left
                               << std::setw(6)
                               << p      << " "
                               << radius << " "
                               << type   << " "
                               << vol    << " "
                               << stress_ << std::endl;
                    }
            }

        }
};

}}}}

#endif
