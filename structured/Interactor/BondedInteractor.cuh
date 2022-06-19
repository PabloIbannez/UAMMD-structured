#ifndef __BONDED_INTERACTOR__
#define __BONDED_INTERACTOR__

#include"Interactor/Interactor.cuh"

namespace uammd{
namespace structured{
namespace Interactor{
    namespace BondedInteractor_ns{

        template<typename Topology,class BondType>
        class BondReaderFromFile{

          public:

              using InputType = Topology;
          
          private:

            using Bond = typename BondType::Bond;
            
            std::shared_ptr<System> sys;
            std::shared_ptr<InputType> top;

            typename InputType::fileBlockIterator rawBondItr;

            std::string bondName;
            
            int nbonds;
        
          public:
        
            struct Parameters{
                std::string bondName;
            };
        
            BondReaderFromFile(std::shared_ptr<System> sys,
                               std::shared_ptr<InputType> top,
                               Parameters param):sys(sys),
                                                 top(top),
                                                 rawBondItr(top->getFileBlockIterator(param.bondName)),
                                                 nbonds(0){
                
                std::string line;
                while(rawBondItr.next(line)){
                    nbonds++;
                }
                rawBondItr.reset();
            }
        
            int getNumberBonds() {
                return nbonds;
            }
        
            Bond readNextBond(){
                
                std::string line;
                if(!rawBondItr.next(line)){
                    sys->log<System::CRITICAL>("There are not more bonds"
                                                "in the current file block");
                }
        
                return BondType::readBond(sys,line);
        
            }
        
        };
        
        template<class BondType>
        class BondReaderFromVector{
          
          public:

              using InputType = typename std::vector<typename BondType::Bond>;
          
          private:

            using Bond = typename BondType::Bond;
            
            std::shared_ptr<System> sys;

            std::shared_ptr<InputType> bondVector;

            std::string bondName;
            
            int nbonds;

            int bondIndex;
        
          public:
        
            struct Parameters{
                std::string bondName;
            };
        
            BondReaderFromVector(std::shared_ptr<System>    sys,
                                 std::shared_ptr<InputType> bondVector,
                                 Parameters param):sys(sys),
                                                   bondVector(bondVector),
                                                   nbonds(0),
                                                   bondIndex(-1){
                nbonds = bondVector->size();
            }
        
            int getNumberBonds() {
                return nbonds;
            }
        
            Bond readNextBond(){

                bondIndex++;
                
                if(bondIndex >= nbonds){
                    sys->log<System::CRITICAL>("There are not more bonds"
                                                "in the current bond vector");
                }
        
                return (*bondVector)[bondIndex];
        
            }
        
        };
        
        template<class BondType>
        class BondProcessor{
            
            std::shared_ptr<System> sys;
            std::shared_ptr<ParticleGroup> pg;

            using Bond = typename BondType::Bond;

            std::vector<Bond> bondList;
            
            std::vector<std::vector<int>> bondsPerParticle;
          
            public:
            
            BondProcessor(std::shared_ptr<ParticleGroup> pg):sys(pg->getParticleData()->getSystem()),
                                                             pg(pg){
                bondsPerParticle.resize(pg->getNumberParticles());
            }
        
            void registerBond(Bond b) {
        
                //TODO check if all bonds particle are in group
                int bondIndex = bondList.size();
                bondList.push_back(b);

                BondType::registerBond(sys,
                                       bondsPerParticle,
                                       bondIndex,
                                       b);
            }
        
            std::vector<int> getParticlesWithBonds() const {
            
                std::vector<int> particlesWithBonds;

                fori(0,bondsPerParticle.size()){
                    if(bondsPerParticle[i].size() > 0){
                        particlesWithBonds.push_back(i);
                    }
                }
                
                std::sort(particlesWithBonds.begin(),particlesWithBonds.end());

                return particlesWithBonds;
            }

            std::vector<Bond> getBondListOfParticle(int index) const {
                std::vector<Bond> bondListOfParticle;
                bondListOfParticle.resize(bondsPerParticle[index].size());
                fori(0, bondListOfParticle.size()) {
                    bondListOfParticle[i] = bondList[bondsPerParticle[index][i]];
                }
                return bondListOfParticle;
            }
        };

    template<class BondType_      ,
             class BondProcessor_ ,
             class BondReader_    >
        class BondList{
            
            private:

                std::shared_ptr<System>       sys;
                std::shared_ptr<ParticleGroup> pg;
                
                using BondType    = BondType_;
                using Bond        = typename BondType::Bond;
                
                using BondProcessor = BondProcessor_;
                using BondReader    = BondReader_;
                
                using InputType     = typename BondReader::InputType;

                std::shared_ptr<InputType> input;
                std::string bondName;
                
                std::shared_ptr<BondProcessor> bondProcessor;
            
                int nbonds;
                thrust::device_vector<Bond> bondList; 
                thrust::device_vector<int>  bondStart, bondEnd; 
                thrust::device_vector<int>  particlesWithBonds;

                void initBondProcessor(){

                    //Load top to bondProcessor
                    //----Init bondReader
                    typename BondReader::Parameters parBR;

                    parBR.bondName = this->bondName;

                    BondReader bondReader(this->sys,
                                          this->input,
                                          parBR);

                    this->nbonds = bondReader.getNumberBonds();

                    //----Init bondProcessor

                    this->bondProcessor = std::make_shared<BondProcessor>(sys,pg);

                    //----Register bonds in bonds processor
                    Bond bond;
                    for(int b=0;b<this->nbonds;b++){

                        bond = bondReader.readNextBond();
                        this->bondProcessor->registerBond(bond);
                    }

                }

                void initBondList(){

                    auto h_particlesWithBonds = bondProcessor->getParticlesWithBonds();

                    sys->log<System::DEBUG>("[BondList] Generating list");

                    sys->log<System::MESSAGE>("[BondList] Detected: %i bonds", 
                            this->nbonds);

                    sys->log<System::MESSAGE>("[BondList] %d particles"
                            " are involved in at least one bond.", 
                            h_particlesWithBonds.size());

                    const int NparticleswithBonds = h_particlesWithBonds.size();

                    std::vector<int> h_bondStart, h_bondEnd;

                    h_bondStart.resize(NparticleswithBonds, 0xffFFffFF);
                    h_bondEnd.resize(NparticleswithBonds, 0);

                    thrust::host_vector<Bond> bondListCPU(BondType::nPart*this->nbonds);

                    std::vector<Bond> blst;
                    fori(0, NparticleswithBonds) {

                        const int index = h_particlesWithBonds[i];
                        blst = bondProcessor->getBondListOfParticle(index);

                        const int nbondsi = blst.size();

                        int offset;
                        if(i>0){offset = h_bondEnd[i-1];}
                        else   {offset = 0;}

                        forj(0,nbondsi) {
                            bondListCPU[offset+j] = blst[j];
                        }

                        h_bondEnd[i]   = offset + nbondsi;
                        h_bondStart[i] = offset;
                    }

                    sys->log<System::DEBUG>("[BondList] Uploading list");

                    this->bondList = bondListCPU;

                    this->bondStart = h_bondStart;
                    this->bondEnd   = h_bondEnd;

                    this->particlesWithBonds = h_particlesWithBonds;
                }

            public:
                
                struct bondListInfo{
                    int   N;
                    int*  particlesWithBonds_ptr;
                    Bond* bondList_ptr;    
                    int*  bondStart_ptr;
                    int*  bondEnd_ptr; 
                };

                BondList(std::shared_ptr<ParticleGroup> pg,
                         std::shared_ptr<InputType>  input,
                         std::string bondName):sys(pg->getParticleData()->getSystem()),pg(pg),
                                               input(input),
                                               bondName(bondName){
                    this->initBondProcessor();
                    this->initBondList();
                }

                bondListInfo getBondListInfo(){
                    bondListInfo bndLst;
                
                    bndLst.N                      = particlesWithBonds.size();
                    bndLst.particlesWithBonds_ptr = thrust::raw_pointer_cast(particlesWithBonds.data());
                    bndLst.bondStart_ptr          = thrust::raw_pointer_cast(bondStart.data());
                    bndLst.bondEnd_ptr            = thrust::raw_pointer_cast(bondEnd.data());
                    bndLst.bondList_ptr           = thrust::raw_pointer_cast(bondList.data());

                    return bndLst;
                }
        };

        template<typename Bond,class BondedPotTransverser>
        __global__ void transverseBondListThreadPerParticle(const int   nParticlesWithBonds,
                                                            const int*  __restrict__ particlesWithBonds,
                                                            const int*  __restrict__ bondStart,
                                                            const int*  __restrict__ bondEnd,
                                                            const Bond* __restrict__ bondList,
                                                            const int * __restrict__ id2index,
                                                            BondedPotTransverser bPotTransverser) {

            int index = blockIdx.x*blockDim.x + threadIdx.x;

            if(index>=nParticlesWithBonds){return;}

            const int id_i = particlesWithBonds[index];
            
            const int index_global_i = id2index[id_i];
            
            const int first = bondStart[index];
            const int last  = bondEnd[index];
            
            typename BondedPotTransverser::resultType quantityLocal = bPotTransverser.zero();

            for(int b = first; b < last; b++) {

                Bond bond = bondList[b];

                bPotTransverser.accumulate(quantityLocal,bPotTransverser.compute(index_global_i,bond));
            }
            
            bPotTransverser.set(index_global_i,quantityLocal);
        }
    }

    template<class BondType_      ,
             class BondProcessor_ ,
             class BondReader_    >
    class BondedInteractor: public Interactor,
                            public ParameterUpdatableDelegate<BondType_> {
    
        private:
    
            using BondType    = BondType_;
            using Bond        = typename BondType::Bond;
            
            using BondProcessor = BondProcessor_;
            using BondReader    = BondReader_;
            
            using InputType     = typename BondReader::InputType;

            using BondListType  = BondedInteractor_ns::BondList<BondType,
                                                                BondProcessor,
                                                                BondReader>;

            using ForceTransverser  = typename BondType::ForceTransverser;
            using EnergyTransverser = typename BondType::EnergyTransverser;
            using VirialTransverser = typename BondType::VirialTransverser;
            
            std::shared_ptr<InputType> input;
            std::string bondName;
            
            std::shared_ptr<BondType> bondType;
    
            std::shared_ptr<BondListType> bondList;

        public:
    
    
            struct Parameters {
                std::string bondName;
            };
    
            BondedInteractor(std::shared_ptr<ParticleGroup> pg, 
                             std::shared_ptr<InputType> input,
                             std::shared_ptr<BondType> bondType,
                             Parameters par):Interactor(pg, "BondedInteractor " + par.bondName),
                                             input(input),
                                             bondType(bondType),
                                             bondName(par.bondName){

                this->setDelegate(bondType);

                bondList = std::make_shared<BondListType>(pg,input,bondName);
            }
    
            void sum(Computables comp,cudaStream_t st) override {
                
                auto bondListInfo = bondList->getBondListInfo();

                int N = bondListInfo.N;

                auto particlesWithBonds_ptr = bondListInfo.particlesWithBonds_ptr;
                auto bondStart_ptr          = bondListInfo.bondStart_ptr;
                auto bondEnd_ptr            = bondListInfo.bondEnd_ptr;
                auto bondList_ptr           = bondListInfo.bondList_ptr;

                auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                int Nthreads=128;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                if(comp.force == true){

                    ForceTransverser ft = bondType->getForceTransverser();
               
                    BondedInteractor_ns::transverseBondListThreadPerParticle<
                                         typename BondType::Bond,
                                         ForceTransverser>
                    <<<Nblocks,Nthreads,0,st>>>(N,
                                                particlesWithBonds_ptr,
                                                bondStart_ptr,
                                                bondEnd_ptr,
                                                bondList_ptr,
                                                id2index,
                                                ft);
                    CudaCheckError();
                }

                if(comp.energy == true){

                    EnergyTransverser et = bondType->getEnergyTransverser();
               
                    BondedInteractor_ns::transverseBondListThreadPerParticle<
                                         typename BondType::Bond,
                                         EnergyTransverser>
                    <<<Nblocks,Nthreads,0,st>>>(N,
                                                particlesWithBonds_ptr,
                                                bondStart_ptr,
                                                bondEnd_ptr,
                                                bondList_ptr,
                                                id2index,
                                                et);
                    CudaCheckError();
                }
                
                if(comp.virial == true){

                    VirialTransverser vt = bondType->getVirialTransverser();
               
                    BondedInteractor_ns::transverseBondListThreadPerParticle<
                                         typename BondType::Bond,
                                         VirialTransverser>
                    <<<Nblocks,Nthreads,0,st>>>(N,
                                                particlesWithBonds_ptr,
                                                bondStart_ptr,
                                                bondEnd_ptr,
                                                bondList_ptr,
                                                id2index,
                                                vt);
                    CudaCheckError();
                }
            }

            void updateBox(Box box){
                ParameterUpdatableDelegate<BondType_>::updateBox(box);
            }
    };

}}}

#endif
