#ifndef PARAMETERHANDLER_CUH
#define PARAMETERHANDLER_CUH

#include<thrust/device_vector.h>
namespace uammd{
namespace structured{
        
        template<class Type>
        class TypeParameterHandler{

            public:

            using TypeParameters      = typename Type::TypeParameters;          
            using InputTypeParameters = typename Type::InputTypeParameters;

            private:

            thrust::device_vector<TypeParameters> typeParameters;

            int ntypes;
            std::map<int,InputTypeParameters> idParameters;

            public:
            TypeParameterHandler():ntypes(1){}
            ~TypeParameterHandler(){}

            InputTypeParameters readTypeParameters(std::string line){
                return Type::readTypeParameters(line);
            }

            void add(int t,InputTypeParameters p){
                int new_ntypes = ntypes;
                if(t >= ntypes) new_ntypes = t+1;

                typeParameters.resize(new_ntypes);

                if(new_ntypes != ntypes){
                    auto tmp = typeParameters;
                    fori(0,ntypes){
                        typeParameters[i] = tmp[i];
                    }
                    ntypes = new_ntypes;
                }

                typeParameters[t]=Type::processTypeParameters(p);
                idParameters[t]=p;
            }
            
            int getNumTypes(){
                return ntypes;
            }

            std::vector<int> getTypeIdList(){
                std::vector<int> idList;
                idList.reserve(idParameters.size());

                for(auto const& e : idParameters){
                    idList.push_back(e.first);
                }

                return idList;
            }

            InputTypeParameters getTypeParameters(int id){
                return idParameters[id];
            }
            
            struct TypeIterator{
                TypeParameters * globalMem;
                
                int ntypes;
                TypeIterator(TypeParameters * globalMem, int ntypes): globalMem(globalMem), ntypes(ntypes){}
                
                size_t getSharedMemorySize(){
                    return 0;
                }

                inline __device__ void zero(){	  
                }

                inline __device__ TypeParameters operator()(int t) const{
                    //printf("ntypes:%i type:%i\n",ntypes,t);
                    if(ntypes==1) return this->globalMem[0];


                    #if CUB_PTX_ARCH < 300
                    constexpr auto cubModifier = cub::LOAD_DEFAULT;
                    #else
                    constexpr auto cubModifier = cub::LOAD_CA;
                    #endif

                    cub::CacheModifiedInputIterator<cubModifier, TypeParameters> itr(globalMem);

                    return itr[t];
                }

                };

                TypeIterator getTypeIterator(){
                    auto tp = thrust::raw_pointer_cast(typeParameters.data());
                    return TypeIterator(tp, ntypes);
                }
                
        };

        template<class PairType>
        class PairParameterHandler{

            public:
            
            using PairParameters       = typename PairType::PairParameters;          
            using InputPairParameters  = typename PairType::InputPairParameters;          

            private:

            std::set<std::pair<int,int>> addedPairs;

            thrust::device_vector<PairParameters>  pairParameters;

            int ntypes;

            public:

            PairParameterHandler():ntypes(1){}
            ~PairParameterHandler(){}
            
            InputPairParameters readPairParameters(std::string line){
                return PairType::readPairParameters(line);
            }

            void update(int ti, int tj,PairParameters newp){
                pairParameters[ti+ntypes*tj] = newp;
                if(ti != tj){
                    pairParameters[tj+ntypes*ti] = newp;
                }
            }

            void add(int ti, int tj,InputPairParameters p){
                int new_ntypes = ntypes;
                if(ti >= ntypes) new_ntypes = ti+1;
                if(tj >= ntypes) new_ntypes = tj+1;
                pairParameters.resize(new_ntypes*new_ntypes);

                if(new_ntypes != ntypes){
                    auto tmp = pairParameters;
                    fori(0,ntypes)
                        forj(0,ntypes){
                            pairParameters[i+new_ntypes*j] = tmp[i+ntypes*j];
                        }
                    ntypes = new_ntypes;
                }

                addedPairs.insert(std::make_pair(ti,tj));
                addedPairs.insert(std::make_pair(tj,ti));

                this->update(ti,tj,PairType::processPairParameters(p));

            }
            
            std::set<std::pair<int,int>> getAddedPairs(){
                return addedPairs;
            }
            
            bool isPairAdded(std::pair<int,int> pair){
                return bool(addedPairs.count(pair));
            }

            int getNumTypes(){
                return ntypes;
            }
            
            PairParameters getPairParameters(int ti,int tj){
                if(ti>tj) thrust::swap(ti,tj);
                int index = ti+this->ntypes*tj;	
                
                return pairParameters[index];
            }

            struct PairIterator{
                PairParameters * globalMem;
                int ntypes;
                PairIterator(PairParameters * globalMem, int ntypes): globalMem(globalMem), ntypes(ntypes){}

                size_t getSharedMemorySize(){
                    return 0;                }

                inline __device__ void zero(){	  
                }

                inline __device__ PairParameters operator()(int ti, int tj){
                    if(ntypes==1) return this->globalMem[0];


                    #if CUB_PTX_ARCH < 300
                    constexpr auto cubModifier = cub::LOAD_DEFAULT;
                    #else
                    constexpr auto cubModifier = cub::LOAD_CA;
                    #endif

                    cub::CacheModifiedInputIterator<cubModifier, PairParameters> itr(globalMem);

                    if(ti>tj) thrust::swap(ti,tj);

                    int typeIndex = ti+this->ntypes*tj;	
                    if(ti >= ntypes || tj >= ntypes) typeIndex = 0;

                    return itr[typeIndex];
                }

            };

            PairIterator getPairIterator(){
                auto tp = thrust::raw_pointer_cast(pairParameters.data());
                return PairIterator(tp, ntypes);
            }

        };
}}

#endif
