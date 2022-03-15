#ifndef __COMMON_PARAMETERS__
#define __COMMON_PARAMETERS__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace CommonParameters{

    namespace StatisticalPotential{
    
        struct InteractionParameters{
        
            struct InputPairParameters{
                real epsilon;
            };
        
            struct PairParameters{
                real epsilon;
            };
        
            static inline __host__ InputPairParameters readPairParameters(const std::string& line){
            
                std::stringstream ss;
                
                InputPairParameters param;

                ss.str(line);

                ss >> param.epsilon;
                
                return param;
            
            }
        
            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){
        
                PairParameters params;
        
                params.epsilon = in_par.epsilon;
        
                return params;
            }
        };

        template<class Topology>
        class StatisticalPotential: public ParameterUpdatable{

            public:

                using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;

            private:

                bool parametersRequested = false;
                
                std::shared_ptr<System>        sys;
                std::shared_ptr<ParticleData>  pd;
                std::shared_ptr<ParticleGroup> pg;
                std::shared_ptr<Topology>      top;

                std::shared_ptr<ParameterPairsHandler> interParam;

                real kBT;

            public:

                StatisticalPotential(std::shared_ptr<System>        sys,
                                     std::shared_ptr<ParticleData>  pd,
                                     std::shared_ptr<ParticleGroup> pg,
                                     std::shared_ptr<Topology>      top,
                                     real refTemperature):sys(sys),
                                                          pd(pd),
                                                          pg(pg),
                                                          top(top),
                                                          kBT(Topology::Units::KBOLTZ*refTemperature){

                    interParam = std::make_shared<ParameterPairsHandler>();
                        
                    auto typesParam    = this->top->getTypes();
                    auto typeList = typesParam->getTypeIdList();
                        
                    for(int type_i=0;type_i<typeList.size();type_i++){
                        for(int type_j=type_i;type_j<typeList.size();type_j++){
                            
                            int type_i_id = typeList[type_i];
                            int type_j_id = typeList[type_j];
                            
                            typename InteractionParameters::InputPairParameters initParam;
                            initParam.epsilon = real(0.0);

                            interParam->add(type_i_id,type_j_id,initParam);
                    }}
                }

                void loadParameters(std::string label,real epsilon_0,real lambda,bool checkSameNumber = false){
                
                    auto typesParam    = this->top->getTypes();
                    auto tmpInterParam = this->top->template readPairs<InteractionParameters>(label);

                    if((typesParam->getNumTypes() != tmpInterParam->getNumTypes()) and checkSameNumber){
                        this->sys->template log<System::CRITICAL>("[StatisticalPotential] The number of types of the types handler (%i) "
                                                                  "does not match with the number of types of the pair handler (%i), pair label: %s",
                                                                  typesParam->getNumTypes(),tmpInterParam->getNumTypes(),label.c_str());

                    } else {
                        //Scale according to model type
                        auto typeList = typesParam->getTypeIdList();

                        auto addedPairsTmp = tmpInterParam->getAddedPairs();

                        for(int type_i=0;type_i<typeList.size();type_i++){
                            for(int type_j=type_i;type_j<typeList.size();type_j++){

                                int type_i_id = typeList[type_i];
                                int type_j_id = typeList[type_j];

                                if(addedPairsTmp.count(std::make_pair(type_i_id,type_j_id))!=0){
                                    real epsilon       = tmpInterParam->getPairParameters(type_i_id,type_j_id).epsilon;
                                    real epsilonScaled = lambda*(epsilon-epsilon_0)*kBT;

                                    typename InteractionParameters::InputPairParameters scaledParam;
                                    scaledParam.epsilon = epsilonScaled;

                                    interParam->add(type_i_id,type_j_id,scaledParam);

                                    std::string type_i_name = typesParam->getTypeParameters(type_i_id).name;
                                    std::string type_j_name = typesParam->getTypeParameters(type_j_id).name;

                                    this->sys->template log<System::DEBUG>("[StatisticalPotential] "
                                                                            "(%s), (%s %s) Pair param: %f, scaled to %f",
                                                                            label.c_str(),
                                                                            type_i_name.c_str(),type_j_name.c_str(),
                                                                            epsilon,interParam->getPairParameters(type_i_id,type_j_id).epsilon);
                                }

                            }
                        }
                    }
                }

                std::shared_ptr<ParameterPairsHandler> getParameters(){
                    if(!parametersRequested){
                    
                        auto typesParam    = this->top->getTypes();
                        auto typeList = typesParam->getTypeIdList();
                        
                        for(int type_i=0;type_i<typeList.size();type_i++){
                            for(int type_j=type_i;type_j<typeList.size();type_j++){
                                
                            int type_i_id = typeList[type_i];
                            int type_j_id = typeList[type_j];
                                    
                            std::string type_i_name = typesParam->getTypeParameters(type_i_id).name;
                            std::string type_j_name = typesParam->getTypeParameters(type_j_id).name;
                        
                            this->sys->template log<System::MESSAGE>("[StatisticalPotential] "
                                                                     "(%s %s) Pair param: %f",
                                                                      type_i_name.c_str(),type_j_name.c_str(),
                                                                      interParam->getPairParameters(type_i_id,type_j_id).epsilon);
                        }}
                        
                        parametersRequested=true;
                    }

                    return interParam;
                };
        };
    }
    

}}}}

#endif
