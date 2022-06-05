#ifndef __COMMON_PARAMETERS__
#define __COMMON_PARAMETERS__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace CommonParameters{
    
    namespace LennardJones{
    
        template<class Topology>
        class LennardJones: public ParameterUpdatable{
    
            public:

                struct InteractionParameters{
                
                    struct InputPairParameters{
                        real epsilon;
                        real sigma;
                    };
                
                    struct PairParameters{
                        real epsilon;
                        real sigma;
                    };
                
                    static inline __host__ InputPairParameters readPairParameters(std::string& line){
                    
                        std::stringstream ss;
                        
                        InputPairParameters param;

                        ss.str(line);

                        ss  >> param.epsilon >> param.sigma;
                        
                        return param;
                    
                    }
                
                    static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){
                
                        PairParameters params;
                
                        params.epsilon = in_par.epsilon;
                        params.sigma   = in_par.sigma;
                
                        return params;
                    }
                };
                
                using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;

            private:

                bool parametersRequested = false;
                
                std::shared_ptr<System>        sys;
                std::shared_ptr<ParticleData>  pd;
                std::shared_ptr<ParticleGroup> pg;
                std::shared_ptr<Topology>      top;

                std::shared_ptr<ParameterPairsHandler> interParam;
                
                std::string label;

                bool checkSameNumber=true;

            public:

                struct Parameters{
                    std::string label;

                    bool checkSameNumber=true;
                };

                LennardJones(std::shared_ptr<System>        sys,
                             std::shared_ptr<ParticleData>  pd,
                             std::shared_ptr<ParticleGroup> pg,
                             std::shared_ptr<Topology>      top,
                             Parameters par):sys(sys),
                                             pd(pd),
                                             pg(pg),
                                             top(top),
                                             label(par.label),
                                             checkSameNumber(par.checkSameNumber){
                    
                    interParam = this->top->template readPairs<InteractionParameters>(label);
                        
                    auto typesParam    = this->top->getTypes();
                    auto typeList = typesParam->getTypeIdList();
                        
                    if((typesParam->getNumTypes() != interParam->getNumTypes()) and checkSameNumber){
                        this->sys->template log<System::CRITICAL>("[LennardJonesParameters] The number of types of the types handler (%i) "
                                                                  "does not match with the number of types of the pair handler (%i), pair label: %s",
                                                                  typesParam->getNumTypes(),interParam->getNumTypes(),label.c_str());
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
                        
                            this->sys->template log<System::DEBUG1>("[LennardJonesParameters] "
                                                                     "(%s %s) Pair param, eps:%f and sigma:%f",
                                                                      type_i_name.c_str(),type_j_name.c_str(),
                                                                      interParam->getPairParameters(type_i_id,type_j_id).epsilon,
                                                                      interParam->getPairParameters(type_i_id,type_j_id).sigma);
                        }}
                        
                        parametersRequested=true;
                    }

                    return interParam;
                };

        };
    }

    namespace StatisticalPotential{
    
        template<class Topology>
        class StatisticalPotential: public ParameterUpdatable{

            public:

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

                using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;

            private:

                bool parametersRequested = false;
                
                std::shared_ptr<System>        sys;
                std::shared_ptr<ParticleData>  pd;
                std::shared_ptr<ParticleGroup> pg;
                std::shared_ptr<Topology>      top;

                std::shared_ptr<ParameterPairsHandler> interParam;
                    
                std::string label;

                real refTemperature;
                real epsilon_0;
                real lambda;

                bool checkSameNumber;

                real kBT;

            public:

                struct Parameters{
                    std::string label;

                    real refTemperature;
                    real epsilon_0;
                    real lambda;

                    bool checkSameNumber=true;
                };

                StatisticalPotential(std::shared_ptr<System>        sys,
                                     std::shared_ptr<ParticleData>  pd,
                                     std::shared_ptr<ParticleGroup> pg,
                                     std::shared_ptr<Topology>      top,
                                     Parameters par):sys(sys),
                                                     pd(pd),
                                                     pg(pg),
                                                     top(top),
                                                     label(par.label),
                                                     refTemperature(par.refTemperature),
                                                     epsilon_0(par.epsilon_0),
                                                     lambda(par.lambda),
                                                     checkSameNumber(par.checkSameNumber),
                                                     kBT(Topology::Units::KBOLTZ*par.refTemperature){

                    interParam = this->top->template readPairs<InteractionParameters>(label);
                        
                    auto typesParam    = this->top->getTypes();
                    auto typeList = typesParam->getTypeIdList();
                        
                    //Scaling

                    if((typesParam->getNumTypes() != interParam->getNumTypes()) and checkSameNumber){
                        this->sys->template log<System::CRITICAL>("[StatisticalPotential] The number of types of the types handler (%i) "
                                                                  "does not match with the number of types of the pair handler (%i), pair label: %s",
                                                                  typesParam->getNumTypes(),interParam->getNumTypes(),label.c_str());

                    } else {

                        auto addedPairs = interParam->getAddedPairs();

                        for(int type_i=0;type_i<typeList.size();type_i++){
                            for(int type_j=type_i;type_j<typeList.size();type_j++){

                                int type_i_id = typeList[type_i];
                                int type_j_id = typeList[type_j];

                                if(addedPairs.count(std::make_pair(type_i_id,type_j_id))!=0){
                                    real epsilon       = interParam->getPairParameters(type_i_id,type_j_id).epsilon;
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
