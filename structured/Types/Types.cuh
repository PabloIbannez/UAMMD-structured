#ifndef __TYPES__
#define __TYPES__

namespace uammd{
namespace structured{
namespace Types{
    
    struct BASIC{

        //Parameter handler structs
        struct InputTypeParameters{

            std::string name;

            real mass=-1;
            real radius=-1;
            real charge=0;
        };

        struct TypeParameters{

            real mass;
            real radius;
            real charge;
        };
        
        static inline __host__ InputTypeParameters readTypeParameters(std::string& line){
            
            std::stringstream ss;
            
            InputTypeParameters param;

            ss.str(line);

            ss >> param.name          >>
                  param.mass          >>
                  param.radius        >>
                  param.charge        ;

            return param;

        }

        static inline __host__ TypeParameters processTypeParameters(InputTypeParameters in_par){
            TypeParameters param;

            param.mass             = in_par.mass;
            param.radius           = in_par.radius;
            param.charge           = in_par.charge;

            return param;

        }

        static void load(std::shared_ptr<ParticleData> pd,const int& index,const InputTypeParameters& tpParam){

            auto mass   = pd->getMass(access::location::cpu, access::mode::write);
            auto radius = pd->getRadius(access::location::cpu, access::mode::write);
            auto charge = pd->getCharge(access::location::cpu, access::mode::write);

            TypeParameters typeParam = processTypeParameters(tpParam);

            mass[index]   = typeParam.mass;
            radius[index] = typeParam.radius;
            charge[index] = typeParam.charge;

            /*
            std::cout << index         << " " << 
                         mass[index]   << " " <<
                         radius[index] << " " <<
                         charge[index] << std::endl;*/
                    
        }
    };
    
    struct SASA{

        //Parameter handler structs
        struct InputTypeParameters{

            std::string name;

            real mass=-1;
            real radius=-1;
            real charge=0;
            real SASArc=-1;
        };

        struct TypeParameters{

            real mass;
            real radius;
            real charge;
            real SASArc;
        };
        
        static inline __host__ InputTypeParameters readTypeParameters(std::string& line){
            
            std::stringstream ss;
            
            InputTypeParameters param;

            ss.str(line);

            ss >> param.name          >>
                  param.mass          >>
                  param.radius        >>
                  param.charge        >>
                  param.SASArc        ;
            
            return param;

        }

        static inline __host__ TypeParameters processTypeParameters(InputTypeParameters in_par){
            TypeParameters param;

            param.mass     = in_par.mass;
            param.radius   = in_par.radius;
            param.charge   = in_par.charge;
            param.SASArc   = in_par.SASArc;

            return param;

        }

        static void load(std::shared_ptr<ParticleData> pd,const int& index,const InputTypeParameters& tpParam){

            auto mass   = pd->getMass(access::location::cpu, access::mode::write);
            auto radius = pd->getRadius(access::location::cpu, access::mode::write);
            auto charge = pd->getCharge(access::location::cpu, access::mode::write);

            TypeParameters typeParam = processTypeParameters(tpParam);

            mass[index]   = typeParam.mass;
            radius[index] = typeParam.radius;
            charge[index] = typeParam.charge;

            /*
            std::cout << index         << " " << 
                         mass[index]   << " " <<
                         radius[index] << " " <<
                         charge[index] << std::endl;*/
                    
        }
    };
    
    struct SURF{

        //Parameter handler structs
        struct InputTypeParameters{

            std::string name;

            real mass=-1;
            real radius=-1;
            real charge=0;
            real surf=0;
        };

        struct TypeParameters{

            real mass;
            real radius;
            real charge;
            real surf;
        };
        
        static inline __host__ InputTypeParameters readTypeParameters(std::string& line){
            
            std::stringstream ss;
            
            InputTypeParameters param;

            ss.str(line);

            ss >> param.name        >>
                  param.mass        >>
                  param.radius      >>
                  param.charge      >>
                  param.surf        ;
            
            return param;

        }

        static inline __host__ TypeParameters processTypeParameters(InputTypeParameters in_par){
            TypeParameters param;

            param.mass     = in_par.mass;
            param.radius   = in_par.radius;
            param.charge   = in_par.charge;
            param.surf   = in_par.surf;

            return param;

        }

        static void load(std::shared_ptr<ParticleData> pd,const int& index,const InputTypeParameters& tpParam){

            auto mass   = pd->getMass(access::location::cpu, access::mode::write);
            auto radius = pd->getRadius(access::location::cpu, access::mode::write);
            auto charge = pd->getCharge(access::location::cpu, access::mode::write);
            auto surf   = pd->getSurface(access::location::cpu, access::mode::write);

            TypeParameters typeParam = processTypeParameters(tpParam);

            mass[index]   = typeParam.mass;
            radius[index] = typeParam.radius;
            charge[index] = typeParam.charge;
            surf[index]   = typeParam.surf;

            /*
            std::cout << index         << " " << 
                         mass[index]   << " " <<
                         radius[index] << " " <<
                         charge[index] << " " <<
                         surf[index]   << std::endl;*/
                    
        }
    };

}}}

#endif
