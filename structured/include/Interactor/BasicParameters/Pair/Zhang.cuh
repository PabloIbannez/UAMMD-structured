#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Pairs{

        struct Zhang {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real radius;
                real epsilon;
                real mu;
                real chi;
                real theta;
            };

            struct PairParameters{
                real rb;
                real rc;
                real epsilon;
                real stheta;
                real mu;
                real chi;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.radius  = real(info.at("radius"));
                param.epsilon = real(info.at("epsilon"));
                param.mu      = real(info.at("mu"));
                param.chi     = real(info.at("chi"));
                param.theta   = real(info.at("theta"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.rb       = in_par.radius*pow(2,7.0/6.0);
                param.rc       = in_par.radius*5.2;
                param.epsilon  = in_par.epsilon;
                param.stheta   = sin(in_par.theta);
                param.mu       = in_par.mu;
                param.chi      = in_par.chi;

                return param;
            }
        };

    }

}}}}


