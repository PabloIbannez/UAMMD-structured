#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Pairs{

        struct DistanceSwitchCosine {

            static constexpr bool symmetric = true;

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real E;
                real K;
                real rc;
            };

            struct PairParameters{
                real E;
                real K;
                real rc;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.E       = real(info.at("E"));
                param.K       = real(info.at("K"));
                param.rc      = real(info.at("rc"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.E        = in_par.E;
                param.K        = in_par.K;
                param.rc       = in_par.rc;

                return param;
            }
        };

    }

}}}}
