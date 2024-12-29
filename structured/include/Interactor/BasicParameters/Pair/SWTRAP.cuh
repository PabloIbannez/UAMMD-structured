#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Pairs{

        struct SWTRAP{

            static constexpr bool symmetric = false;

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real  E;
                real  rc;
                real4 R;
                real  Kswt;
                real  Krap;
            };

            struct PairParameters{
                real  E;
                real  rc;
                real4 R;
                real  Kswt;
                real  Krap;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.E    = real(info.at("E"));
                param.rc   = real(info.at("rc"));
                param.R    = real4(info.at("R"));
                param.Kswt = real(info.at("Kswt"));
                param.Krap = real(info.at("Krap"));

                real norm = sqrt(dot(param.R, param.R));
                if (abs(norm - 1.0) > 1e-6){
                    System::log<System::CRITICAL>("[SWTRAP Parameters] The R (%f, %f, %f, %f) given for the pair %s-%s is not a unit quaternion. The norm is %f.\n",
                                                  param.R.x, param.R.y, param.R.z, param.R.w, param.name_i.c_str(), param.name_j.c_str(), norm);
                }

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.E  = in_par.E;
                param.rc = in_par.rc;
                param.R  = in_par.R;
                param.Kswt = in_par.Kswt;
                param.Krap = in_par.Krap;

                return param;
            }
        };
    }

}}}}
