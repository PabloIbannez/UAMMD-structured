#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Pairs{

        template<typename potential>
        struct Helix {

            static inline __host__ tensor3 loadR_H(real3 e_x,
                                                   real3 e_y,
                                                   real3 e_z){

                // Check e_x, e_y, e_z are orthonormal
                real e_x_norm = sqrt(dot(e_x,e_x));
                real e_y_norm = sqrt(dot(e_y,e_y));
                real e_z_norm = sqrt(dot(e_z,e_z));

                if (abs(e_x_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x is not a unit vector");
                }
                if (abs(e_y_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_y is not a unit vector");
                }
                if (abs(e_z_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_z is not a unit vector");
                }

                real e_xy = dot(e_x,e_y);
                real e_xz = dot(e_x,e_z);
                real e_yz = dot(e_y,e_z);

                if (abs(e_xy) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x and e_y are not orthogonal");
                }
                if (abs(e_xz) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x and e_z are not orthogonal");
                }
                if (abs(e_yz) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_y and e_z are not orthogonal");
                }

                tensor3 R_H;

                R_H.xx = e_x.x;
                R_H.yx = e_x.y;
                R_H.zx = e_x.z;

                R_H.xy = e_y.x;
                R_H.yy = e_y.y;
                R_H.zy = e_y.z;

                R_H.xz= e_z.x;
                R_H.yz= e_z.y;
                R_H.zz= e_z.z;

                return R_H;
            }

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real Eb;
                typename potential::params helixParams;

                tensor3 R_H;
            };

            struct PairParameters{

                real Eb;
                typename potential::params helixParams;

                tensor3 R_H;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.Eb  = real(info.at("Eb"));

                param.helixParams = potential::readParamsMap(info);

                real3 e_x = real3(info.at("e_x"));
                real3 e_y = real3(info.at("e_y"));
                real3 e_z = real3(info.at("e_z"));

                param.R_H = loadR_H(e_x,e_y,e_z);

                System::log<System::MESSAGE>("[HelixParameters] Added parameter Eb: %f",param.Eb);

                System::log<System::MESSAGE>("              %f,%f,%f",param.R_H.xx,param.R_H.xy,param.R_H.xz);
                System::log<System::MESSAGE>("[Helix] R_H = %f,%f,%f",param.R_H.yx,param.R_H.yy,param.R_H.yz);
                System::log<System::MESSAGE>("              %f,%f,%f",param.R_H.zx,param.R_H.zy,param.R_H.zz);

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.Eb  = in_par.Eb;

                param.helixParams = in_par.helixParams;

                param.R_H = in_par.R_H;

                return param;
            }
        };

    }

}}}}
