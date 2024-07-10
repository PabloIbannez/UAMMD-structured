#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/DistanceSwitch.cuh"
#include "Interactor/BasicPotentials/Harmonic.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    namespace Helix{

        static inline __device__ real3 ex_optimal(const Quat& qs,
                                                  const tensor3 R_H){
            const real3 ex_opt = rotateVector(qs,make_real3(R_H.xx,R_H.yx,R_H.zx));
            return ex_opt;
        }

        static inline __device__ real3 ey_optimal(const Quat& qs,
                                                  const tensor3 R_H){
            const real3 ey_opt = rotateVector(qs,make_real3(R_H.xy,R_H.yy,R_H.zy));
            return ey_opt;
        }

        static inline __device__ real cos_theta(const Quat& qs, const Quat& qe,
                                                const tensor3 R_H){

            const real3 ex     = qe.getVx();
            const real3 ex_opt = ex_optimal(qs,R_H);

            real ctheta = dot(ex,ex_opt);

            ctheta = fminf(ctheta, real(1.0));
            ctheta = fmaxf(ctheta,-real(1.0));

            return ctheta;
        }

        static inline __device__ real cos_phi(const Quat& qs, const Quat& qe,
                                              const tensor3 R_H){

            const real3 ey     = qe.getVy();
            const real3 ey_opt = ey_optimal(qs,R_H);

            real cphi = dot(ey,ey_opt);

            cphi = fminf(cphi, real(1.0));
            cphi = fmaxf(cphi,-real(1.0));

            return cphi;
        }

        struct exponential{

            struct potential{

                struct cosParameters{
                    real K;
                };

                static inline __device__ real cosEnergy(const Quat& qs, const Quat& qe,
                                                        const real& cosine,const cosParameters& params){

                    const real K = params.K;
                    if(K == real(0.0)) return real(0.0);

                    const real A   = expf(-real(2.0)*K);
                    const real den = real(1.0)-A;        // Denominator

                    const real e = (expf(-K*(real(1.0)-cosine))-A)/den;

                    return e;
                }

                static inline __device__ real cosEnergyDerivative(const Quat& qs, const Quat& qe,
                                                                  const real& cosine,const cosParameters& params){

                    const real K = params.K;
                    if(K == real(0.0)) return real(0.0);

                    const real A   = expf(-real(2.0)*K);
                    const real den = real(1.0)-A;        // Denominator
                    const real dU_dcos_theta = K*expf(-K*(real(1.0)-cosine))/den;

                    return dU_dcos_theta;
                }

                struct dstParameters{
                    real rc;
                    real K;
                };

                static inline __device__ real dstEnergy(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchExponential::energy(dr,r2,real(-1.0),params.rc,params.K);
                }

                static inline __device__ real3 dstForce(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchExponential::force(dr,r2,real(-1.0),params.rc,params.K);
                }

            };

            struct params{
                typename potential::dstParameters dstParams;
                typename potential::cosParameters thetaParams;
                typename potential::cosParameters phiParams;
            };

            static params readParams(DataEntry& data, bool loadDst = false){
                params p;

                p.thetaParams.K = data.getParameter<real>("Ka");
                p.phiParams.K   = data.getParameter<real>("Kd");

                System::log<System::MESSAGE>("[Helix] Ka = %f", p.thetaParams.K);
                System::log<System::MESSAGE>("[Helix] Kd = %f", p.phiParams.K);

                if(loadDst){
                    p.dstParams.rc = data.getParameter<real>("rc");
                    p.dstParams.K  = data.getParameter<real>("Kb");

                    System::log<System::MESSAGE>("[Helix] rc = %f", p.dstParams.rc);
                    System::log<System::MESSAGE>("[Helix] Kb = %f", p.dstParams.K);
                }

                return p;
            }

            template <typename T>
            static params readParamsMap(std::map<std::string,T>& info,std::string suffix = ""){
                params p;

                p.dstParams.rc = real(info.at("rc"+suffix));
                p.dstParams.K  = real(info.at("Kb"+suffix));

                p.thetaParams.K = real(info.at("Ka"+suffix));
                p.phiParams.K   = real(info.at("Kd"+suffix));

                System::log<System::MESSAGE>("[Helix] rc%s = %f", suffix.c_str(), p.dstParams.rc);
                System::log<System::MESSAGE>("[Helix] Kb%s = %f", suffix.c_str(), p.dstParams.K);

                System::log<System::MESSAGE>("[Helix] Ka%s = %f", suffix.c_str(), p.thetaParams.K);
                System::log<System::MESSAGE>("[Helix] Kd%s = %f", suffix.c_str(), p.phiParams.K);

                return p;
            }
        };

        struct cosine{

            struct potential{

                struct cosParameters{
                    real cos_angle_start;
                    real cos_angle_end;
                };

                static inline __device__ real cosEnergy(const Quat& qs, const Quat& qe,
                                                        const real& cosine,const cosParameters& params){

                    const real cos_angle_start = params.cos_angle_start;
                    const real cos_angle_end = params.cos_angle_end;

                    real swt;
                    if(cosine > cos_angle_start){
                        swt = real(-1.0);
                    } else if(cosine < cos_angle_end){
                        swt = real(1.0);
                    } else {
                        real norm = (cosine-cos_angle_start)/(cos_angle_end-cos_angle_start);

                        swt = cosf(norm*real(M_PI));

                        swt = fminf(swt, real(1.0));
                        swt = fmaxf(swt,-real(1.0));

                        swt = -swt;
                    }

                    swt = (real(1.0)-swt)*real(0.5);

                    return swt;
                }

                static inline __device__ real cosEnergyDerivative(const Quat& qs, const Quat& qe,
                                                                  const real& cosine,const cosParameters& params){

                    const real cos_angle_start = params.cos_angle_start;
                    const real cos_angle_end = params.cos_angle_end;

                    real dswt_dcos;

                    if(cosine > cos_angle_start){
                        dswt_dcos = real(0.0);
                    } else if(cosine < cos_angle_end){
                        dswt_dcos = real(0.0);
                    } else {
                        real norm = (cosine-cos_angle_start)/(cos_angle_end-cos_angle_start);

                        dswt_dcos = sinf(norm*real(M_PI));

                        dswt_dcos = fminf(dswt_dcos, real(1.0));
                        dswt_dcos = fmaxf(dswt_dcos,-real(1.0));

                        dswt_dcos = dswt_dcos*real(M_PI)/(cos_angle_end-cos_angle_start);
                    }

                    dswt_dcos = -dswt_dcos*real(0.5);

                    return dswt_dcos;
                }

                struct dstParameters{
                    real r_start;
                    real rc;
                };

                static inline __device__ real dstEnergy(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchCosine::energy(dr,r2,real(-1.0),params.r_start,params.rc);
                }

                static inline __device__ real3 dstForce(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchCosine::force(dr,r2,real(-1.0),params.r_start,params.rc);
                }

            };

            struct params{
                typename potential::dstParameters dstParams;
                typename potential::cosParameters thetaParams;
                typename potential::cosParameters phiParams;
            };

            static params readParams(DataEntry& data, bool loadDst = false){
                params p;

                p.thetaParams.cos_angle_start = data.getParameter<real>("theta_start");
                p.thetaParams.cos_angle_end   = data.getParameter<real>("theta_end");

                p.phiParams.cos_angle_start = data.getParameter<real>("phi_start");
                p.phiParams.cos_angle_end   = data.getParameter<real>("phi_end");

                System::log<System::MESSAGE>("[Helix] theta_start = %f", p.thetaParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] theta_end   = %f", p.thetaParams.cos_angle_end);

                System::log<System::MESSAGE>("[Helix] phi_start = %f", p.phiParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] phi_end   = %f", p.phiParams.cos_angle_end);

                p.thetaParams.cos_angle_start = cos(p.thetaParams.cos_angle_start);
                p.thetaParams.cos_angle_start = min(p.thetaParams.cos_angle_start, real(1.0));
                p.thetaParams.cos_angle_start = max(p.thetaParams.cos_angle_start,-real(1.0));

                p.thetaParams.cos_angle_end   = cos(p.thetaParams.cos_angle_end);
                p.thetaParams.cos_angle_end   = min(p.thetaParams.cos_angle_end, real(1.0));
                p.thetaParams.cos_angle_end   = max(p.thetaParams.cos_angle_end,-real(1.0));

                p.phiParams.cos_angle_start = cos(p.phiParams.cos_angle_start);
                p.phiParams.cos_angle_start = min(p.phiParams.cos_angle_start, real(1.0));
                p.phiParams.cos_angle_start = max(p.phiParams.cos_angle_start,-real(1.0));

                p.phiParams.cos_angle_end   = cos(p.phiParams.cos_angle_end);
                p.phiParams.cos_angle_end   = min(p.phiParams.cos_angle_end, real(1.0));
                p.phiParams.cos_angle_end   = max(p.phiParams.cos_angle_end,-real(1.0));

                if(loadDst){
                    p.dstParams.r_start = data.getParameter<real>("r_start");
                    p.dstParams.rc      = data.getParameter<real>("rc");

                    System::log<System::MESSAGE>("[Helix] r_start = %f", p.dstParams.r_start);
                    System::log<System::MESSAGE>("[Helix] rc      = %f", p.dstParams.rc);
                }

                return p;
            }

            template <typename T>
            static params readParamsMap(std::map<std::string,T>& info,std::string suffix = ""){
                params p;

                p.dstParams.rc = real(info.at("rc"+suffix));
                p.dstParams.r_start = real(info.at("r_start"+suffix));

                p.thetaParams.cos_angle_start = real(info.at("theta_start"+suffix));
                p.thetaParams.cos_angle_end   = real(info.at("theta_end"+suffix));

                p.phiParams.cos_angle_start = real(info.at("phi_start"+suffix));
                p.phiParams.cos_angle_end   = real(info.at("phi_end"+suffix));

                System::log<System::MESSAGE>("[Helix] rc%s = %f", suffix.c_str(), p.dstParams.rc);
                System::log<System::MESSAGE>("[Helix] r_start%s = %f", suffix.c_str(), p.dstParams.r_start);

                System::log<System::MESSAGE>("[Helix] theta_start%s = %f", suffix.c_str(), p.thetaParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] theta_end%s   = %f", suffix.c_str(), p.thetaParams.cos_angle_end);

                System::log<System::MESSAGE>("[Helix] phi_start%s = %f", suffix.c_str(), p.phiParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] phi_end%s   = %f", suffix.c_str(), p.phiParams.cos_angle_end);

                p.thetaParams.cos_angle_start = cos(p.thetaParams.cos_angle_start);
                p.thetaParams.cos_angle_start = min(p.thetaParams.cos_angle_start, real(1.0));
                p.thetaParams.cos_angle_start = max(p.thetaParams.cos_angle_start,-real(1.0));

                p.thetaParams.cos_angle_end   = cos(p.thetaParams.cos_angle_end);
                p.thetaParams.cos_angle_end   = min(p.thetaParams.cos_angle_end, real(1.0));
                p.thetaParams.cos_angle_end   = max(p.thetaParams.cos_angle_end,-real(1.0));

                p.phiParams.cos_angle_start = cos(p.phiParams.cos_angle_start);
                p.phiParams.cos_angle_start = min(p.phiParams.cos_angle_start, real(1.0));
                p.phiParams.cos_angle_start = max(p.phiParams.cos_angle_start,-real(1.0));

                p.phiParams.cos_angle_end   = cos(p.phiParams.cos_angle_end);
                p.phiParams.cos_angle_end   = min(p.phiParams.cos_angle_end, real(1.0));
                p.phiParams.cos_angle_end   = max(p.phiParams.cos_angle_end,-real(1.0));

                return p;
            }
        };

        template <typename potential>
        static inline __device__ real energyTheta(const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H, const typename potential::cosParameters& params){



            real ctheta = cos_theta(qs,qe,R_H);
            return potential::cosEnergy(qs,qe,ctheta,params);

        }

        template <typename potential>
        static inline __device__ real energyThetaForwardBackward(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                                 const typename potential::cosParameters& params){

            const real energyForward  = energyTheta<potential>(qs,qe,R_H,params);
            const real energyBackward = energyTheta<potential>(qe,qs,R_H.transpose(),params);

            return (energyForward+energyBackward)*real(0.5);
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativeTheta(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                          const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dsx
            // Second row: dU/dsy
            // Third row:  dU/dsz

            const real ctheta = cos_theta(qs,qe,R_H);
            const real dU_dcos_theta = potential::cosEnergyDerivative(qs,qe,ctheta,params);

            const real3 ex = qe.getVx();
            const real3 dcos_theta_dsx = ex*R_H.xx;
            const real3 dcos_theta_dsy = ex*R_H.yx;
            const real3 dcos_theta_dsz = ex*R_H.zx;

            const real3 dU_dsx = dU_dcos_theta*dcos_theta_dsx;
            const real3 dU_dsy = dU_dcos_theta*dcos_theta_dsy;
            const real3 dU_dsz = dU_dcos_theta*dcos_theta_dsz;

            tensor3 dU_dS;

            dU_dS.xx = dU_dsx.x;
            dU_dS.xy = dU_dsx.y;
            dU_dS.xz = dU_dsx.z;

            dU_dS.yx = dU_dsy.x;
            dU_dS.yy = dU_dsy.y;
            dU_dS.yz = dU_dsy.z;

            dU_dS.zx = dU_dsz.x;
            dU_dS.zy = dU_dsz.y;
            dU_dS.zz = dU_dsz.z;

            return dU_dS;
        }

        template<typename potential>
        static inline __device__ tensor3 eDerivativeTheta(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                          const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dex
            // Second row: dU/dey
            // Third row:  dU/dez

            const real ctheta = cos_theta(qs,qe,R_H);
            const real dU_dcos_theta = potential::cosEnergyDerivative(qs,qe,ctheta,params);

            const real3 dcos_theta_dex = ex_optimal(qs,R_H);
            //const real3 dcos_theta_dey = make_real3(0.0);
            //const real3 dcos_theta_dez = make_real3(0.0);

            const real3 dU_dex = dU_dcos_theta*dcos_theta_dex;

            tensor3 dU_dE;

            dU_dE.xx = dU_dex.x;
            dU_dE.xy = dU_dex.y;
            dU_dE.xz = dU_dex.z;

            dU_dE.yx = real(0.0);
            dU_dE.yy = real(0.0);
            dU_dE.yz = real(0.0);

            dU_dE.zx = real(0.0);
            dU_dE.zy = real(0.0);
            dU_dE.zz = real(0.0);

            return dU_dE;
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativeThetaForwardBackward(const Quat& qs, const Quat& qe,
                                                                         const tensor3 R_H,
                                                                         const typename potential::cosParameters& params){

            const tensor3 dU_ds_forward  = sDerivativeTheta<potential>(qs,qe,R_H,params);
            const tensor3 dU_ds_backward = eDerivativeTheta<potential>(qe,qs,R_H.transpose(),params);

            const tensor3 dU_ds = (dU_ds_forward + dU_ds_backward)*real(0.5);

            return dU_ds;
        }

        // PHI

        template <typename potential>
        static inline __device__ real energyPhi(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                const typename potential::cosParameters& params){

            real cphi = cos_phi(qs,qe,R_H);
            return potential::cosEnergy(qs,qe,cphi,params);
        }

        template <typename potential>
        static inline __device__ real energyPhiForwardBackward(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                               const typename potential::cosParameters& params){

            const real energyForward  = energyPhi<potential>(qs,qe,R_H,params);
            const real energyBackward = energyPhi<potential>(qe,qs,R_H.transpose(),params);

            return (energyForward+energyBackward)*real(0.5);
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativePhi(const Quat& qs, const Quat& qe,
                                                        const tensor3 R_H,
                                                        const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dsx
            // Second row: dU/dsy
            // Third row:  dU/dsz

            const real cphi = cos_phi(qs,qe,R_H);
            const real dU_dcos_phi = potential::cosEnergyDerivative(qs,qe,cphi,params);

            const real3 ey = qe.getVy();
            const real3 dcos_phi_dsx = ey*R_H.xy;
            const real3 dcos_phi_dsy = ey*R_H.yy;
            const real3 dcos_phi_dsz = ey*R_H.zy;

            const real3 dU_dsx = dU_dcos_phi*dcos_phi_dsx;
            const real3 dU_dsy = dU_dcos_phi*dcos_phi_dsy;
            const real3 dU_dsz = dU_dcos_phi*dcos_phi_dsz;

            tensor3 dU_dS;

            dU_dS.xx = dU_dsx.x;
            dU_dS.xy = dU_dsx.y;
            dU_dS.xz = dU_dsx.z;

            dU_dS.yx = dU_dsy.x;
            dU_dS.yy = dU_dsy.y;
            dU_dS.yz = dU_dsy.z;

            dU_dS.zx = dU_dsz.x;
            dU_dS.zy = dU_dsz.y;
            dU_dS.zz = dU_dsz.z;

            return dU_dS;
        }

        template <typename potential>
        static inline __device__ tensor3 eDerivativePhi(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                        const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dex
            // Second row: dU/dey
            // Third row:  dU/dez

            const real cphi = cos_phi(qs,qe,R_H);
            const real dU_dcos_phi = potential::cosEnergyDerivative(qs,qe,cphi,params);

            //const real3 dcos_phi_dex = make_real3(0.0);
            const real3 dcos_phi_dey   = ey_optimal(qs,R_H);
            //const real3 dcos_phi_dez = make_real3(0.0);

            const real3 dU_dey = dU_dcos_phi*dcos_phi_dey;

            tensor3 dU_dE;

            dU_dE.xx = real(0.0);
            dU_dE.xy = real(0.0);
            dU_dE.xz = real(0.0);

            dU_dE.yx = dU_dey.x;
            dU_dE.yy = dU_dey.y;
            dU_dE.yz = dU_dey.z;

            dU_dE.zx = real(0.0);
            dU_dE.zy = real(0.0);
            dU_dE.zz = real(0.0);

            return dU_dE;
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativePhiForwardBackward(const Quat& qs, const Quat& qe,
                                                                       const tensor3 R_H,
                                                                       const typename potential::cosParameters& params){


            const tensor3 dU_ds_forward  = sDerivativePhi<potential>(qs,qe,R_H,params);
            const tensor3 dU_ds_backward = eDerivativePhi<potential>(qe,qs,R_H.transpose(),params);

            const tensor3 dU_ds = (dU_ds_forward + dU_ds_backward)*real(0.5);

            return dU_ds;
        }

        // Combined theta and phi

        template <typename potential>
        static inline __device__ real orientationEnergy(const Quat& qs, const Quat& qe,
                                                        const tensor3 R_H,
                                                        const typename potential::cosParameters& paramsTheta,
                                                        const typename potential::cosParameters& paramsPhi){
            //The orientation energy is given by:
            // U = UthetaFB*UphiFB

            return energyThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*energyPhiForwardBackward<potential>(qs,qe,R_H,paramsPhi);
        }

        //For orientation only no forces are present

        template <typename potential>
        static inline __device__ real3 orientationTorque(const Quat& qs, const Quat& qe,
                                                         const tensor3 R_H,
                                                         const typename potential::cosParameters& paramsTheta,
                                                         const typename potential::cosParameters& paramsPhi){

            //The orientation energy is given by:
            // U = UthetaFB*UphiFB
            // This function returns the total torque over the particle s


            // First we compute derivates of the potential energy respect s
            // dU/ds = dUthetaFB/ds*UphiFB + UthetaFB*dUphiFB/ds

            const tensor3 dU_ds = sDerivativeThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*energyPhiForwardBackward<potential>(qs,qe,R_H,paramsPhi) +
                                  energyThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*sDerivativePhiForwardBackward<potential>(qs,qe,R_H,paramsPhi);

            // Then we compute the torque as:
            // T = -sx X dU/dsx - sy X dU/dsy - sz X dU/dsz

            real3 T = make_real3(0.0);

            const real3 dUdsx = make_real3(dU_ds.xx,dU_ds.xy,dU_ds.xz);
            const real3 dUdsy = make_real3(dU_ds.yx,dU_ds.yy,dU_ds.yz);
            const real3 dUdsz = make_real3(dU_ds.zx,dU_ds.zy,dU_ds.zz);

            const real3 sx = qs.getVx();
            const real3 sy = qs.getVy();
            const real3 sz = qs.getVz();

            T += -cross(sx,dUdsx);
            T += -cross(sy,dUdsy);
            T += -cross(sz,dUdsz);

            return T;

        }

        namespace Fixed {

            template <typename potential>
            static inline __device__ real energy(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const real& Kb,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){

                // dr is expected to be dr = e - s

                //The fixed energy is given by:
                // U = -E*UthetaFB*UphiFB+0.5*Kb*(r)^2
                // Where r is the distance between the points s and e
                // r^2 = dot(dr,dr)

                const real r2 = dot(dr,dr);

                return -E*orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi)
                        + Harmonic::energy(dr,r2,Kb,real(0.0));
            }

            static inline __device__ real3 force(const real3& dr,
                                                 const real& Kb){

                // Only the distance dependent of the fixed energy contributes to the force

                const real r2 = dot(dr,dr);

                if (r2 < real(1e-6)) {
                    return make_real3(0.0);
                }

                return Harmonic::force(dr,r2,Kb,real(0.0));
            }

            template <typename potential>
            static inline __device__ real3 torque(const real3& dr,
                                                  const real3& lps, // Local positions of s
                                                  const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H,
                                                  const real& Kb,
                                                  const typename potential::potential::cosParameters& paramsTheta,
                                                  const typename potential::potential::cosParameters& paramsPhi,
                                                  const real& E){

                // lps is a vector from the particle center to the point s

                // We first compute the torque due to the orientation

                const real3 T  = -E*orientationTorque<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
                const real3 fs =    force(dr,Kb);

                return T + cross(lps,fs);
            }

        } // namespace Fixed

        namespace Dynamic {

            template <typename potential>
            static inline __device__ real energy(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const typename potential::potential::dstParameters& paramsDst,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){
                // dr is expected to be dr = e - s

                //The dynamic energy is given by:
                // U = -E*Ubond*UthetaFB*UphiFB
                const real r2 = dot(dr,dr);

                const real Ub = -E*potential::potential::dstEnergy(dr,r2,paramsDst);
                const real Uo = orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);

                return Ub*Uo;
            }

            template <typename potential>
            static inline __device__ real3 force(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const typename potential::potential::dstParameters& paramsDst,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){

                // dr is expected to be dr = e - s

                //The dynamic energy is given by:
                // U = -E*Ubond*UthetaFB*UphiFB
                // Only the distance dependent (Ubond) of the dynamic energy contributes to the force
                // fs = -dU/dr = -E*(-dUbond/dr)*UthetaFB*UphiFB = -E*dstSwitch::force*orientationEnergy

                const real r2 = dot(dr,dr);

                return -E*potential::potential::dstForce(dr,r2,paramsDst)*orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
            }

            template <typename potential>
            static inline __device__ real3 torque(const real3& dr,
                                                  const real3& lps, // Local positions of s
                                                  const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H,
                                                  const typename potential::potential::dstParameters& paramsDst,
                                                  const typename potential::potential::cosParameters& paramsTheta,
                                                  const typename potential::potential::cosParameters& paramsPhi,
                                                  const real& E){

                // lps is a vector from the particle center to the point s

                // U = -E*Ubond*UthetaFB*UphiFB

                // dU/ds = -EUbond*(dUthetaFB/ds*UphiFB + UthetaFB*dUphiFB/ds)

                const real r2 = dot(dr,dr);

                const real Ubond = -E*potential::potential::dstEnergy(dr,r2,paramsDst);

                const real3 T  =  Ubond*orientationTorque<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
                const real3 fs =  force<potential>(dr,qs,qe,R_H,paramsDst,paramsTheta,paramsPhi,E);

                return T + cross(lps,fs);
            }
        }
    }


}}}}
