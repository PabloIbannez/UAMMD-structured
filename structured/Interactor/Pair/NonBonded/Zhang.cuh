#ifndef __ZHANG_POTENTIAL__
#define __ZHANG_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    namespace Zhang_ns{

        __device__ inline real3 getCurrentAxis(Quat q, const real3& axis){
            real3 currentAxis = make_real3(0.0);

            if (axis.x != real(0.0)){
                currentAxis += q.getVx()*axis.x;
            }

            if (axis.y != real(0.0)){
                currentAxis += q.getVy()*axis.y;
            }

            if (axis.z != real(0.0)){
                currentAxis += q.getVz()*axis.z;
            }

            return currentAxis;
        }

        /* Term of the force that depend on the distance used when the particles are near*/
        __device__ real repulsiveForceModulus(const real& rb, const real& eps, const real& dist){

            real rbdivr  = rb/dist;
            real rbdivr2 = rbdivr*rbdivr;
            real rbdivr3 = rbdivr2*rbdivr;
            real rbdivr5 = rbdivr3*rbdivr2;

            real force   = -real(4.0)*eps*(rbdivr5-rbdivr3)/rb;

            return force;
        }

        /* Function a(rij, ni, nj) of the references [1] and [2]*/
        __device__ real potential_a(const real3& rij, const real3& ni, const real3& nj,
                                    const real& stheta){

            real t1 =  dot(cross(ni, rij), cross(nj, rij));
            real t2 = -stheta*dot((nj-ni), rij);
            real t3 = -stheta*stheta;

            return t1+t2+t3;
        }

        /* Function \phi(rij, ni, nj) of the references*/
        __device__ real potential_phi(const real3& rij, const real3& ni, const real3& nj,
                                      const real& mu, const real& stheta){

            return real(1.0)+mu*(potential_a(rij, ni, nj, stheta)-real(1.0));
        }

        /* Gradient respect of the position coordinates of the function a(rij, ni, nj)*/
        __device__ real3 gradr_a(const real3& rij, const real3& ni, const real3& nj,
                                 const real& stheta, const real& dist){

            real ni_dot_rij = dot(ni,rij);
            real nj_dot_rij = dot(nj,rij);

            real3 t1 = real(2.0)*ni_dot_rij*nj_dot_rij*rij;
            real3 t2 = -ni_dot_rij*nj;
            real3 t3 = -nj_dot_rij*ni;
            real3 t4 = -stheta*(nj-ni-nj_dot_rij*rij+ni_dot_rij*rij);

            real3 result = (t1+t2+t3+t4)/dist;

            return result;
        }

        /* Gradient respect of the director vector of the function a(rij, ni, nj)*/
        __device__ real3 gradn_a(const real3& rij_unit, const real3& ni, const real3& nj,
                                 const real& stheta){

            real3 result = nj-(dot(nj,rij_unit)-stheta)*rij_unit;

            return result;
        }

        /* Interaction energy between two blobs when r<rb */
        __device__ real repulsiveEnergy(const Box& box,
                                        const real& dist2, const real& rb, const real& eps){

            if (dist2==real(0.0)) return real(0.0);

            real rb2_div_r2 = rb*rb/dist2;
            real energy = eps*rb2_div_r2*(rb2_div_r2-real(2.0));

            return energy;
        }

        __device__ real attractiveEnergy(const Box& box,
                                         const real& dist,
                                         const real& rb, const real& rc, const real& eps, const real& chi){

            if (dist<rb) return real(0.0);
            if (dist>rc) return real(0.0);

            real ang = real(0.5)*(dist-rb)/(rc-rb);
            real cang = cospif(ang);
            real energy = -eps*pow(cang,real(2.0)*chi);

            return energy;
        }

    }

    struct Zhang_{

        using ParametersPairsType   = typename BasicParameters::Pairs::Zhang;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersPairsType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;
        using PairParameters          = typename ParameterPairsHandler::PairParameters;

        //Computational data
        struct ComputationalData{

            real4* pos;
            real4* dir;

            Box box;

            ParametersPairsIterator paramPairIterator;

            real cutOff;

            real3 axis;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> atzParam;

            real cutOff;

            real3 axis;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.paramPairIterator = storage.atzParam->getPairIterator();

            computational.axis = storage.axis;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.atzParam = std::make_shared<ParameterPairsHandler>(gd, pg,
                                                                     data);

            storage.axis = data.getParameter<real3>("axis", {0.0, 0.0, 1.0});

            //Check axis is normalized

            real l = length(storage.axis);

            if (l<real(1.0)+real(1e-6) && l>real(1.0)-real(1e-6)){
                storage.axis = storage.axis/l;
                System::log<System::MESSAGE>("[Zhang] Axis: (%f, %f, %f)",
                                              storage.axis.x, storage.axis.y, storage.axis.z);
            } else{
                System::log<System::CRITICAL>("[Zhang] Axis is not normalized");
            }

            ///////////////////////////////////////////////////////////

            auto pairsParam = storage.atzParam->getPairParameters();

            real cutOff = 0.0;
            for(auto p : pairsParam){
                cutOff = std::max(cutOff, p.second.rc);
            }

            storage.cutOff = cutOff;
            System::log<System::MESSAGE>("[Zhang] cutOff: %f" ,storage.cutOff);

            return storage;
        }


        static inline __device__ real energy(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const auto param = computational.paramPairIterator(index_i, index_j);

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            Quat qi = computational.dir[index_i];
            Quat qj = computational.dir[index_j];

            const real3 ni = Zhang_ns::getCurrentAxis(qi, computational.axis);
            const real3 nj = Zhang_ns::getCurrentAxis(qj, computational.axis);

            const real r2 = dot(rij, rij);

            real e = real(0.0);

            const real r = sqrtf(r2);

            if(r<param.rc){

                real3 rij_unit = rij/r;

                real energyOrientation = Zhang_ns::potential_phi(rij_unit, ni, nj, param.mu, param.stheta);

                if (r>param.rb){
                    real energyPosition = Zhang_ns::attractiveEnergy(computational.box, r, param.rb, param.rc, param.epsilon, param.chi);
                    e += energyPosition*energyOrientation;
                } else {
                    real energyPosition = Zhang_ns::repulsiveEnergy(computational.box, r2, param.rb, param.epsilon);
                    e += energyPosition+(real(1.0)-energyOrientation)*param.epsilon;
                }
            }

            return e;

        }

        /* Function that computes the forces and torques between two particles when the particles are near (r<rb)*/
        static inline __device__ ForceTorque nearForcesAndTorques(const real3& rij_unit, const real3& ni, const real3& nj,
                                                                  const real& dist, const PairParameters& par){
            ForceTorque fandt;

            real3 f = make_real3(0);
            real3 t = make_real3(0);

            fandt.force  = make_real4(f,0.0);
            fandt.torque = make_real4(t,0.0);

            if (dist == real(0.0) ) return fandt;

            f = Zhang_ns::repulsiveForceModulus(par.rb, par.epsilon, dist)*rij_unit;

            real epsmu = par.epsilon*par.mu;

            f -=  epsmu*Zhang_ns::gradr_a(rij_unit, ni, nj, par.stheta, dist);
            t  = -epsmu*Zhang_ns::gradn_a(rij_unit, ni, nj, par.stheta);
            t  =  cross(t, ni);

            fandt.force  = make_real4(f,0.0);
            fandt.torque = make_real4(t,0.0);

            return fandt;
        }

        /* Function that computes the forces and torques between two particles when the particles are at a medium distance (rb<r<rc)*/
        static inline __device__ ForceTorque mediumForcesAndTorques(const real3& rij_unit, const real3& ni, const real3& nj,
                                                                    const real& dist, const PairParameters& par){

            ForceTorque fandt;

            real3 f = make_real3(0);
            real3 t = make_real3(0);

            real ang = real(0.5)*(dist-par.rb)/(par.rc-par.rb);
            real cang, sang;
            //sincospif(ang, &sang, &cang);
            sincospi(ang, &sang, &cang);

            real cangpow2chi        = powf(cang,real(2.0)*par.chi);
            real3 forceAtractive    = par.epsilon*real(M_PI)*par.chi/(par.rc-par.rb)*cangpow2chi*sang/cang*rij_unit;
            real potentialAtractive = -par.epsilon*cangpow2chi;

            f  = forceAtractive*Zhang_ns::potential_phi(rij_unit, ni, nj, par.mu, par.stheta);
            f += potentialAtractive*par.mu*Zhang_ns::gradr_a(rij_unit, ni, nj, par.stheta, dist);

            t = par.mu*potentialAtractive*Zhang_ns::gradn_a(rij_unit, ni, nj, par.stheta);
            t = cross(t, ni);

            fandt.force  = make_real4(f,0.0);
            fandt.torque = make_real4(t,0.0);

            return fandt;
        }

        static inline __device__ ForceTorque forceTorque(const int index_i,const int index_j,
                                                         const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const auto param = computational.paramPairIterator(index_i, index_j);

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            Quat qi = computational.dir[index_i];
            Quat qj = computational.dir[index_j];

            const real3 ni = Zhang_ns::getCurrentAxis(qi, computational.axis);
            const real3 nj = Zhang_ns::getCurrentAxis(qj, computational.axis);

            const real r2 = dot(rij, rij);

            ForceTorque forceTorque;

            forceTorque.force  = make_real4(0.0);
            forceTorque.torque = make_real4(0.0);

            const real r = sqrtf(r2);

            if(r<param.rc){

                if (r == real(0.0)) return forceTorque;

                if       (r<param.rb){
                    const real3 rij_unit = rij/r;
                    forceTorque = nearForcesAndTorques(rij_unit, ni, nj, r, param);
                }else if (r<param.rc){
                    const real3 rij_unit = rij/r;
                    forceTorque = mediumForcesAndTorques(rij_unit, ni, nj, r, param);
                }
            }

            return forceTorque;

        }
    };

    using Zhang     = NonBondedTorque_<Zhang_>;

}}}}

#endif
