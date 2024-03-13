#ifndef MEASURES_TRANSFORMS_CUH
#define MEASURES_TRANSFORMS_CUH

namespace uammd{
namespace structured{
namespace Measures{
namespace MeasuresTransforms{

    struct posX : public thrust::unary_function<int,real>
    {
        real4* pos;
        Box    box;

        posX(real4* pos,
             Box box):pos(pos),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                return box.apply_pbc(make_real3(pos[index])).x;
            }
    };

    struct posY : public thrust::unary_function<int,real>
    {
        real4* pos;
        Box    box;

        posY(real4* pos,
             Box box):pos(pos),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                return box.apply_pbc(make_real3(pos[index])).y;
            }
    };

    struct posZ : public thrust::unary_function<int,real>
    {
        real4* pos;
        Box    box;

        posZ(real4* pos,
             Box box):pos(pos),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                return box.apply_pbc(make_real3(pos[index])).z;
            }
    };

    struct totalDistance : public thrust::unary_function<int,real>
    {
        real4* pos;
        real3  point;
        Box    box;

        totalDistance(real4* pos,real3 point,Box box):pos(pos),point(point),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                real3 dr = box.apply_pbc(make_real3(pos[index])-point);
                return sqrt(dot(dr,dr));
            }
    };

    struct totalSquareDistance : public thrust::unary_function<int,real>
    {
        real4* pos;
        real3  point;
        Box    box;

        totalSquareDistance(real4* pos,real3 point,Box box):pos(pos),point(point),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                real3 dr = box.apply_pbc(make_real3(pos[index])-point);
                return dot(dr,dr);
            }
    };

    struct computeVirial : public thrust::unary_function<int,real>
    {
        real4* pos;
        real4* force;
        Box    box;

        computeVirial(real4* pos,real4* force,Box box):pos(pos),force(force),box(box){}

        __host__ __device__
            real operator()(int index) const
            {
                return dot(box.apply_pbc(make_real3(pos[index])),make_real3(force[index]));
            }
    };

    struct totalVirial : public thrust::unary_function<int,real>
    {
        real* virial;

        totalVirial(real* virial):virial(virial){}

        __host__ __device__
            real operator()(int index) const
            {
                return virial[index];
            }
    };

    struct totalStress : public thrust::unary_function<int,tensor3>
    {
        tensor3* stress;

        totalStress(tensor3* stress):stress(stress){}

        __host__ __device__
            tensor3 operator()(int index) const
            {
                return stress[index];
            }
    };

    struct kineticPressure : public thrust::unary_function<int,tensor3>
    {
        real* mass;
        real3* vel;

        kineticPressure(real* mass,
                real3* vel):mass(mass),
        vel(vel){}

        __host__ __device__
            tensor3 operator()(int index) const
            {
                return mass[index]*outer(vel[index],vel[index]);
            }
    };

    struct kineticEnergy : public thrust::unary_function<int,real>
    {
        real* mass;
        real3* vel;

        kineticEnergy(real* mass,
                real3* vel):mass(mass),
        vel(vel){}

        __host__ __device__
            real operator()(int index) const
            {
                return real(0.5)*mass[index]*dot(vel[index],vel[index]);
            }
    };

    struct potentialEnergy : public thrust::unary_function<int,real>
    {
        real* energy;

        potentialEnergy(real* energy):energy(energy){}

        __host__ __device__
            real operator()(int index) const
            {
                return energy[index];
            }
    };

    struct lambdaDerivative : public thrust::unary_function<int,real>
    {
        real* lambdaD;

        lambdaDerivative(real* lambdaD):lambdaD(lambdaD){}

        __host__ __device__
            real operator()(int index) const
            {
                return lambdaD[index];
            }
    };


    struct totalPos : public thrust::unary_function<int,real4>
    {
        real4* pos;

        totalPos(real4* pos):pos(pos){}

        __host__ __device__
            real4 operator()(int index) const
            {
                return pos[index];
            }
    };

    struct maxForce : public thrust::unary_function<int,real>
    {
        real4* force;

        maxForce(real4* force):force(force){}

        __host__ __device__
            real operator()(int index) const
            {
                real3 f = make_real3(force[index]);
                return sqrt(dot(f,f));
            }
    };

    struct totalForce : public thrust::unary_function<int,real4>
    {
        real4* force;

        totalForce(real4* force):force(force){}

        __host__ __device__
            real4 operator()(int index) const
            {
                return force[index];
            }
    };

    struct totalCharge : public thrust::unary_function<int,real>
    {
        real* charge;

        totalCharge(real* charge):charge(charge){}

        __host__ __device__
            real operator()(int index) const
            {
                return charge[index];
            }
    };

    struct totalMass : public thrust::unary_function<int,real>
    {
        real* mass;

        totalMass(real* mass):mass(mass){}

        __host__ __device__
            real operator()(int index) const
            {
                return mass[index];
            }
    };

    struct magneticMoment_vec : public thrust::unary_function<int,real3>
    {
        real4* dir;
        real4* magnetization;

        magneticMoment_vec(real4* dir, real4* magnetization):
            dir(dir), magnetization(magnetization){}

        __host__ __device__
            real3 operator()(int index) const
            {
                Quat diri = dir[index];
                real4 m_and_M = magnetization[index];
                real3 m = make_real3(m_and_M);
                // m is stored in the particle's frame, but we want it in the lab frame
                return m_and_M.w*rotateVector(diri, make_real3(m_and_M));
            }
    };

    struct magneticMoment_mod : public thrust::unary_function<int,real>
    {
        real4* magnetization;

        magneticMoment_mod(real4* magnetization):magnetization(magnetization){}

        __host__ __device__
            real operator()(int index) const
            {
                return magnetization[index].w;
            }
    };

    template<class targetType>
        struct weightedSum : public thrust::unary_function<int,targetType>
    {
        real* weights;
        targetType* target;

        weightedSum(real* weights,
                targetType* target):weights(weights),
        target(target){}

        __host__ __device__
            targetType operator()(int index) const
            {
                return weights[index]*target[index];
            }
    };

    struct angularMomentum : public thrust::unary_function<int,real3>
    {
        real* mass;
        real4* pos;
        real3* vel;

        real3 refp;
        real3 refv;

        Box box;

        angularMomentum(real* mass,
                        real4* pos,
                        real3* vel,
                        real3 refp,
                        real3 refv,
                        Box box):mass(mass),
                                 pos(pos),
                                 vel(vel),
                                 refp(refp),
                                 refv(refv),
                                 box(box){}

        __host__ __device__
            real3 operator()(int index) const
            {
                real3 p = box.apply_pbc(make_real3(pos[index])-refp);
                real3 v = vel[index]-refv;

                real3 L ={p.y*v.z-p.z*v.y,
                    p.z*v.x-p.x*v.z,
                    p.x*v.y-p.y*v.x};

                return L*mass[index];
            }
    };


    struct inertiaTensor : public thrust::unary_function<int,tensor3>
    {
        real* mass;
        real4* pos;

        real3 refp;
        Box box;

        inertiaTensor(real* mass,
                      real4* pos,
                      real3 refp,
                      Box box):mass(mass),
                               pos(pos),
                               refp(refp),
                               box(box){}

        __host__ __device__
            tensor3 operator()(int index) const
            {
                real3 p = box.apply_pbc(make_real3(pos[index])-refp);

                real x = p.x;
                real y = p.y;
                real z = p.z;

                real xx = x*x;
                real xy = x*y;
                real xz = x*z;
                real yy = y*y;
                real yz = y*z;
                real zz = z*z;

                tensor3 L =tensor3(yy+zz,-xy  ,-xz,
                                   -xy  ,xx+zz,-yz,
                                   -xz  ,-yz  ,xx+yy);

                return L*mass[index];
            }
    };


}}}}


#endif
