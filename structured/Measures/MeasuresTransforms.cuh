#ifndef MEASURES_TRANSFORMS_CUH
#define MEASURES_TRANSFORMS_CUH

namespace uammd{
namespace structured{
namespace Measures{
namespace MeasuresTransforms{
    
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

        angularMomentum(real* mass,
                        real4* pos, 
                        real3* vel,
                        real3 refp,
                        real3 refv):mass(mass),
                                    pos(pos),
                                    vel(vel),
                                    refp(refp),
                                    refv(refv){}

        __host__ __device__
        real3 operator()(int index) const
        {
            real3 p = make_real3(pos[index])-refp;
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

        inertiaTensor(real* mass,
                      real4* pos,
                      real3 refp):mass(mass),
                                  pos(pos),
                                  refp(refp){}

        __host__ __device__
        tensor3 operator()(int index) const
        {
            real3 p = make_real3(pos[index])-refp;

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
