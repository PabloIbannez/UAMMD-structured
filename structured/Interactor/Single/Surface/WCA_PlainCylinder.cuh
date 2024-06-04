#ifndef __WCA_PLAIN_CYLINDER_POT__
#define __WCA_PLAIN_CYLINDER_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

    struct WCA_PlainCylinder_{

        using ParametersSingleType   = typename BasicParameters::Single::LennardJones;
        using ParameterSingleHandler = typename structured::SingleParameterHandler<ParametersSingleType>;

        using ParametersSingleIterator = typename ParameterSingleHandler::SingleIterator;

        struct ComputationalData{
            real4* pos;
            ParametersSingleIterator paramSingleIterator;
            real plainPosition;
	    real cylinderRadius;
        };

        //Potential parameters
        struct StorageData{
            std::shared_ptr<ParameterSingleHandler> ljParam;
            real plainPosition;
	    real cylinderRadius;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>                   gd,
                                                                 std::shared_ptr<ParticleGroup>  pg,
                                                                 DataEntry& data){

            StorageData storage;

            storage.ljParam = std::make_shared<ParameterSingleHandler>(gd,pg,
                                                                     data);

            storage.plainPosition = data.getParameter<real>("plainPosition",0);
	    storage.cylinderRadius = data.getParameter<real>("cylinderRadius");

            return storage;
        }


        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData& storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st) {

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.paramSingleIterator = storage.ljParam->getSingleIterator();

            computational.plainPosition = storage.plainPosition;

	    computational.cylinderRadius = storage.cylinderRadius;

            return computational;
        }

        //Storage data reader


        static inline __device__  real3 force(const int& index_i,
                                              ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;
	    const real z0      = computational.plainPosition;
	    const real Rcyl    = computational.cylinderRadius;


	    real r = sqrt(posi.x*posi.x+posi.y*posi.y);

	    real r_cyl   = Rcyl - r;
	    real r_plain = posi.z-z0;
	    real r2wall;
	    real dr;
	    real dz;
	    real rwall;

	    if (z0>posi.z && Rcyl>r){ //Under the plain Inside the cylinder
		    r2wall = r_cyl*r_cyl;
		    dr = real(-1.0);
		    dz = real(0.0);

	    } else if (z0<posi.z && Rcyl<r){ //Above the plain outside the cylinder
		    r2wall = r_plain*r_plain;
		    dr = real(0.0);
		    dz = real(1.0);

	    } else if (z0<posi.z && Rcyl>r){ //Above the plain insider the cylinder (corner)
		    r2wall = r_cyl*r_cyl+r_plain*r_plain;
		    rwall = sqrt(r2wall);
		    dr = (r-Rcyl)/rwall;
		    dz = (posi.z-z0)/rwall;


	    } else { //Particles never would eneter the wall so no force is defined here
		    return make_real3(0.0);
	    }

	    if(r2wall > sigma*sigma*real(1.259921)){ return make_real3(0.0);}

	    real invr2wall = sigma*sigma/r2wall;
	    real invr6wall = invr2wall*invr2wall*invr2wall;
	    real invrwall = 1/sqrt(r2wall);

	    real f_pre = -real(4.0)*epsilon*(real(6.0)*invr6wall-real(12.0)*invr6wall*invr6wall)*invrwall;

	    real fx = f_pre*dr*posi.x/r;
	    real fy = f_pre*dr*posi.y/r;
	    real fz = f_pre*dz;

            real3 f = make_real3(fx,fy,fz);

            return f;
        }

        static inline __device__  real energy(const int& index_i,
                                              ComputationalData computational){


	    real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;
	    const real z0      = computational.plainPosition;
	    const real Rcyl    = computational.cylinderRadius;

	    real r = sqrt(posi.x*posi.x+posi.y*posi.y);

	    real r_cyl   = Rcyl - r;
	    real r_plain = posi.z-z0;
	    real r2wall;


	    if (z0 > posi.z && Rcyl>r){ //Under the plain Inside the cylinder
		    r2wall = r_cyl*r_cyl;

	    } else if (z0 < posi.z && Rcyl<r){ //Above the plain outside the cylinder
		    r2wall = r_plain*r_plain;

	    } else if (z0 < posi.z && Rcyl>r){ //Above the plain insider the cylinder (corner)
		    r2wall = r_cyl*r_cyl+r_plain*r_plain;

	    } else { //Inside the wall (set r2 bigger than rcut so there is no force.
		    return real(0);
	    }

	    if(r2wall > sigma*sigma*real(1.259921)) return real(0);

	    real invr2wall = sigma*sigma/r2wall;
	    real invr6wall = invr2wall*invr2wall*invr2wall;

            real e = real(4.0)*epsilon*(invr6wall*invr6wall-invr6wall)+epsilon;

            return e;
        }
    };

    using WCA_PlainCylinder = Surface_<WCA_PlainCylinder_>;

}}}}

#endif
