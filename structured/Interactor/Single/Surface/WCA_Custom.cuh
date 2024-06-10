#ifndef __WCA_PLAIN_CYLINDER_POT__
#define __WCA_PLAIN_CYLINDER_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

    template <class Geometry>
    struct wcaCustomBase_{

        using ParametersSingleType   = typename BasicParameters::Single::LennardJones;
        using ParameterSingleHandler = typename structured::SingleParameterHandler<ParametersSingleType>;
        using ParametersSingleIterator = typename ParameterSingleHandler::SingleIterator;


        struct StorageData: public Geometry::StorageData{
            std::shared_ptr<ParameterSingleHandler> ljParam;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>                   gd,
                                                                 std::shared_ptr<ParticleGroup>  pg,
                                                                 DataEntry& data){
            StorageData storage;
	    
	    static_cast<typename Geometry::StorageData&>(storage) =
            Geometry::getStorageData(gd, pg, data);
            
	    storage.ljParam = std::make_shared<ParameterSingleHandler>(gd,pg,
                                                                     data);
            return storage;
        }

	struct ComputationalData: public Geometry::ComputationalData{
            real4* pos;
            ParametersSingleIterator paramSingleIterator;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData& storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){
            ComputationalData computational;

	    static_cast<typename Geometry::ComputationalData&>(computational) =
            Geometry::getComputationalData(gd, pg, storage, comp, st);

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.paramSingleIterator = storage.ljParam->getSingleIterator();

            return computational;
        }

        static inline __device__  real3 force(const int& index_i,
                                              ComputationalData computational){

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

	    real4 gradDistance = Geometry::gradientDistance(computational.pos[index_i],computational);
	    real r2wall = gradDistance.w;

	    if(r2wall > sigma*sigma*real(1.259921)){ return make_real3(0.0);}

	    real invr2wall = sigma*sigma/r2wall;
	    real invr6wall = invr2wall*invr2wall*invr2wall;
	    real invrwall = 1/sqrt(r2wall);

	    real f_pre = -real(4.0)*epsilon*(real(6.0)*invr6wall-real(12.0)*invr6wall*invr6wall)*invrwall;

	    real fx = f_pre*gradDistance.x;
	    real fy = f_pre*gradDistance.y;
	    real fz = f_pre*gradDistance.z;

            real3 f = make_real3(fx,fy,fz);

            return f;
        }

        static inline __device__  real energy(const int& index_i,
                                              ComputationalData computational){

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

	    real4 gradDistance = Geometry::gradientDistance(computational.pos[index_i],computational);
 	    real r2wall = gradDistance.w;

	    if(r2wall > sigma*sigma*real(1.259921)) return real(0);

	    real invr2wall = sigma*sigma/r2wall;
	    real invr6wall = invr2wall*invr2wall*invr2wall;

            real e = real(4.0)*epsilon*(invr6wall*invr6wall-invr6wall)+epsilon;

            return e;
        }
    };

    struct PlainCylinder_{

	struct StorageData{
		real z0;
		real Rcyl;
	};

	struct ComputationalData{
		real z0;
		real Rcyl;
	};

	static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>     gd,
                                                   std::shared_ptr<ParticleGroup>  pg,
                                                   DataEntry& data){
		StorageData storage;
		storage.z0   = data.getParameter<real>("plainPosition",0);
		storage.Rcyl = data.getParameter<real>("cylinderRadius");
		return storage;
	}

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData& storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){
                ComputationalData computational;
		computational.z0   = storage.z0;
		computational.Rcyl = storage.Rcyl;
		return computational;
	}


	//gradientDistance calculates the distance to the wall and the gradient of r_wall in cartessian coordinates (gradient needed to apply chain rule in the Force)
	//gradientDistance.w = Distance to the wall Squared
	//gradientDistance.(x,y,z) = partial derivation with respect (x,y,z)
	static inline __device__ real4 gradientDistance(real4 pos,
							ComputationalData computational){

	    //Custom parameters reading
	    real Rcyl = computational.Rcyl;
	    real z0 = computational.z0;
				
	    real r =  sqrt(pos.x*pos.x+pos.y*pos.y); //This Surface is defined in Cylindric coordinates
	    real3 drWall;

	    real r_cyl   = Rcyl - r;
	    real r_plain = pos.z-z0;

	    real dr;
	    real dz;
	    real r_wall_square;

	    if (pos.z<z0 && r<Rcyl){ //Under the plain Inside the cylinder
		    r_wall_square = r_cyl*r_cyl;
		    dr = real(-1.0);
		    dz = real(0.0);

	    } else if (pos.z>z0 && r>Rcyl){ //Above the plain outside the cylinder
		    r_wall_square = r_plain*r_plain;
		    dr = real(0.0);
		    dz = real(1.0);

	    } else if (pos.z>z0 && r<Rcyl){ //Above the plain insider the cylinder (corner)
		    r_wall_square = (r_cyl*r_cyl+r_plain*r_plain);
		    real rwall = sqrt(r_wall_square);

		    dr = (r-Rcyl)/rwall;
		    dz = (pos.z-z0)/rwall;

	    } else { //Particles never would eneter the wall so no force is defined here
		    return make_real4(0.0);
	    }

	    real dx = dr*pos.x/r;
	    real dy = dr*pos.y/r;

	    return make_real4(dx,dy,dz,r_wall_square);
	}



};

    using WCA_PlainCylinder = Surface_<wcaCustomBase_<PlainCylinder_>>;

}}}}

#endif
