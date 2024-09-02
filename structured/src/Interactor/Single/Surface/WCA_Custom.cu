#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/Surface/Surface.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicParameters/Single/LennardJones.cuh"
#include "Utils/ParameterHandler/SingleParameterHandler.cuh"

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


                System::log<System::MESSAGE>("[WCA Custom] Loading WCA parameters");
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
                Box box;
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
                computational.box = gd->getEnsemble()->getBox();

                return computational;
            }

            static inline __device__  real3 force(const int& index_i,
                    ComputationalData computational){

                const real epsilon = computational.paramSingleIterator(index_i).epsilon;
                const real sigma   = computational.paramSingleIterator(index_i).sigma;
                const real3 posi   = computational.box.apply_pbc(make_real3(computational.pos[index_i]));

                real4 gradDistance = Geometry::gradientDistance(posi,computational);
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
                const real3 posi   = computational.box.apply_pbc(make_real3(computational.pos[index_i]));

                real4 gradDistance = Geometry::gradientDistance(posi,computational);
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
        static inline __device__ real4 gradientDistance(real3 pos,
                ComputationalData computational){

            //Custom parameters reading
            const real Rcyl = computational.Rcyl;
            const real z0 = computational.z0;

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

    struct Polyhedron_{

        struct StorageData{
            std::shared_ptr<thrust::device_vector<real3>> vertex;
            std::shared_ptr<thrust::device_vector<real3>> faces;
            int facesNumber;
        };

        struct ComputationalData{
            real3* vertex;
            real3* faces;
            int facesNumber;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>     gd,
                std::shared_ptr<ParticleGroup>  pg,
                DataEntry& data){

            //data structure is [[[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]],[[],[],[]],...]
            //labels = ["vertex1","vertex2","vertex3"] each three vertex of the same row form a face, vertex should be repeated, order define outside/inside

            System::log<System::MESSAGE>("[WCA Custom] Geometry loaded: Polyhedron");
            std::vector<real3> host_vertex = data.getParameter<std::vector<real3>>("vertex");
            std::vector<real3>  host_faces  = data.getParameter<std::vector<real3>>("faces");
            for(int i=0;i<host_vertex.size();++i){
                System::log<System::MESSAGE>("[WCA Polyhedron] readed vertex: [%f,%f,%f]",host_vertex[i].x,host_vertex[i].y,host_vertex[i].z);
            }
            for(int i=0;i<host_faces.size();++i){ //Comprobar que FACES no itera fuera de vertex
                System::log<System::MESSAGE>("[WCA Polyhedron] readed faces: [%d,%d,%d]",int(host_faces[i].x),int(host_faces[i].y),int(host_faces[i].z));
            }
            // Crear vectores de dispositivo a partir de vectores de host
            //thrust::device_vector<real3> device_vertex(host_vertex.begin(), host_vertex.end());
            //thrust::device_vector<real3>  device_faces (host_faces.begin() , host_faces.end());


            System::log<System::MESSAGE>("[WCA Polyhedron] Vertex and faces copied to GPU succesfully");
            StorageData storage;
            storage.facesNumber = host_faces.size();
            storage.vertex      = std::make_shared<thrust::device_vector<real3>>(host_vertex.begin(),host_vertex.end());
            storage.faces       = std::make_shared<thrust::device_vector<real3>>(host_faces.begin(),host_faces.end());

            // Obtener punteros a los datos en el dispositivo
            //storage.vertex = thrust::raw_pointer_cast(device_vertex.data());
            //storage.faces  = thrust::raw_pointer_cast(device_faces.data());

            return storage;
        }

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                std::shared_ptr<ParticleGroup> pg,
                const StorageData& storage,
                const Computables& comp,
                const cudaStream_t& st){
            ComputationalData computational;

            computational.vertex      = thrust::raw_pointer_cast(storage.vertex->data());
            computational.faces       = thrust::raw_pointer_cast(storage.faces->data());
            computational.facesNumber = storage.facesNumber;

            return computational;
        }


        static inline __device__ real anglePointVertices(real3 posi,
                real3 v0,
                real3 v1){
            real3 vector1 = v0 - posi;
            real3 vector2 = v1 - posi;

            real norm_vector2 = sqrt(vector2.x*vector2.x+vector2.y*vector2.y+vector2.z*vector2.z);
            real norm_vector1 = sqrt(vector1.x*vector1.x+vector1.y*vector1.y+vector1.z*vector1.z);

            real cos = dot(vector1,vector2)/(norm_vector1*norm_vector2);
            cos = std::min(real(1.0),cos);
            cos = std::max(real(-1.0),cos);

            return acos(cos);
        }


        static inline __device__ real4 distancePointToPlane(real3 posi,
                real3 v0,
                real3 v1,
                real3 v2){

            real distance;
            real3 gradient;

            //Distance to plane
            real3 normal = cross(v1-v0,v2-v0);
            normal = normal/sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
            real3 posi_to_face_vector = posi-v0;
            real distance_to_plane = dot(posi_to_face_vector,normal);
            real3 projection = posi - distance_to_plane*normal;

            //Cehck if projection is on the face
            real angle = 0;
            angle += anglePointVertices(projection,v0,v1);
            angle += anglePointVertices(projection,v1,v2);
            angle += anglePointVertices(projection,v2,v0);
            bool isPointInFace = abs(angle-real(2.0)*real(M_PI))<real(1e-4); //to avoid float error
                                                                             //bool isPointInFace = true;
            if(isPointInFace){
                distance = abs(distance_to_plane);
            }

            else{
                distance = std::numeric_limits<real>::max(); //Meter el float mas grande que exista
            }

            //Gradient
            gradient = normal*distance_to_plane/distance;

            return make_real4(gradient.x,gradient.y,gradient.z,distance*distance);


        }

        static inline __device__ real4 distancePointToEdge(real3 posi,
                real3 v0,
                real3 v1){

            real3 gradient;
            real distance;

            real3 edge_vector = v1 - v0;
            real3 posi_vector = posi - v0;
            real3 nearest_posi = make_real3(0.0);
            real t = dot(posi_vector,edge_vector)/(dot(edge_vector,edge_vector));
            if(t<real(0.0)){
                nearest_posi = v0;
            }
            else if(t>real(1.0)){
                nearest_posi = v1;
            }
            else{
                nearest_posi = v0 + t*edge_vector;
            }

            real3 dif = posi - nearest_posi;
            distance = dif.x*dif.x+dif.y*dif.y+dif.z*dif.z;
            //distance = std::numeric_limits<real>::max();
            gradient = dif/sqrt(distance);

            return make_real4(gradient.x,gradient.y,gradient.z,distance);
        }



        //gradientDistance calculates the distance to the wall and the gradient of r_wall in cartessian coordinates (gradient needed to apply chain rule in the Force)
        //gradientDistance.w = Distance to the wall Squared
        //gradientDistance.(x,y,z) = partial derivation with respect (x,y,z)
        static inline __device__ real4 gradientDistance(real3 pos,
                ComputationalData computational){

            real4 gradient_distance = make_real4(0.0,0.0,0.0,std::numeric_limits<real>::max());
            real4 aux;

            for (int i=0;i<computational.facesNumber;++i){
                int3 face = make_int3(computational.faces[i]);
                real3 v0 = computational.vertex[face.x];
                real3 v1 = computational.vertex[face.y];
                real3 v2 = computational.vertex[face.z];

                //Distance to faces
                aux = distancePointToPlane(pos,v0,v1,v2);
                if(aux.w<gradient_distance.w){
                    gradient_distance = aux;
                }

                //Distance to edges
                aux = distancePointToEdge(pos,v0,v1);
                if(aux.w<gradient_distance.w){
                    gradient_distance = aux;
                }
                aux = distancePointToEdge(pos,v1,v2);
                if(aux.w<gradient_distance.w){
                    gradient_distance = aux;
                }
                aux = distancePointToEdge(pos,v2,v0);
                if(aux.w<gradient_distance.w){
                    gradient_distance = aux;
                }
            }

            return gradient_distance;
        }
    };

    using WCA_Polyhedron = Surface_<wcaCustomBase_<Polyhedron_>>;
    using WCA_PlainCylinder = Surface_<wcaCustomBase_<PlainCylinder_>>;


}}}}

REGISTER_SINGLE_INTERACTOR(
    Surface,WCA_Polyhedron,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::Surface::WCA_Polyhedron>
)

REGISTER_SINGLE_INTERACTOR(
    Surface,WCA_PlainCylinder,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::Surface::WCA_PlainCylinder>
)
