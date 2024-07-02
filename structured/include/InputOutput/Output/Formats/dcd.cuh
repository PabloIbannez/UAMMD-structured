#ifndef __DCD__
#define __DCD__

namespace uammd{
namespace structured{
namespace dcd{

    namespace dcd_ns{

        inline
        void write_int(std::ofstream& out,int i){
            out.write(reinterpret_cast<char*>(&i),sizeof(int));
        }

        inline
        void write_float(std::ofstream& out,float f){
            out.write(reinterpret_cast<char*>(&f),sizeof(float));
        }

        inline
        void pad(char *s,int len){
            int curlen;
            int i;

            curlen = strlen(s);
            for (i = curlen; i < len; i++){
                s[i] = ' ';
            }

            s[len] = '\0';
        }
    }

    inline
    void WriteDCDheader(std::shared_ptr<ParticleGroup> pg,
                        int start,
                        int interval,
                        std::ofstream& out){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem();

        int N = pg->getNumberParticles();

        float float_buffer;
        char  char_buffer[255];

        dcd_ns::write_int(out,84);

        sprintf(char_buffer,"CORD");
        out.write(reinterpret_cast<char*>(&char_buffer),sizeof(char)*4);

        dcd_ns::write_int(out, 0);        /* NFILE - Number of frames  */
        dcd_ns::write_int(out, start   ); /* NPRIV - Starting timestep of DCD file */
        dcd_ns::write_int(out, interval); /* NSAVC - Timesteps between DCD saves */
        dcd_ns::write_int(out, 0);        /* NSTEP - Number of timesteps */

        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);

        float_buffer = 1.0;
        out.write(reinterpret_cast<char*>(&float_buffer),sizeof(float_buffer));

        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);

        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);

        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);
        dcd_ns::write_int(out, 0);

        dcd_ns::write_int(out, 24);
        dcd_ns::write_int(out, 84);
        dcd_ns::write_int(out, 164);
        dcd_ns::write_int(out, 2);

        sprintf(char_buffer,"REMARK");
        dcd_ns::pad(char_buffer,80);
        out.write(reinterpret_cast<char*>(&char_buffer),sizeof(char)*80);

	    sprintf(char_buffer,"REMARKS DATE: %i", -1);
        dcd_ns::pad(char_buffer,80);
        out.write(reinterpret_cast<char*>(&char_buffer),sizeof(char)*80);

        dcd_ns::write_int(out, 164);
        dcd_ns::write_int(out, 4);
        dcd_ns::write_int(out, N);
        dcd_ns::write_int(out, 4);
    }

    inline
    void WriteDCDstep(std::shared_ptr<ParticleGroup> pg,
                      Box box,
                      int frame,
                      int step,
                      std::ofstream& out){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem();

        int N = pg->getNumberParticles();

        auto id    = pd->getId(access::location::cpu,
                               access::mode::read);

        auto pos   = pd->getPos(access::location::cpu,
                                access::mode::read);

        auto groupIndex  = pg->getIndexIterator(access::location::cpu);
        auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

        std::map<int,int> id_index;
        fori(0,pg->getNumberParticles()){
            int id_   = id[groupIndex[i]];
            int index = sortedIndex[id_];

            id_index[id_]=index;
        }

        dcd_ns::write_int(out, N*4);

        for(const auto& ii : id_index){

            int index = ii.second;

            real4 pc = pos[index];
            real3 p  = box.apply_pbc(make_real3(pc));
            dcd_ns::write_float(out, p.x);

        }
        dcd_ns::write_int(out, N*4);

        dcd_ns::write_int(out, N*4);
        for(const auto& ii : id_index){

            int index = ii.second;

            real4 pc = pos[index];
            real3 p  = box.apply_pbc(make_real3(pc));

            dcd_ns::write_float(out, p.y);

        }
        dcd_ns::write_int(out, N*4);

        dcd_ns::write_int(out, N*4);
        for(const auto& ii : id_index){

            int index = ii.second;

            real4 pc = pos[index];
            real3 p  = box.apply_pbc(make_real3(pc));

            dcd_ns::write_float(out, p.z);

        }
        dcd_ns::write_int(out, N*4);

        //Update dcd_head
        out.seekp(8L , std::ios::beg);
        dcd_ns::write_int(out, frame);
        out.seekp(20L, std::ios::beg);
        dcd_ns::write_int(out, step);

        out.seekp(0, std::ios::end);
    }

}}}

#endif
