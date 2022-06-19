#ifndef __OUTPUT_STATE__
#define __OUTPUT_STATE__

namespace uammd{
namespace structured{
namespace InputOutput{
namespace Output{

void WriteCoord(std::shared_ptr<ParticleGroup> pg,
                Box box,
                std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();
    
    auto id = pd->getId(access::location::cpu, 
                        access::mode::read);

    auto pos   = pd->getPos(access::location::cpu, 
                            access::mode::read);
    
    auto vel   = pd->getVel(access::location::cpu, 
                            access::mode::read);
    
    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
    
    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;

        real4 pc = pos[index];
        //Check nan,inf ...
        if(std::isfinite(pc.x) and 
           std::isfinite(pc.y) and 
           std::isfinite(pc.z) and 
           std::isfinite(pc.w)){

            real3 p = box.apply_pbc(make_real3(pc));

            real3 v  = vel[index];

            int idc  =  id[index];
            
            out << idc << " " 
                << p   << " "
                << v   << std::endl;
        } else {
            sys->log<System::CRITICAL>("[WriteCoord] Nan or inf values found");
        }
    }
}
        
void WriteSP(std::shared_ptr<ParticleGroup> pg,
             Box box,
             std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();
    
    auto id   = pd->getId(access::location::cpu, 
                          access::mode::read);

    auto pos   = pd->getPos(access::location::cpu, 
                            access::mode::read);
    
    auto rad   = pd->getRadius(access::location::cpu, 
                               access::mode::read);
    
    /*
    auto chg   = pd->getCharge(access::location::cpu, 
                               access::mode::read);
                               */
    
    /*
    auto mdl   = pd->getModelId(access::location::cpu, 
                                access::mode::read);
    */

    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

    out <<"#Lx="<<box.boxSize.x*0.5
        <<";Ly="<<box.boxSize.y*0.5
        <<";Lz="<<box.boxSize.z*0.5<<";"<< std::endl;
    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;

        real4 pc = pos[index];
        real3 p = box.apply_pbc(make_real3(pc));
        int  type   = int(pc.w);
        real radius = rad[index];
        //int  m      = mdl[index];
            
        out << std::left 
            << std::setw(6) 
            << p      << " " 
            << radius << " " 
            << type   << std::endl;

        /*
        if(type==0){

            if(p.x > real(0.0)) continue;

            out << std::left 
                << std::setw(6) 
                << p      << " " 
                << radius << " " 
                << type   << std::endl;

        } else {
            out << std::left 
                << std::setw(6) 
                << p      << " " 
                << radius << " " 
                << type   << std::endl;
        }*/
    }
}

void WriteVelocity(std::shared_ptr<ParticleGroup> pg,
                   std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();
    
    auto id   = pd->getId(access::location::cpu, 
                          access::mode::read);
    
    auto vel   = pd->getVel(access::location::cpu, 
                            access::mode::read);
    
    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;

        real3 velocity = vel[index];

        out << std::left 
            << std::setw(6) 
            << velocity.x   << " " 
            << std::setw(6) 
            << velocity.y   << " " 
            << std::setw(6) 
            << velocity.z   << " "  
            << std::setw(6) 
            << sqrt(dot(velocity,velocity)) << std::endl;
    }
}
        
void WriteXYZ(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();
    
    auto id   = pd->getId(access::location::cpu, 
                          access::mode::read);

    auto pos   = pd->getPos(access::location::cpu, 
                            access::mode::read);

    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

    out << pg->getNumberParticles() << std::endl;
    out << std::endl;

    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;

        real4 pc = pos[index];
        real3 p = box.apply_pbc(make_real3(pc));
        int  type   = int(pc.w);

        out << type   << " " 
            << p      << std::endl;
    }

}
        
void WriteITPV(std::shared_ptr<ParticleGroup> pg,
               Box box,
               std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();

    auto id  = pd->getId(access::location::cpu, 
                         access::mode::read);
    
    auto pos = pd->getPos(access::location::cpu, 
                            access::mode::read);
    
    auto vel = pd->getVel(access::location::cpu, 
                          access::mode::read);
    
    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
    
    out << pg->getNumberParticles() << std::endl;
    out << std::endl;
    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;
    
        int   id_ = id[index];
        real4 pc = pos[index];
        real3 p = box.apply_pbc(make_real3(pc));
        int  type   = int(pc.w);
        real3 v   = vel[index];
        
        out << id_   << " " 
                   << type  << " " 
                   << p     << " " 
                   << v     << std::endl;
    }
}

void WriteITPD(std::shared_ptr<ParticleGroup> pg,
               Box box,
               std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();

    auto id  = pd->getId(access::location::cpu, 
                         access::mode::read);
    
    auto pos = pd->getPos(access::location::cpu, 
                            access::mode::read);
    
    auto dir = pd->getDir(access::location::cpu, 
                          access::mode::read);
    
    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
    
    out << pg->getNumberParticles() << std::endl;
    out << std::endl;
    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;
    
        int   id_ = id[index];
        real4 pc  = pos[index];
        real3 p   = box.apply_pbc(make_real3(pc));
        int  type = int(pc.w);
        real4 d   = dir[index];
        
        out << id_   << " " 
                   << type  << " " 
                   << p     << " " 
                   << d     << std::endl;
    }
}

void WritePDB(std::shared_ptr<ParticleGroup> pg,
              Box box,
              int frame,
              std::ofstream& out){
    
    auto pd  = pg->getParticleData(); 
    auto sys = pd->getSystem();

    auto id   = pd->getId(access::location::cpu, 
                          access::mode::read);

    auto pos   = pd->getPos(access::location::cpu, 
                            access::mode::read);

    auto chainId = pd->getChainId(access::location::cpu, 
                                  access::mode::read);
    auto resId   = pd->getResId(access::location::cpu, 
                                access::mode::read);

    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

    out   << std::left  << std::fixed <<
        std::setw(6) << "MODEL"    <<
        "    "                     <<
        std::setw(4) << frame  <<
        std::endl;

    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    for(const auto& ii : id_index){

        int index = ii.second;

        real4 pc = pos[index];
        real3 p = box.apply_pbc(make_real3(pc));

        int res = resId[index];
        int ch  = chainId[index];

        int id_ = ii.first;
        if(id_<0){id_=0;}
        if(id_>99999){id_=99999;}

        if(res<0){res=0;}
        if(res>9999){res=9999;}

        if(ch<0){ch=0;}
        if(ch>9){ch=9;}

        std::string name = "CA";

        out << std::left << std::fixed           <<
            std::setw(6) << "ATOM"               <<
            std::right                           <<
            std::setw(5) << id_                  <<
            " "                                   ;

        out << std::left << std::fixed <<" "     <<
            std::setw(3) << name                ;

        out << std::left << std::fixed                  <<
            std::setw(1) << " "        <<
            std::setw(3) << std::to_string(int(pc.w))<<
            " "                                      <<
            std::setw(1) << ch                       <<
            std::right                               <<
            std::setw(4) << res                      <<
            std::setw(1) << " "                      <<
            "   "                                    <<
            std::setprecision(3)                     <<
            std::setw(8) << p.x                      <<
            std::setw(8) << p.y                      <<
            std::setw(8) << p.z                      <<
            std::setprecision(2)                     <<
            std::setw(6) << 1.0                      <<
            std::setw(6) << 0.0                      <<
            "          "                             <<
            std::setw(2) << ""    ;

        out << std::left << std::fixed <<
            std::setw(2) << int(1);

        out << std::endl;
    }

    out << std::left    << std::fixed <<
        std::setw(6) << "ENDMDL"   <<
        std::endl;
}

void WritePDB(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out){

    WritePDB(pg,box,0,out);

}
        
        
void WriteDCDheader(std::shared_ptr<ParticleGroup> pg,
                    int start,
                    int interval,
                    std::ofstream& out){
    
    dcd::WriteDCDheader(pg,
                        start,
                        interval,
                        out);
}

void WriteDCD(std::shared_ptr<ParticleGroup> pg,
              Box box,
              int frame,int step,
              std::ofstream& out){

            dcd::WriteDCDstep(pg,
                              box,
                              frame,step,
                              out); 
}

void WriteLAMMPS(std::shared_ptr<ParticleGroup> pg,
                 Box box,
                 real t,
                 std::ofstream& out){

    out<<"ITEM: TIMESTEP"<<std::endl;
    out<<t<<std::endl;

    out<<"ITEM: NUMBER OF ATOMS"<<std::endl;
    out<<pg->getNumberParticles()<<std::endl;

    out<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl;
    out<<-0.5*box.boxSize.x<<" "<<0.5*box.boxSize.x<<std::endl;
    out<<-0.5*box.boxSize.y<<" "<<0.5*box.boxSize.y<<std::endl;
    out<<-0.5*box.boxSize.z<<" "<<0.5*box.boxSize.z<<std::endl;
        
    out<<"ITEM: ATOMS id type ";

    if(box.isPeriodicX()){
        out << " x";
    } else {
        out << " xu";
    }
    
    if(box.isPeriodicY()){
        out << " y";
    } else {
        out << " yu";
    }
    
    if(box.isPeriodicZ()){
        out << " z";
    } else {
        out << " zu";
    }

    out << std::endl;

    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
    
    auto id   = pd->getId(access::location::cpu, 
                          access::mode::read);
    
    auto pos   = pd->getPos(access::location::cpu, 
                            access::mode::read);

    std::map<int,int> id_index;
    fori(0,pg->getNumberParticles()){
        int id_   = id[groupIndex[i]];
        int index = sortedIndex[id_];

        id_index[id_]=index;
    }

    int i=0;
    for(const auto& ii : id_index){
        i++;

        int index = ii.second;

        real4 pc = pos[index];

        out << i << " " << int(pc.w)+1 << " " << box.apply_pbc(make_real3(pc)) << std::endl;
    
    }

}

}}}}

#endif
