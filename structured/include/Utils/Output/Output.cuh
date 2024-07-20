#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Formats/psf.cuh"
#include "Formats/dcd.cuh"

#include "utils/quaternion.cuh"

namespace uammd{
namespace structured{
namespace Output{


void WriteCoord(std::shared_ptr<ParticleGroup> pg,
                Box box,
                std::ofstream& out);


void WriteSP(std::shared_ptr<ParticleGroup> pg,
             Box box,
             std::ofstream& out);


void WriteSPO(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out);


void WriteSPF(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out);


void WriteVelocity(std::shared_ptr<ParticleGroup> pg,
                   std::ofstream& out);


void WriteXYZ(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out);


void WriteITPV(std::shared_ptr<ParticleGroup> pg,
               Box box,
               std::ofstream& out);


void WriteITPD(std::shared_ptr<ParticleGroup> pg,
               Box box,
               std::ofstream& out);


void WritePDB(std::shared_ptr<ParticleGroup> pg,
              Box box,
              int frame,
              std::ofstream& out);


void WritePDB(std::shared_ptr<ParticleGroup> pg,
              Box box,
              std::ofstream& out);


void WriteDCDheader(std::shared_ptr<ParticleGroup> pg,
                    int start,
                    int interval,
                    std::ofstream& out);

void WriteDCD(std::shared_ptr<ParticleGroup> pg,
              Box box,
              int frame,int step,
              std::ofstream& out);

void WriteLAMMPS(std::shared_ptr<ParticleGroup> pg,
                 Box box,
                 real t,
                 std::ofstream& out);

void WriteMagnetization(std::shared_ptr<ParticleGroup> pg,
		                std::ofstream& out);


void WriteXYZMagnetization(std::shared_ptr<ParticleGroup> pg,
  	                       std::ofstream& out);


void WriteSPM(std::shared_ptr<ParticleGroup> pg,
  	          Box box,
  	          std::ofstream& out);

//Writes the position and the orientation of each particle in order to be visualized
//with SVV3D (similar to spunto but includes also arrows.)
void WriteSVV(std::shared_ptr<ParticleGroup> pg,
	          Box box,
	          std::ofstream& out);

//Writes the position and the magnetization of each particle in order to be visualized
//with SVV3D (similar to spunto but includes also arrows.)
void WriteSVVM(std::shared_ptr<ParticleGroup> pg,
	           Box box,
	           std::ofstream& out);

//Writes the position, the magnetization and the anisotropy axis of each particle
//in order to be visualized
//with SVV3D (similar to spunto but includes also arrows.)
void WriteSVVMA(std::shared_ptr<ParticleGroup> pg,
	            Box box,
	            std::ofstream& out);


}}}
