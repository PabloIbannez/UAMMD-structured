#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"

#include <vector>
#include <string>
#include <map>
#include <algorithm>

namespace uammd{
namespace structured{

void stateLoader(ParticleData* pd,DataEntry& data);
void updateState(ParticleData* pd,DataEntry& data);

}}
