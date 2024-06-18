#pragma once
#ifndef __SCATTERED_HPP__
#define __SCATTERED_HPP__

#include "Ray.hpp"

namespace ConductorPathTracer
{
    struct Scattered
    {
        Ray ray = {};
        Vec3 attenuation = {};
        Vec3 emitted = {};
        float pdf = {0.f};
        Ray refraction = {};
        Vec3 refractionRate = {};
    };
    
}

#endif