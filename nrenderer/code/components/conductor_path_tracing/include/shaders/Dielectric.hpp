#pragma once
#ifndef __DIELECTRIC_HPP__
#define __DIELECTRIC_HPP__

#include "Shader.hpp"

namespace ConductorPathTracer
{
    class Dielectric : public Shader
    {
        //绝缘体材质
    private:
        Vec3 absorbed; //金属原色
        float ior;//折射率

    public:
        Dielectric(Material& material, vector<Texture>& textures);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}

#endif