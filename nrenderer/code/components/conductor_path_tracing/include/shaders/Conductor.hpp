#pragma once

#ifndef __CONDUCTOR_HPP__
#define __CONDUCTOR_HPP__

#include "Shader.hpp"

namespace ConductorPathTracer
{
    class Conductor : public Shader
    {
        //ȫ����
    private:
        Vec3 reflect; //����

    public:
        Conductor(Material& material, vector<Texture>& textures);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}

#endif