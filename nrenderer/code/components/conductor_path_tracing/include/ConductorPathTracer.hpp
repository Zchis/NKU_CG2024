#pragma once
#ifndef __CONDUCTOR_PATH_TRACER_HPP__
#define __CONDUCTOR_PATH_TRACER_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"
#include "BVH.hpp"
#include "shaders/ShaderCreator.hpp"

#include <tuple>
namespace ConductorPathTracer
{
    using namespace NRenderer;
    using namespace std;

    class ConductorPathTracerRenderer
    {
    public:
    private:
        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;

        using SCam = ConductorPathTracer::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        ConductorPathTracerRenderer(SharedScene spScene)
            : spScene               (spScene)
            , scene                 (*spScene)
            , camera                (spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;
        }
        ~ConductorPathTracerRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);
        BVHTree* bvhtree;
        vector<AABB> aabb;
        void buildAABB() {
            BVHTree* bvhtree = new BVHTree(this->spScene);
            BVHNode* n = new BVHNode();
            vector<Node> objects = spScene->nodes;
            n->buildAABB(spScene, objects);
            aabb = n->bounds;
        }
    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);
    };
}

#endif