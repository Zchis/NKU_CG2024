#include "server/Server.hpp"

#include "ConductorPathTracer.hpp"
#include "BVH.hpp"
#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"

namespace ConductorPathTracer
{
    RGB ConductorPathTracerRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void ConductorPathTracerRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        for(int i=off; i<height; i+=step) {
            for (int j=0; j<width; j++) {
                Vec3 color{0, 0, 0};
                for (int k=0; k < samples; k++) {
                    auto r = defaultSamplerInstance<UniformInSquare>().sample2d();
                    float rx = r.x;
                    float ry = r.y;
                    float x = (float(j)+rx)/float(width);
                    float y = (float(i)+ry)/float(height);
                    auto ray = camera.shoot(x, y);
                    color += trace(ray, 0);
                }
                color /= samples;
                color = gamma(color);
                pixels[(height-i-1)*width+j] = {color, 1};
            }
        }
    }

    auto ConductorPathTracerRenderer::render() -> RenderResult {
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(m, scene.textures));
        }

        RGBA* pixels = new RGBA[width*height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);
        buildAABB();
        bvhtree = new BVHTree(spScene);
        bvhtree->root = bvhtree->buildBVH(aabb);
        const auto taskNums = 24;
        thread t[taskNums];
        for (int i=0; i < taskNums; i++) {
            t[i] = thread(&ConductorPathTracerRenderer::renderTask,
                this, pixels, width, height, i, taskNums);
        }
        for(int i=0; i < taskNums; i++) {
            t[i].join();
        }
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void ConductorPathTracerRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord ConductorPathTracerRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
        HitRecord hitRecord = bvhtree->treeIntersect(r, closest);
        if (hitRecord && hitRecord->t < closest) {
            closest = hitRecord->t;
            closestHit = hitRecord;
        }
        /*for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }*/
        return closestHit; 
    }
    
    tuple<float, Vec3> ConductorPathTracerRenderer::closestHitLight(const Ray& r) {
        Vec3 v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v };
    }

    RGB ConductorPathTracerRenderer::trace(const Ray& r, int currDepth) {
        if (currDepth == depth) return scene.ambient.constant;
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);
        // hit object
        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            if (spScene->materials[mtlHandle.index()].type ==0) {
                auto scatteredRay = scattered.ray;
                auto attenuation = scattered.attenuation;
                auto emitted = scattered.emitted;
                auto next = trace(scatteredRay, currDepth + 1);
                float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
                float pdf = scattered.pdf;
                /**
                 * emitted      - Le(p, w_0)
                 * next         - Li(p, w_i)
                 * n_dot_in     - cos<n, w_i>
                 * atteunation  - BRDF
                 * pdf          - p(w)
                 **/
                return emitted + attenuation * next * n_dot_in / pdf;
            }
            else if (spScene->materials[mtlHandle.index()].type == 2) {
                auto r = scattered.ray;
                auto R = scattered.attenuation;
                auto t = scattered.refraction;
                auto T = scattered.refractionRate;
                auto rRGB = R == Vec3(0.f) ? RGB(0.f) : trace(r, currDepth + 1); 
                auto tRGB = trace(t, currDepth + 1); 
                return rRGB * R + tRGB * T;
            }
            else if (spScene->materials[mtlHandle.index()].type == 3) {
                auto L = scattered.ray;
                auto attenuation = scattered.attenuation;
                auto next = trace(L, currDepth + 1);
                //getServer().logger.log("111111111");
                return attenuation * next;
            }
        }
        // 
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{0};
        }
    }
}