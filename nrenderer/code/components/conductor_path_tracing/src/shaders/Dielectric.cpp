#include "server/Server.hpp"
#include "shaders/Dielectric.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace ConductorPathTracer {
	Dielectric::Dielectric(Material& material, vector<Texture>& textures) :Shader(material, textures) {
		auto absorbed = material.getProperty<Property::Wrapper::RGBType>("absorbed");
		if (absorbed) this->absorbed = (*absorbed).value;
		else this->absorbed = { 1,1,1 };
		auto ior = material.getProperty<Property::Wrapper::FloatType>("ior");
		if (ior) this->ior = (*ior).value;
		else this->ior = 1.5f;

	}
	Scattered Dielectric::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal)const {
		Vec3 begin = hitPoint;
		Vec3 V = glm::normalize(ray.direction);//计算入射方向
		Vec3 N = glm::normalize(normal);
		Vec3 R = glm::reflect(V, N);//反射
		float trueIor = 1.0f / ior;
		//计算导体的菲涅尔系数
		float cosTheta = glm::dot(N, V);
		
		float R0 = (1 - ior) / (1 + ior);
		R0 = R0 * R0;   //计算反射率
		float F0=1; // 
		F0 -=  R0; // 
		F0 *= pow(1 - abs(cosTheta), 5);
		F0 += R0;  //反射率

		Vec3 Reflex = F0 * absorbed;//反射率和自身颜色
		Vec3 Refraction = (1-F0) * absorbed;//除去反射和自身颜色
		if (cosTheta > 0) {//里面射到外面,所以是锐角,法向量朝向外
			N = -N;
			trueIor = ior;
		}
		
		Vec3 T; // 折射光线

			// 光线从外部射入
		T = glm::refract(V, N, trueIor);


		if (!T.x && !T.y && !T.z) { // 折射光线是零向量
			Refraction = Vec3(0, 0, 0); // 没有折射光线
		}
		Scattered s;
		s.ray = Ray(begin, R);
		s.attenuation = Reflex;
		s.refraction = Ray{ begin, T };
		s.refractionRate = Refraction;
		
		//getServer().logger.log("111111111");
		return s;
	}
}