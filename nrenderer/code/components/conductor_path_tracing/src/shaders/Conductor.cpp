#include "server/Server.hpp"
#include "shaders/Conductor.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace ConductorPathTracer {
	Conductor::Conductor(Material& material, vector<Texture>& textures) :Shader(material, textures) {
		auto reflect = material.getProperty<Property::Wrapper::RGBType>("reflect");
		if (reflect) this->reflect = (*reflect).value;
		else this->reflect = { 1,1,1 };
	}
	Scattered Conductor::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal)const {
		//首先计算出射方向
		Vec3 V = glm::normalize(ray.direction);//计算入射方向
		Vec3 N = glm::normalize(normal);
		
		// 计算反射光线方向
		Vec3 R = glm::reflect(V, N);
		//计算导体的菲涅尔系数
		float cosTheta = glm::dot(N, V);
		Vec3 F0(1.0f); // 初始化
		F0 = F0 - this->reflect; // 逐元素减法
		F0*= pow(1 - abs(cosTheta), 5);
		F0 += this->reflect ;
		
		Vec3 BRDF = F0;
		BRDF*= abs(glm::dot(N, V));
		

		Vec3 attenuation = BRDF * abs(glm::dot(R, N)) * this->reflect;
		
		Scattered s;
		s.ray = Ray(hitPoint, R);
		s.attenuation = Vec3(1.0f);
		s.attenuation = attenuation;
		//getServer().logger.log("111111111");
		return s;
	}
}