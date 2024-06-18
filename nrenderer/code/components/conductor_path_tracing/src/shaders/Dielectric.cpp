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
		Vec3 V = glm::normalize(ray.direction);//�������䷽��
		Vec3 N = glm::normalize(normal);
		Vec3 R = glm::reflect(V, N);//����
		float trueIor = 1.0f / ior;
		//���㵼��ķ�����ϵ��
		float cosTheta = glm::dot(N, V);
		
		float R0 = (1 - ior) / (1 + ior);
		R0 = R0 * R0;   //���㷴����
		float F0=1; // 
		F0 -=  R0; // 
		F0 *= pow(1 - abs(cosTheta), 5);
		F0 += R0;  //������

		Vec3 Reflex = F0 * absorbed;//�����ʺ�������ɫ
		Vec3 Refraction = (1-F0) * absorbed;//��ȥ�����������ɫ
		if (cosTheta > 0) {//�����䵽����,���������,������������
			N = -N;
			trueIor = ior;
		}
		
		Vec3 T; // �������

			// ���ߴ��ⲿ����
		T = glm::refract(V, N, trueIor);


		if (!T.x && !T.y && !T.z) { // ���������������
			Refraction = Vec3(0, 0, 0); // û���������
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