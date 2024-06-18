#pragma once

#ifndef __BVH_HPP__
#define __BVH_HPP__
#include <vector>
#include <algorithm>
#include "scene/Scene.hpp"
#include "scene/Model.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"
#include "intersections/intersections.hpp"
#include "shaders/ShaderCreator.hpp"

#include <tuple>
namespace ConductorPathTracer
{
	using namespace NRenderer;
	using namespace std;
	class AABB;
	class AABB {
	public:
		Vec3 min; // Minimum corner of the bounding box
		Vec3 max; // Maximum corner of the bounding box
		enum class Type
		{
			SPHERE,
			TRIANGLE,
			PLANE,
			//MESH,
			EMPTY
		};
		Sphere* sphere = nullptr;
		Triangle* triangle = nullptr;
		Plane* plane = nullptr;
		//Mesh* mesh;

		Type type;
		AABB() : min(-1000000, -1000000, -1000000), max(1000000, 1000000, 1000000) { type = Type::EMPTY; }
		AABB(Sphere* sp) {
			this->sphere = sp;
			type = Type::SPHERE;
			//根据球创建Box
			Vec3 nodesp = sp->position;
			nodesp -= sp->radius;//中心减去半径
			min = nodesp;
			nodesp = sp->position;
			nodesp += sp->radius;
			max = nodesp;
		}
		AABB(Triangle* tr) {
			this->triangle = tr;
			type = Type::TRIANGLE;
			float x_min = std::min(std::min(tr->v1.x, tr->v2.x), tr->v3.x);
			float y_min = std::min(std::min(tr->v1.y, tr->v2.y), tr->v3.y);
			float z_min = std::min(std::min(tr->v1.z, tr->v2.z), tr->v3.z);
			float x_max = std::max(std::max(tr->v1.x, tr->v2.x), tr->v3.x);
			float y_max = std::max(std::max(tr->v1.y, tr->v2.y), tr->v3.y);
			float z_max = std::max(std::max(tr->v1.z, tr->v2.z), tr->v3.z);
			min = Vec3(x_min, y_min, z_min);
			max = Vec3(x_max, y_max, z_max);
		}
		AABB(Plane* pl) {
			this->plane = pl;
			type = Type::PLANE;
			Vec3 pos1 = pl->position;
			Vec3 u = pl->u;
			Vec3 v = pl->v;
			Vec3 pos2 = pos1 + u;
			Vec3 pos3 = pos1 + v;
			Vec3 pos4 = pos1 + v + u;
			Vec3 normal = pl->normal;//法向量 

			pos1 -= normal * 0.01f;
			pos2 -= normal * 0.01f;
			pos3 += normal * 0.01f;
			pos4 += normal * 0.01f;
			float x_min = std::min(pos1.x, std::min(std::min(pos2.x, pos3.x), pos4.x));
			float y_min = std::min(pos1.y, std::min(std::min(pos2.y, pos3.y), pos4.y));
			float z_min = std::min(pos1.z, std::min(std::min(pos2.z, pos3.z), pos4.z));
			float x_max = std::max(pos1.x, std::max(std::max(pos2.x, pos3.x), pos4.x));
			float y_max = std::max(pos1.y, std::max(std::max(pos2.y, pos3.y), pos4.y));
			float z_max = std::max(pos1.z, std::max(std::max(pos2.z, pos3.z), pos4.z));

			min = Vec3(x_min, y_min, z_min);
			max = Vec3(x_max, y_max, z_max);

		}
		//AABB(Vec3 v_1, Vec3 v_2, Vec3 v_3, Handle mat);
		inline bool Intersect(const Ray& ray);
		Vec3 centroid()const {
			return Vec3(min.x + (max.x - min.x) * 0.5f,
				min.y + (max.y - min.y) * 0.5f,
				min.z + (max.z - min.z) * 0.5f);

		}

	};
	class BVHNode {
	public:
		AABB box; // 节点的包围盒
		BVHNode* left; // 左子节点
		BVHNode* right; // 右子节点
		//Entity object;
		vector <AABB> bounds;
		BVHNode() : left(nullptr), right(nullptr) {
			box = AABB();
		}

		BVHNode(AABB b, BVHNode* l, BVHNode* r) : box(b), left(l), right(r) {}

		/*void setObject(SharedScene spscene, Node object) {
			auto node=object
		}*/
		bool isleaf() {
			if (this->left == nullptr && this->right == nullptr)
				return true;
			return false;
		}
		void buildAABB(SharedScene spscene, vector<Node> objects) {
			for (int i = 0; i < objects.size(); i++) {
				auto node = objects[i];
				auto model = spscene->models[node.model];
				if (objects[i].type == Node::Type::SPHERE) {
					//printf("dwa");
					AABB s(&spscene->sphereBuffer[node.entity]);
					bounds.push_back(s);
				}
				if (objects[i].type == Node::Type::TRIANGLE) {

					AABB t(&spscene->triangleBuffer[node.entity]);
					bounds.push_back(t);

				}
				if (objects[i].type == Node::Type::PLANE) {
					//printf("dwdda ");
					AABB p(&spscene->planeBuffer[node.entity]);
					bounds.push_back(p);

				}
			}
		}
		bool Intersect(const Ray& ray) {
			return box.Intersect(ray);
		}
	};
	/*struct AABBcomparex {
		bool operator()(const AABB& b1, const AABB& b2){ return b1.centroid().x < b2.centroid().x; }
	};
	struct AABBcomparey {
		bool operator()(const AABB& b1, const AABB& b2) { return b1.centroid().y < b2.centroid().y; }
	};
	struct AABBcomparez {
		bool operator()(const AABB& b1, const AABB& b2) { return b1.centroid().z < b2.centroid().z; }
	};*/
	class BVHTree {
	public:
		BVHNode* root;
		SharedScene spscene;
		vector<AABB> bounds;
		BVHTree() {};
		BVHTree(SharedScene spScene) {
			spscene = spScene;
			root = nullptr;
		}
		inline AABB Merge(AABB& box1, AABB& box2);
		inline AABB Merge(AABB& box1, Vec3& p);
		inline BVHNode* buildBVH(vector<AABB> bounds) {//建树
			BVHNode* node = new BVHNode;
			node->bounds = bounds;
			if (bounds.size() == 1) {
				node->box = bounds[0];
				node->left = nullptr;
				node->right = nullptr;
				return node;
			}
			if (bounds.size() == 2) {
				node->box = Merge(bounds[1], bounds[0]);
				node->left = buildBVH({ bounds[0] });
				node->right = buildBVH({ bounds[1] });
				return node;
			}
			if (bounds.size() > 2) {
				auto max_x = bounds.begin()->centroid().x;
				auto min_x = bounds.begin()->centroid().x;
				auto max_y = bounds.begin()->centroid().y;
				auto min_y = bounds.begin()->centroid().y;
				auto max_z = bounds.begin()->centroid().z;
				auto min_z = bounds.begin()->centroid().z;
				/* auto max_x = std::max(bounds.begin(), bounds.end(), AABBcomparex());
				 auto min_x = std::min(bounds.begin(), bounds.end(), AABBcomparex() );
				 auto max_y = std::max(bounds.begin(), bounds.end(), AABBcomparey() );
				 auto min_y = std::min(bounds.begin(), bounds.end(), AABBcomparey() );
				 auto max_z = std::max(bounds.begin(), bounds.end(), AABBcomparez() );
				 auto min_z = std::min(bounds.begin(), bounds.end(), AABBcomparez() );*/
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					max_x = std::max(max_x, it->centroid().x);
				}
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					min_x = std::min(min_x, it->centroid().x);
				}
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					max_y = std::max(max_y, it->centroid().y);
				}
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					min_y = std::min(min_y, it->centroid().y);
				}
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					max_z = std::max(max_z, it->centroid().z);
				}
				for (auto it = bounds.begin(); it != bounds.end(); ++it) {
					min_z = std::min(min_z, it->centroid().z);
				}
				float stride_x = max_x - min_x;
				float stride_y = max_y - min_y;
				float stride_z = max_z - min_z;

				if (stride_x == std::max({ stride_x ,stride_y ,stride_z })) {
					std::sort(bounds.begin(), bounds.end(), []
					(const AABB& b1, const AABB& b2) {
							return b1.centroid().x < b2.centroid().x;
						}
					);
				}
				else if (stride_y == std::max({ stride_x ,stride_y ,stride_z })) {
					std::sort(bounds.begin(), bounds.end(), []
					(const AABB& b1, const AABB& b2) {
							return b1.centroid().y < b2.centroid().y;
						}
					);
				}
				else if (stride_z == std::max({ stride_x ,stride_y, stride_z })) {
					std::sort(bounds.begin(), bounds.end(), []
					(const AABB& b1, const AABB& b2) {
							return b1.centroid().z < b2.centroid().z;
						}
					);
				}

				vector<AABB> boundsleft(bounds.begin(), bounds.begin() + (bounds.size() / 2) + 1);
				vector<AABB> boundsright(bounds.begin() + (bounds.size() / 2) + 1, bounds.end());
				node->left = buildBVH(boundsleft);
				node->right = buildBVH(boundsright);
				node->box = Merge(node->left->box, node->right->box);
				return node;
			}
		}
		inline HitRecord Intersect(const Ray& ray, BVHNode* node, float closest) {
			HitRecord Hit;
			if (!node->box.Intersect(ray)) {

				return getMissRecord();
			}
			//printf("ddwdw3333a");
			if (node->isleaf()) {
				//printf("ddwdw3333a");
				//printf("ddwdw3333a");
				if (node->box.type == AABB::Type::SPHERE) {
					//printf("ddwa");
					auto hitRecord = Intersection::xSphere(ray, *node->box.sphere, 0.000001, closest);
					if (hitRecord && hitRecord->t < closest) {
						closest = hitRecord->t;
						Hit = hitRecord;
					}
				}
				else if (node->box.type == AABB::Type::TRIANGLE) {
					auto hitRecord = Intersection::xTriangle(ray, *node->box.triangle, 0.000001, closest);
					if (hitRecord && hitRecord->t < closest) {
						closest = hitRecord->t;
						Hit = hitRecord;
					}
				}
				else if (node->box.type == AABB::Type::PLANE) {
					// printf("  dd    wsd        a  ");
					auto hitRecord = Intersection::xPlane(ray, *node->box.plane, 0.000001, closest);
					if (hitRecord && hitRecord->t < closest) {
						closest = hitRecord->t;
						Hit = hitRecord;
					}
				}
				return Hit;
			}
			else {
				//printf("ddwdw3333a");
				HitRecord hitleft = Intersect(ray, node->left, closest);
				HitRecord hitright = Intersect(ray, node->right, closest);
				if (hitleft == nullopt && hitright == nullopt) {
					return nullopt;
				}
				else if (hitleft != nullopt && hitright == nullopt) {
					return hitleft;
				}
				else if (hitright != nullopt && hitleft == nullopt) {
					return hitright;
				}
				else {
					if (hitleft->t == std::min(hitleft->t, hitright->t))
						return hitleft;
					return hitright;
				}
			}
		}
		inline HitRecord treeIntersect(const Ray& ray, float closest) {
			if (BVHTree::root == nullptr) {
				return getMissRecord();
			}
			//printf("123456");
			root->bounds = bounds;
			HitRecord hitrecord = Intersect(ray, root, closest);
			return hitrecord;
		}
	};
	inline AABB BVHTree::Merge(AABB& box1, AABB& box2) {
		AABB box;
		box.min = Vec3(std::min(box1.min.x, box2.min.x), std::min(box1.min.y, box2.min.y), std::min(box1.min.z, box2.min.z));
		box.max = Vec3(std::max(box1.max.x, box2.max.x), std::max(box1.max.y, box2.max.y), std::max(box1.max.z, box2.max.z));
		return box;
	}
	inline AABB BVHTree::Merge(AABB& box1, Vec3& p)
	{
		AABB box;
		box.min = Vec3(std::min(box1.min.x, p.x), std::min(box1.min.y, p.y), std::min(box1.min.z, p.z));
		box.max = Vec3(std::max(box1.max.x, p.x), std::max(box1.max.y, p.y), std::max(box1.max.z, p.z));
		return box;
	}

	inline bool AABB::Intersect(const Ray& ray) {
		Vec3 pos = ray.origin;
		Vec3 dir = ray.direction;
		dir = 1.0f / dir;
		Vec3 dis_max = max;
		dis_max -= pos;
		Vec3 dis_min = min;
		dis_min -= pos;
		float t_xmin = dis_min.x * dir.x;//时间比
		float t_ymin = dis_min.y * dir.y;//时间比
		float t_zmin = dis_min.z * dir.z;//时间比
		float t_xmax = dis_max.x * dir.x;//时间比
		float t_ymax = dis_max.y * dir.y;//时间比
		float t_zmax = dis_max.z * dir.z;//时间比
		if (t_xmin > t_xmax)std::swap(t_xmin, t_xmax);
		if (t_ymin > t_ymax)std::swap(t_ymin, t_ymax);
		if (t_zmin > t_zmax)std::swap(t_zmin, t_zmax);
		float t_enter = std::max(t_xmin, std::max(t_ymin, t_zmin));
		float t_exit = std::min(t_xmax, std::min(t_ymax, t_zmax));
		if (t_exit < 0) {//说明射线不打向盒子
			return false;
		}
		else if (t_exit >= 0 && t_enter < 0) {//说明在内部
			return true;
		}
		else if (t_exit > t_enter) {
			return true;
		}
		else {
			return false;
		}
	}
}
#endif