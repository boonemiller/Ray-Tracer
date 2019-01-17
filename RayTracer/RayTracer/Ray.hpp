//
//  Ray.hpp
//  RayTracer
//
//  Created by Bo Miller on 1/2/19.
//  Copyright © 2019 Bo Miller. All rights reserved.
//

#include "SceneObjects.hpp"

#ifndef Ray_hpp
#define Ray_hpp

#include <stdio.h>
void startRayTracing(float width, float height, float (&pixelcolorBuffer)[360][720][3],glm::vec3 cameraPosition, glm::vec3 cameraDirection, std::vector<SceneObject>& scene, std::vector<Light>& lights);
void* CastRays(void *arguments);
//bool intersectSphere(glm::vec3 position, glm::vec3 direction, std::vector<SceneObject>& scene, std::vector<Light>& lights, glm::vec3& color, bool shadowTest,int numBounces,glm::vec3& interesectPoint);
bool intersectSphere(glm::vec3 position, glm::vec3 direction, SceneObject s, float& iTime,glm::vec3& normal, glm::vec3& intersection);
glm::vec3 checkLights(glm::vec3 position, glm::vec3 direction, glm::vec3 normal, glm::vec3 intersection, std::vector<Light>& lights, float t, SceneObject s, std::vector<SceneObject>& scene);
bool rayPlaneIntersection(glm::vec3 position, glm::vec3 direction, glm::vec3& color, std::vector<SceneObject>& scene, std::vector<Light>& lights, int numBounces);
bool intersectTriangle(glm::vec3 v1,glm::vec3 v2,glm::vec3 v3,glm::vec3 position, glm::vec3 direction, glm::vec3& n, glm::vec3& intersection, float& time);
bool intersectObjects(glm::vec3 position, glm::vec3 direction, std::vector<SceneObject>& scene, std::vector<Light>& lights, glm::vec3& color, bool shadowTest,int numBounces, float& t);
#endif /* Ray_hpp */
