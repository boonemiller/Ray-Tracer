//
//  Ray.cpp
//  RayTracer
//
//  Created by Bo Miller on 1/2/19.
//  Copyright © 2019 Bo Miller. All rights reserved.
//
#include "glm/glm/glm.hpp"
#include "glm/glm/gtx/io.hpp"
#include <iostream>
#include "Ray.hpp"
#include <vector>
#include <math.h>
#define PI 3.14159265359
#include <pthread.h>
#include <chrono>
#include "bvh.hpp"
#include <random>

float RAY_EPSILON = 0.000000001;
bool antialiasing = true;
int numBounces = 5;
int numThreads = 4;
int SampPerPix = 4;
float dummyt;
int root = 0;

struct arguments
{
    float chunkSize;
    float width;
    float height;
    float pix[360][720][3];
    glm::vec3 cameraPosition;
    glm::vec3 cameraDirection;
    std::vector<SceneObject>* scene;
    std::vector<Light>* lights;
    int threadNum;
    long time;
    std::vector<Node>* tree;
};


void startRayTracing(float width, float height, float (&pixelcolorBuffer)[360][720][3],glm::vec3 cameraPosition, glm::vec3 cameraDirection, std::vector<SceneObject>& scene, std::vector<Light>& lights)
{
    //*****Builds bvh tree*****
    std::vector<Node> nodes;
    Node start;
    start.maxX = std::numeric_limits<float>::min();
    start.minX = std::numeric_limits<float>::max();
    start.maxY = std::numeric_limits<float>::min();
    start.minY = std::numeric_limits<float>::max();;
    start.maxZ = std::numeric_limits<float>::min();
    start.minZ = std::numeric_limits<float>::max();;
    
    //creates world bounding box information, includes all the objects
    for(SceneObject obj: scene)
    {
        if(obj.position[0]-obj.radius < start.minX)
            start.minX = obj.position[0]-obj.radius;
        if(obj.position[1]-obj.radius < start.minY)
            start.minY = obj.position[1]-obj.radius;
        if(obj.position[2]-obj.radius < start.minZ)
            start.minZ = obj.position[2]-obj.radius;
        
        if(obj.position[0]+obj.radius > start.maxX)
            start.maxX = obj.position[0]+obj.radius;
        if(obj.position[1]+obj.radius > start.maxY)
            start.maxY = obj.position[1]+obj.radius;
        if(obj.position[2]+obj.radius > start.maxZ)
            start.maxZ = obj.position[2]+obj.radius;
    }
    
    if(start.maxX-start.minX > start.maxY-start.minY)
    {
        if(start.maxX-start.minX > start.maxZ - start.minZ)
        {
            start.longestAxis = 0;
            start.midpoint = (start.maxX+start.minX)/2;
        }
    }
    if(start.maxY-start.minY > start.maxX-start.minX)
    {
        if(start.maxY-start.minY > start.maxZ-start.minZ)
        {
            start.longestAxis = 1;
            start.midpoint = (start.maxY+start.minY)/2;
        }
    }
    if(start.maxZ-start.minZ > start.maxX-start.minX)
    {
        if(start.maxZ-start.minZ > start.maxY-start.minY)
        {
            start.longestAxis = 2;
            start.midpoint = (start.maxZ+start.minZ)/2;
        }
    }
    root = constructTree(scene, start, nodes);
    
    //*************************
    pthread_t thread[numThreads];
    int chunkSize = width/numThreads;
    struct arguments *arr[numThreads];
    for(int n = 0; n < numThreads; n++)
    {
        struct arguments *newThread = (struct arguments *)malloc(sizeof(struct arguments));
        newThread->width = width;
        newThread->height = height;
        newThread->cameraPosition = cameraPosition;
        newThread->cameraDirection = cameraDirection;
        newThread->scene = &scene;
        newThread->lights = &lights;
        newThread->chunkSize = chunkSize;
        newThread->threadNum = n;
        newThread->tree = &nodes;
        pthread_create(&thread[n], NULL, &CastRays,(void *)newThread);
        arr[n] = newThread;
        
    }
    
    for(int n = 0; n < numThreads; n++)
    {
       pthread_join(thread[n], NULL);
    }
    
    //fills image buffer
    for(int n = 0; n<numThreads;n++)
    {
        for(int i = arr[n]->threadNum*chunkSize; i<(arr[n]->threadNum+1)*chunkSize;i++)
        {
            for(int j = 0; j<height; j++)
            {
                pixelcolorBuffer[j][i][0] = arr[n]->pix[j][i][0];
                pixelcolorBuffer[j][i][1] = arr[n]->pix[j][i][1];
                pixelcolorBuffer[j][i][2] = arr[n]->pix[j][i][2];
            }
        }
    }
    for(int n = 0; n<numThreads;n++)
        free(arr[n]);
}



void* CastRays(void *arguments)
{
    struct arguments* args = (struct arguments *)arguments;
    
    glm::vec3 n = glm::normalize(args->cameraPosition-args->cameraDirection);
    glm::vec3 u = glm::normalize(glm::cross(glm::vec3(0,1,0),n));
    glm::vec3 v = glm::cross(n,u);
    float fov = 45/(180.0 / PI);
    float d = (args->height/tan(fov/2))/2;
    glm::vec3 L = (args->cameraPosition-n*d) - u * (args->width/2) - v*(args->height/2);
    
    for(int i = args->chunkSize*args->threadNum; i < args->chunkSize*(args->threadNum + 1.0); i++)
    {
        for(int j = 0; j<args->height; j++)
        {
        
            glm::vec3 rayPosition = args->cameraPosition;
            glm::vec3 color = glm::vec3(0.0,0.0,0.0);
            
            if(antialiasing)
            {
                
                glm::vec3 pix = (L+u*float(i)+v*float(j));
                for(int samp = 0; samp<SampPerPix;samp++)
                {
                    float xoffset = ((float) rand()) / (float) RAND_MAX;
                    float yoffset = ((float) rand()) / (float) RAND_MAX;
                    glm::vec3 sample = glm::normalize(glm::vec3(pix[0]+xoffset,pix[1]+yoffset,pix[2])-args->cameraPosition);
                    glm::vec3 colorTest = glm::vec3(0,0,0);
                    if(!intersectObjects(rayPosition, sample, *args->scene, *args->lights, colorTest, false, numBounces, dummyt,*args->tree))
                        rayPlaneIntersection(rayPosition, sample, colorTest,*args->scene,*args->lights,0,*args->tree);
                    color+= colorTest;
                }
                color/=float(SampPerPix);
            }
            else
            {
                glm::vec3 rayDirection = glm::normalize((L+u*float(i)+v*float(j))-args->cameraPosition);
                
                
                if(!intersectObjects(rayPosition, rayDirection, *args->scene, *args->lights, color, false, numBounces, dummyt,*args->tree))
                    rayPlaneIntersection(rayPosition, rayDirection, color,*args->scene,*args->lights,0,*args->tree);
            }
            color = glm::clamp(color,glm::vec3(0,0,0),glm::vec3(1,1,1));
            args->pix[359-j][i][0] = color[0];
            args->pix[359-j][i][1] = color[1];
            args->pix[359-j][i][2] = color[2];
            
        }
    }
    pthread_exit(NULL);
    return NULL;
}

//returns the indecies of the objects to check based on ray
//I think the problem is from overlapping bounding boxes and im not checking both boxes, maybe I need to return a list of possible bounding boxes
void bvhTraverse(glm::vec3 position, glm::vec3 direction, std::vector<Node>& tree, int currentNode,std::vector<int>& boxes)
{
    int retValue = -1;
    if(tree[currentNode].isleaf)
    {
        if(boundingBoxIntersection(position, direction, tree[currentNode]))
        {
            retValue = currentNode;
            boxes.push_back(retValue);
        }
    }
    else
    {
        if(boundingBoxIntersection(position, direction, tree[tree[currentNode].left]))
            bvhTraverse(position, direction, tree, tree[currentNode].left,boxes);
        if(boundingBoxIntersection(position, direction, tree[tree[currentNode].right]))
            bvhTraverse(position, direction, tree, tree[currentNode].right,boxes);
    }
}
//function taken and adapted from
//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
bool boundingBoxIntersection(glm::vec3 position, glm::vec3 direction, Node& box)
{
    float tmin = (box.minX-position[0])/direction[0];
    float tmax = (box.maxX-position[0])/direction[0];
    
    if(tmin>tmax)
        std::swap(tmin,tmax);
    
    float tymin = (box.minY-position[1])/direction[1];
    float tymax = (box.maxY-position[1])/direction[1];
    
    if(tymin>tymax)
        std::swap(tymin,tymax);
    
    if((tmin > tymax) || (tymin > tmax))
        return false;
    
    if (tymin > tmin)
        tmin = tymin;
    
    if (tymax < tmax)
        tmax = tymax;
    
    float tzmin = (box.minZ-position[2])/direction[2];
    float tzmax = (box.maxZ-position[2])/direction[2];
    
    if (tzmin > tzmax)
        std::swap(tzmin, tzmax);
    
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    
    if (tzmin > tmin)
        tmin = tzmin;
    
    if (tzmax < tmax)
        tmax = tzmax;
    
    return true;
}

bool intersectObjects(glm::vec3 position, glm::vec3 direction, std::vector<SceneObject>& scene, std::vector<Light>& lights, glm::vec3& color, bool shadowTest,int numBounces, float& t,std::vector<Node>& tree)
{
    /*Set up for future additions
    
     1. Can add in other types of objects such as triangles
    
     */
    float minT = std::numeric_limits<float>::max();
    SceneObject intersectObj;
    glm::vec3 minTnormal;
    glm::vec3 minTintersection;
    bool intersect = false;
    std::vector<int> boundingBoxes;
    bvhTraverse(position, direction, tree,root,boundingBoxes);
    if(boundingBoxes.size() == 0)
    {
        color = glm::vec3(0,0,0);
        return false;
    }
    
    for(int box: boundingBoxes)
    {
        for(int i = 0; i<tree[box].numObjs; i++)
        {
            if(scene[tree[box].objs[i]].sphere)
            {
                float intersectT;
                glm::vec3 normal;
                glm::vec3 intersection;
                if(intersectSphere(position, direction, scene[tree[box].objs[i]], intersectT,normal,intersection))
                {
                    if(intersectT<minT)
                    {
                        minTnormal = normal;
                        minTintersection = intersection;
                        intersectObj = scene[tree[box].objs[i]];
                        minT = intersectT;
                        intersect = true;
                    }
                }
            }
        }
    }
    
    if(intersect)
    {
        minTintersection = minTintersection+minTnormal*RAY_EPSILON;
        
        if(shadowTest)
        {
            t = minT;
            return true;
        }
        
        color += checkLights(position,direction, minTnormal, minTintersection, lights, minT, intersectObj,scene,tree) + 0.2f * intersectObj.ambient;
        
        if(numBounces != 0)
        {
            glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, minTnormal));
            
            glm::vec3 reflectColor = glm::vec3(0.0,0.0,0.0);;
            float dummyt;
            if(!intersectObjects(minTintersection, reflectedRay, scene, lights, reflectColor, false, numBounces-1,dummyt,tree))
            {
                rayPlaneIntersection(minTintersection, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                color += glm::min(intersectObj.reflective * reflectColor,.2f * reflectColor);
            }
            else
                color += glm::min(intersectObj.reflective * reflectColor,.2f * reflectColor);
        }
        return true;
    }
    color = glm::vec3(0,0,0);
    return false;
}

glm::vec3 checkLights(glm::vec3 position, glm::vec3 direction, glm::vec3 normal, glm::vec3 intersection, std::vector<Light>& lights, float t, SceneObject s, std::vector<SceneObject>& scene,std::vector<Node>& tree)
{
    glm::vec3 color;
    direction = glm::normalize(direction);
    float distance;
    glm::vec3 toLight;
    glm::vec3 reflectFromLight;
    
    for(Light l : lights){
        if(l.area)
        {
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, l.radius);
            glm::vec3 avgcolor;
            int lightSamples = 8;
            for(int i = 0 ;i<lightSamples;i++)
            {
                float r = dis(gen);
                float theta = dis(gen);
                //printf("r = %f theta = %f",randomsamples[i],randomsamples[i+1]);
                
                float x = r * cos(theta);
                float z = r * sin(theta);
                l.position[0] += x;
                l.position[2] += z;
                toLight = glm::normalize(l.position-intersection);
                reflectFromLight = -toLight;
                glm::vec3 dummyC;
                distance = 1.0f;
                float t;
                if(intersectObjects(intersection, toLight, scene, lights, dummyC, true, 0,t,tree))
                {
                    glm::vec3 ipoint = intersection+t*toLight;
                    float dtoLight = sqrt(pow(intersection[0]-l.position[0],2)+pow(intersection[1]-l.position[1],2)+pow(intersection[2]-l.position[2],2));
                    float dtoLightIntersection = sqrt(pow(ipoint[0]-intersection[0],2)+pow(ipoint[1]-intersection[1],2)+pow(ipoint[2]-intersection[2],2));
                    if(dtoLight>dtoLightIntersection)
                        distance = distance * 0;
                }
                
                avgcolor += distance * l.color * ( .6f * s.diffuse * glm::max(glm::dot(toLight,normal),0.0f) + .2f * s.specular * glm::pow(glm::dot(glm::reflect(reflectFromLight, normal), -direction),s.shininess));
            }
            color += avgcolor/(float)lightSamples;
            
            
        }
        else if(l.point)
        {
            float d = sqrt(pow(intersection[0]-l.position[0],2)+pow(intersection[1]-l.position[1],2)+pow(intersection[2]-l.position[2],2));
            distance = 1.0f/(l.constantTerm + l.linearTerm * d + l.quadraticTerm * pow(d,2));
            
            if(distance>1.5)
                distance = .5;
            toLight = glm::normalize(l.position-intersection);
            reflectFromLight = -toLight;
            glm::vec3 dummyC;
            float t;
            if(intersectObjects(intersection, toLight, scene, lights, dummyC, true, 0,t,tree))
            {
                glm::vec3 ipoint = intersection+t*toLight;
                float dtoLight = sqrt(pow(intersection[0]-l.position[0],2)+pow(intersection[1]-l.position[1],2)+pow(intersection[2]-l.position[2],2));
                float dtoLightIntersection = sqrt(pow(ipoint[0]-intersection[0],2)+pow(ipoint[1]-intersection[1],2)+pow(ipoint[2]-intersection[2],2));
                if(dtoLight>dtoLightIntersection)
                    distance = distance * 0;
            }
            color += distance * l.color * ( .6f * s.diffuse * glm::max(glm::dot(toLight,normal),0.0f) + .2f * s.specular * glm::pow(glm::dot(glm::reflect(reflectFromLight, normal), -direction),s.shininess));
        }
        else
        {
            distance = 1.0f;
            toLight = -glm::normalize(l.direction);
            reflectFromLight = glm::normalize(l.direction);
            glm::vec3 dummyC;
            float dummyT;
            if(intersectObjects(intersection, toLight, scene, lights, dummyC, true, 0,dummyT,tree))
            {
                distance = distance * 0;
            }
            color += distance * l.color * ( .6f * s.diffuse * glm::max(glm::dot(toLight,normal),0.0f) + .2f * s.specular * glm::pow(glm::dot(glm::reflect(reflectFromLight, normal), -direction),s.shininess));
        }
        
        
    }
    return color;
}

/*Checks intersection of a ray with a sphere*/
bool intersectSphere(glm::vec3 position, glm::vec3 direction, SceneObject s, float& iTime,glm::vec3& normal, glm::vec3& intersection)
{

    float a = glm::dot(direction, direction);
    float b = 2 * glm::dot(direction,position-s.center);
    float c = glm::dot(s.center,s.center) + glm::dot(position,position) + (-2 * glm::dot(s.center,position)) - pow(s.radius,2);
 
    float discriminant = b*b - 4*a*c;
 
    if(discriminant > 0.0+RAY_EPSILON)
    {
 
        float t = (-b - sqrt(discriminant))/(2*a);
 
        float t2 = (-b + sqrt(discriminant))/(2*a);
 
 
        if(t2>RAY_EPSILON)
        {
            //we know we have some intersection
 
            if( t > RAY_EPSILON )
            {
                iTime = t;
            }
            else
            {
                iTime = t2;
            }
            
            intersection = position+t*direction;
            normal = glm::normalize((intersection-s.center)/s.radius);
            return true;
        }
    }
    return false;
}


/*Used to get a ray-wall intersection point*/
bool rayPlaneIntersection(glm::vec3 position, glm::vec3 direction, glm::vec3& color, std::vector<SceneObject>& scene, std::vector<Light>& lights, int numBounces,std::vector<Node>& tree)
{
    //floor
    glm::vec3 up = glm::vec3(0,1,0);
    float denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(0,-2,0)-position),up)/denom;
        if(t >= 0.0-.0001f)
        {
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(1.0, 0.2, 0.2);
            wall.diffuse = glm::vec3(1.0, 0.2, 0.2);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            color = wall.ambient * 0.2f + checkLights(position, direction, up, intersect, lights, t, wall, scene, tree);
            
            if(intersect[2]>-15 && intersect[2]<13 && intersect[0]>-6 && intersect[0] < 6)
            {
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt, tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                    
                }
                return true;
            }
        }
    }
    
    //left wall
    up = glm::vec3(1,0,0);
    denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(-6,0,0)-position),up)/denom;
        if(t >= 0.0-.0001f)
        {
            
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(0.2, 0.2, 1.0);
            wall.diffuse = glm::vec3(0.2, 0.2, 1.0);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            color = wall.ambient * 0.2f+checkLights(position, direction, up, intersect, lights, t, wall, scene,tree);
            if(intersect[2]>-15 && intersect[2] < 13 && intersect[1] < 9 && intersect[1]>-2)
            {
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt,tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                }
                return true;
            }
        }
    }
    
    //front wall, green wall in front of camera
    up = glm::vec3(0,0,1);
    denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(0,0,-15)-position),up)/denom;
        if(t >= 0.0-.0001f)
        {
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(0.0, 1.0, 0.0);
            wall.diffuse = glm::vec3(0.0, 1.0, 0.0);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            color = wall.ambient *.2f + checkLights(position, direction, up, intersect, lights, t, wall, scene,tree);
            
            if(intersect[0] < 6 && intersect[0] > -6 && intersect[1] < 9 && intersect[1] > -2 )
            {
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt,tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                }
                return true;
            }
        }
    }
    
    //back wall, yellow wall behind camera
    up = glm::vec3(0,0,-1);
    denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(0,0,13)-position),up)/denom;
        if(t >= 0.0-.0001f)
        {
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(1.0, 1.0, 0.0);
            wall.diffuse = glm::vec3(1.0, 1.0, 0.0);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            color = wall.ambient *.2f + checkLights(position, direction, up, intersect, lights, t, wall, scene,tree);
            
            if(intersect[0] < 6 && intersect[0] > -6 && intersect[1] < 9 && intersect[1] > -2 )
            {
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt,tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                }
                return true;
            }
        }
    }
    
    //right wall
    up = glm::vec3(-1,0,0);
    denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(6,0,0)-position),up)/denom;
        if(t >= 0.0-.0001f)
        {
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(0.2, 0.2, 1.0);
            wall.diffuse = glm::vec3(0.2, 0.2, 1.0);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            color = wall.ambient *.2f + checkLights(position, direction, up, intersect, lights, t, wall, scene,tree);
            
            if(intersect[2]>-15 && intersect[2] < 13 && intersect[1] < 9 && intersect[1]>-2)
            {
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt,tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                }
                return true;
            }
        }
    }
    
    //ceiling
    up = glm::vec3(0,-1,0);
    denom = glm::dot(up,direction);
    if(abs(denom) > .0001f)
    {
        float t = glm::dot((glm::vec3(0,9,0)-position), up)/denom;
        if(t >= 0.0-.0001f)
        {
            glm::vec3 intersect = position+t*direction;
            
            SceneObject wall;
            wall.ambient = glm::vec3(.5, .5, .5);
            wall.diffuse = glm::vec3(.9, .9, .9);
            wall.specular = glm::vec3(0.0,0.0,0.0);
            wall.shininess = 2;
            wall.reflective = glm::vec3(0.0,0.0,0.0);
            color = wall.ambient * .2f + checkLights(position, direction, up, intersect, lights, t, wall, scene,tree);
        
            if(intersect[2]>-15 && intersect[2]<13 && intersect[0]>-6 && intersect[0] < 6)
            {
                
                if(numBounces != 0)
                {
                    glm::vec3 reflectedRay = glm::normalize(glm::reflect(direction, up));
                    
                    glm::vec3 reflectColor = glm::vec3(0,0,0);
                    if(intersectObjects(intersect, reflectedRay, scene, lights, reflectColor, false, numBounces-1, dummyt,tree))
                        color += .2f*reflectColor;
                    else
                        rayPlaneIntersection(intersect, reflectedRay, reflectColor, scene, lights, numBounces-1,tree);
                }
                return true;
            }
        }
    }
    color = glm::vec3(0,0,0);
    return false;
}

