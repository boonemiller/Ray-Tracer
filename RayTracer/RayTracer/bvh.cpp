//
//  bvh.cpp
//  RayTracer
//
//  Created by Bo Miller on 2/1/19.
//  Copyright © 2019 Bo Miller. All rights reserved.
//
#include "bvh.hpp"
#include <limits>
int constructTree(std::vector<SceneObject>& objects, Node& currentNode, std::vector<Node>& nodes)
{
    if(objects.size() <= 3)
    {
        for(int i = 0; i<objects.size();i++)
        {
            currentNode.objs[i] = objects[i].objNum;
        }
        currentNode.numObjs = objects.size();
        currentNode.isleaf = true;
        nodes.push_back(currentNode);
        return nodes.size()-1;
    }
    
    Node newLeftNode;
    newLeftNode.left = NULL;
    newLeftNode.right = NULL;
    Node newRightNode;
    newRightNode.left = NULL;
    newRightNode.right = NULL;
    std::vector<int> leftPointers;
    std::vector<SceneObject> leftObjects;
    glm::vec3 midLeft;
    float maxLeftX = std::numeric_limits<float>::min();
    float minLeftX = std::numeric_limits<float>::max();
    float maxLeftY = std::numeric_limits<float>::min();
    float minLeftY = std::numeric_limits<float>::max();
    float maxLeftZ = std::numeric_limits<float>::min();
    float minLeftZ = std::numeric_limits<float>::max();
    glm::vec3 midRight;
    float maxRightX = std::numeric_limits<float>::min();
    float minRightX = std::numeric_limits<float>::max();
    float maxRightY = std::numeric_limits<float>::min();
    float minRightY = std::numeric_limits<float>::max();
    float maxRightZ = std::numeric_limits<float>::min();
    float minRightZ = std::numeric_limits<float>::max();
    std::vector<int> rightPointers;
    std::vector<SceneObject> rightObjects;
    for(int i = 0; i < objects.size(); i++)
    {
        
        if(objects[i].position[currentNode.longestAxis] < currentNode.midpoint)
        {
            if(objects[i].position[0]-objects[i].radius < minLeftX)
                minLeftX = objects[i].position[0]-objects[i].radius;
            if(objects[i].position[1]-objects[i].radius < minLeftY)
                minLeftY = objects[i].position[1]-objects[i].radius;
            if(objects[i].position[2]-objects[i].radius < minLeftZ)
                minLeftZ = objects[i].position[2]-objects[i].radius;
            
            if(objects[i].position[0]+objects[i].radius > maxLeftX)
                maxLeftX = objects[i].position[0]+objects[i].radius;
            if(objects[i].position[1]+objects[i].radius > maxLeftY)
                maxLeftY = objects[i].position[1]+objects[i].radius;
            if(objects[i].position[2]+objects[i].radius > maxLeftZ)
                maxLeftZ = objects[i].position[2]+objects[i].radius;
            
            midLeft += objects[i].position;
            leftObjects.push_back(objects[i]);
        }
        else
        {
            if(objects[i].position[0]-objects[i].radius < minRightX)
                minRightX = objects[i].position[0]-objects[i].radius;
            if(objects[i].position[1]-objects[i].radius < minRightY)
                minRightY = objects[i].position[1]-objects[i].radius;
            if(objects[i].position[2]-objects[i].radius < minRightZ)
                minRightZ = objects[i].position[2]-objects[i].radius;
            
            if(objects[i].position[0]+objects[i].radius > maxRightX)
                maxRightX = objects[i].position[0]+objects[i].radius;
            if(objects[i].position[1]+objects[i].radius > maxRightY)
                maxRightY = objects[i].position[1]+objects[i].radius;
            if(objects[i].position[2]+objects[i].radius > maxRightZ)
                maxRightZ = objects[i].position[2]+objects[i].radius;
            midRight += objects[i].position;
            rightObjects.push_back(objects[i]);
        }
    }
    
    midLeft = glm::vec3(midLeft[0]/leftObjects.size(),midLeft[1]/leftObjects.size(),midLeft[2]/leftObjects.size());
    midRight = glm::vec3(midRight[0]/rightObjects.size(),midRight[1]/rightObjects.size(),midRight[2]/rightObjects.size());
    
    if(maxLeftX-minLeftX > maxLeftY - minLeftY)
    {
        if(maxLeftX-minLeftX > maxLeftZ - minLeftZ)
        {
            newLeftNode.longestAxis = 0;
            newLeftNode.midpoint = midLeft[0];
        }
    }
    if(maxLeftY-minLeftY > maxLeftX - minLeftX)
    {
        if(maxLeftY-minLeftY > maxLeftZ - minLeftZ)
        {
            newLeftNode.longestAxis = 1;
            newLeftNode.midpoint = midLeft[1];
        }
    }
    if(maxLeftZ - minLeftZ > maxLeftX - minLeftX)
    {
        if(maxLeftZ - minLeftZ > maxLeftY-minLeftY)
        {
            newLeftNode.longestAxis = 2;
            newLeftNode.midpoint = midLeft[2];
        }
    }
    
    if(maxRightX-minRightX > maxRightY - minRightY)
    {
        if(maxRightX-minRightX > maxRightZ - minRightZ)
        {
            newRightNode.longestAxis = 0;
            newRightNode.midpoint = midRight[0];
        }
    }
    if(maxRightY-minRightY > maxRightX - minRightX)
    {
        if(maxRightY-minRightY > maxRightZ - minRightZ)
        {
            newRightNode.longestAxis = 1;
            newRightNode.midpoint = midRight[1];
        }
    }
    if(maxRightZ-minRightZ > maxRightX - minRightX)
    {
        if(maxRightZ-minRightZ > maxRightY - minRightY)
        {
            newRightNode.longestAxis = 2;
            newRightNode.midpoint = midRight[2];
        }
    }
    newLeftNode.minX = minLeftX;
    newLeftNode.maxX = maxLeftX;
    newLeftNode.minY = minLeftY;
    newLeftNode.maxY = maxLeftY;
    newLeftNode.minZ = minLeftZ;
    newLeftNode.maxZ = maxLeftZ;
    
    newRightNode.minX = minRightX;
    newRightNode.maxX = maxRightX;
    newRightNode.minY = minRightY;
    newRightNode.maxY = maxRightY;
    newRightNode.minZ = minRightZ;
    newRightNode.maxZ = maxRightZ;
    
    int l = constructTree(leftObjects, newLeftNode,nodes);
    int r = constructTree(rightObjects, newRightNode,nodes);
    
    currentNode.left = l;
    currentNode.right = r;
    
    nodes.push_back(currentNode);
    return nodes.size()-1;
}


#include "bvh.hpp"
