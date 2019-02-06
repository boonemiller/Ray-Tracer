//
//  bvh.hpp
//  RayTracer
//
//  Created by Bo Miller on 2/1/19.
//  Copyright © 2019 Bo Miller. All rights reserved.
//
#include "SceneObjects.hpp"
#include <vector>
#ifndef bvh_hpp
#define bvh_hpp

class Node
{
public:
    int left;
    int right;
    bool isleaf = false;
    int objs[3];
    int numObjs;
    //some bounding box variables
    double minX;
    double maxX;
    double minY;
    double maxY;
    double minZ;
    double maxZ;
    
    double midpoint;
    double longestAxis;
};

int constructTree(std::vector<SceneObject>& objects, Node& currentNode, std::vector<Node>& nodes);


#include <stdio.h>

#endif /* bvh_hpp */
