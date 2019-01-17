//
//  main.cpp
//  RayTracer
//
//  Created by Bo Miller on 1/2/19.
//  Copyright © 2019 Bo Miller. All rights reserved.
//

#include <stdio.h>
#include <vector>

#include "glm/glm/glm.hpp"
#include "SceneObjects.hpp"
#include "Ray.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(int argc, const char * argv[]) {

    
    //*******Manually describing scene************
    float width = 720;
    float height = 360;
    std::vector<SceneObject> scene;
    
    SceneObject s1;
    s1.radius = 1.0;
    s1.ambient = glm::vec3(0,0,1);
    s1.specular = glm::vec3(0.9,0.4,0);
    s1.diffuse = glm::vec3(0.8,0.3,0.1);
    s1.center = glm::vec3(-4,3,0);
    s1.shininess = 64;
    s1.reflective = glm::vec3(.5,.5,.5);
    s1.sphere = true;
    scene.push_back(s1);
    
    SceneObject s2;
    s2.radius = 1.5;
    s2.ambient = glm::vec3(1.0,1,1);
    s2.specular = glm::vec3(.5,.5,.5);
    s2.diffuse = glm::vec3(1,1,1);
    s2.center = glm::vec3(-2,0,0);
    s2.reflective = glm::vec3(.5,.5,.5);
    s2.shininess = 64;
    s2.sphere = true;
    scene.push_back(s2);
    
    SceneObject s3;
    s3.radius = 1;
    s3.ambient = glm::vec3(1.0,1.0,0.0);
    s3.specular = glm::vec3(.5,.5,.5);
    s3.diffuse = glm::vec3(1,1,1);
    s3.center = glm::vec3(1,0,3);
    s3.reflective = glm::vec3(.5,.5,.5);
    s3.shininess = 64;
    s3.sphere = true;
    scene.push_back(s3);
    
    SceneObject s4;
    s4.radius = 1;
    s4.ambient = glm::vec3(1.0,1.0,0.0);
    s4.specular = glm::vec3(.5,.5,.5);
    s4.diffuse = glm::vec3(1,1,1);
    s4.center = glm::vec3(2,0,3);
    s4.reflective = glm::vec3(.5,.5,.5);
    s4.shininess = 64;
    s4.sphere = true;
    scene.push_back(s4);
    
    SceneObject s5;
    s5.radius = .75;
    s5.ambient = glm::vec3(1.0,1.0,0.0);
    s5.specular = glm::vec3(.5,.5,.5);
    s5.diffuse = glm::vec3(1,1,1);
    s5.center = glm::vec3(1.5,1,3);
    s5.reflective = glm::vec3(.5,.5,.5);
    s5.shininess = 64;
    s5.sphere = true;
    scene.push_back(s5);
    
    SceneObject s6;
    s6.radius = .75;
    s6.ambient = glm::vec3(1.0,1.0,0.0);
    s6.specular = glm::vec3(.5,.5,.5);
    s6.diffuse = glm::vec3(1,1,1);
    s6.center = glm::vec3(1.5,-1,3);
    s6.reflective = glm::vec3(.5,.5,.5);
    s6.shininess = 64;
    s6.sphere = true;
    scene.push_back(s6);
    
    SceneObject s7;
    s7.radius = 1.0;
    s7.ambient = glm::vec3(1.0,1,1);
    s7.specular = glm::vec3(.5,.5,.5);
    s7.diffuse = glm::vec3(1,1,1);
    s7.center = glm::vec3(0,2,0);
    s7.reflective = glm::vec3(.5,.5,.5);
    s7.shininess = 64;
    s7.sphere = true;
    scene.push_back(s7);
    
    SceneObject s8;
    s8.radius = 1;
    s8.ambient = glm::vec3(0,1,0);
    s8.specular = glm::vec3(.2,.2,.2);
    s8.diffuse = glm::vec3(1,1,1);
    s8.center = glm::vec3(3,3,4);
    s8.reflective = glm::vec3(.5,.5,.5);
    s8.shininess = 64;
    s8.sphere = true;
    scene.push_back(s8);
    
    SceneObject s9;
    s9.radius = .75;
    s9.ambient = glm::vec3(.9,.9,.9);
    s9.specular = glm::vec3(.2,.2,.2);
    s9.diffuse = glm::vec3(1,1,1);
    s9.center = glm::vec3(3,2,5.5);
    s9.reflective = glm::vec3(.5,.5,.5);
    s9.shininess = 64;
    s9.sphere = true;
    scene.push_back(s9);
    
    SceneObject s10;
    s10.radius = 1;
    s10.ambient = glm::vec3(0,0,1);
    s10.specular = glm::vec3(1,1,1);
    s10.diffuse = glm::vec3(1,1,1);
    s10.center = glm::vec3(-3,-.5,5);
    s10.reflective = glm::vec3(.5,.5,.5);
    s10.shininess = 64;
    s10.sphere = true;
    scene.push_back(s10);
    
    SceneObject s11;
    s11.radius = 1;
    s11.ambient = glm::vec3(1,0,1);
    s11.specular = glm::vec3(1,1,1);
    s11.diffuse = glm::vec3(1,1,1);
    s11.center = glm::vec3(-4.5,2,4);
    s11.reflective = glm::vec3(.5,.5,.5);
    s11.shininess = 64;
    s11.sphere = true;
    scene.push_back(s11);
    
    SceneObject s12;
    s12.radius = 1;
    s12.ambient = glm::vec3(0,1,1);
    s12.specular = glm::vec3(.5,.5,.5);
    s12.diffuse = glm::vec3(1,1,1);
    s12.center = glm::vec3(2.5,3.5,-3.5);
    s12.reflective = glm::vec3(.3,.3,.3);
    s12.shininess = 64;
    s12.sphere = true;
    scene.push_back(s12);
    
    glm::vec3 cameraPosition = glm::vec3(0,4,12);
    glm::vec3 cameraDirection = glm::vec3(0,0,0);
    //************************************
    
    
    //**************Adding Lights*****************
    std::vector<Light> lights;
    Light l1;
    l1.direction = glm::vec3(0, -1, 0);
    l1.color = glm::vec3(1.0, 1.0, 1.0);
    l1.point = false;
    lights.push_back(l1);
    
    Light l2;
    l2.point = true;
    l2.position = glm::vec3(4, 3, -4);
    l2.color = glm::vec3(.5, .5, .5);
    l2.constantTerm = 0.25f;
    l2.linearTerm = 0.003372407f;
    l2.quadraticTerm = 0.000045492f;
    lights.push_back(l2);
    
    Light l3;
    l3.direction = glm::vec3(0,1,0);
    l3.color = glm::vec3(0.2,0.2,0.2);
    l3.point = false;
    lights.push_back(l3);
    //********************************************
    
    
    /*Declare pixel buffer and start ray tracing*/
    float buffer[360][720][3];
    startRayTracing(width, height, buffer, cameraPosition, cameraDirection, scene, lights);
    
    
    
    /*Writes pixel buffer to a .bmp file*/
    {
        int w = (int) width;
        int h = (int) height;
        
        uint8_t pix[w*h*3];
        memset(pix, 0xff, w*h*3);
        for(int i = 0; i<h; i++)
        {
            for(int j = 0;j<w;j++)
            {
                int elem = (i*w*3) + (j*3);
                pix[elem+0] = (uint8_t) (buffer[i][j][0] * 255.0);
                pix[elem+1] = (uint8_t) (buffer[i][j][1] * 255.0);
                pix[elem+2] = (uint8_t) (buffer[i][j][2] * 255.0);

            }
        }
        stbi_write_bmp("/tmp/fb.bmp", w, h, 3, pix);
    }
    
    return 0;
}