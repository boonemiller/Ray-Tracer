# Ray-Tracer
Ray Tracer implementation in C++, Random Sample AA, multi-threading, bvh acceleration, and runtime comparisons on different CPUs

### Introduction

Described in this README and the code in this repo is a personal project on a basic raytracing implementation.
Added features on top of the basic ray tracer include random sample per pixel anti-aliasing, multithreading support, and a bvh acceleration structure. 
I described the features of the ray-tracer and included a performance report of multithreading the ray tracer on three different Intel processors.

![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/16AA.bmp)

### Phong-shading

I implemented standard phong shading in the ray tracer, look at reference 3. Used the following formula for the shading

Ambient * Ka + Distance attenuation * Shadow attenuation * sum across lights(diffuse*Kd + Specular*Ks) + Kr*Results of Secondary rays	 

Where Ka, Kd, Ks, and Kr are all constants that change for each object.

### Relevant Terms

Ambient- the inherent color of the object

Diffuse- a color that changes intensity with the direction of the surface.

Specular- highlights, only show up if there is a particular angle between a reflected light ray and the camera

Distance Attenuation- how the distance from the light source affects the intensity of the light

Shadow Attenuation- if a point in the scene is in shadow or not

Point lights- a particular point in space that is emitting light. Has a color and terms that affect the Distance Attenuation

Directional Lights - a light that comes from an infinite location and only has a direction and color.

Bounce Depth- number of ray bounces we perform.

### Ray-tracer Steps

#### Primary Rays
The ray tracer starts by casting rays from the location of the camera through each pixel. This initial ray is called the primary ray. 
In my implementation, the ray tracer first tests a primary ray through each pixel to see if it intersects any of the spheres in the scene. 
If it does, this is called an intersection point to be used later when casting secondary rays and shadow rays. 
If it doesn’t intersect any of the spheres, we know that it intersects one of the walls, defined as oriented planes within the space. 
The ray-plane intersection tests finds which plane it intersects and finds an intersection point for that plane.   

At the intersection point, either a sphere or a plane intersection point, we can get the initial ambient color of the object. From the intersection point we can cast secondary and shadow rays. 

#### Shadow Rays
At each intersection point, we need to see if that point is in shadow or not. To do this we cast another ray from the intersection point in the direction of each light. If we intersect any of the spheres we know that it is in shadow and the light doesn’t contribute anything to the color of the intersection point. 

If it isn’t in shadow for any particular light than the light contributes its color, a diffuse, and specular term. For a point light we have to determine how far away this light is and take it into account.

For in-shadow intersection points, shadow attenuation is 0, for out of shadow intersection points shadow attenuation is 1. 

For out of shadow intersection points we can also calculate the distance attenuation.
Directional lights- our distance attenuation is always 1. 
Point Lights- our distance attenuation is 1/(ConstantTerm + distance*LinearTerm + (distance^2) * QuadraticTerm)

#### Secondary Rays
Now from the initial intersection point, we can cast secondary rays from the intersection point, based on the incoming vector and the intersection normal. We treat these rays the same way we treated the initial casted primary rays, testing the sphere intersections first, then the wall intersections. We repeat this process bounce depth of times. 


After we have reached our bounce depth we return the value calculated from the phong shading equation listed above.

### Random Sample Anti-Aliasing

To add anti-aliasing into the process described above. Instead of only casting one primary ray through the center of the pixel. We casts some number of primary rays through each pixel at random locations and average the results returned by each of these rays. The results can be seen around the edges of the objects and at the intersection of the walls, and causes them to be smoother and not appear as jagged. This is because if an edge of an object passes through a pixel and we sample the part of pixel that it isn’t in, we don’t get the contribution of the object in the pixel. Random sample anti-aliasing increases the chance of an intersection and including the result of the edge in the pixel.


![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/16AA.bmp)
![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/AAOFF.bmp)

Scene with 16 Random samples per pixel (top image), AA-Off one sample per pixel (bottom image) 

### Noise

Noise is a common problem in ray tracers which is why there are plenty of methods to reduce noise in rendered scenes. The only method I implemented was anti-aliasing as described above. As you can see in the pictures, rendering with anti-aliasing(top image) reduced the noise present in the one sample per pixel render (bottom image).

### Multi-Threading 

Adding multithreading support is fairly simple. Ray tracers are known to be embarrassingly parallel. I used p-threads and divided up the frame into chunks of work based on the width. For example our 720x360 frame across two threads gets divided up into two 360X360 vertical chunks. And I extended this to 4, 6, and 8 threads if the processor could handle it for the performance profiling across different processors.


![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/profiledScene.bmp)

### Processor Comparison on Multithreaded Ray Tracing

I ran the above scenes raytracing code on three different Intel processors. A intel dual-core 3.1 GHz i5, a 2.5 GHz quad-core i7, and a 2.9 Ghz quad-core i7. Rather unsurprisingly the results are what should be expected. The 3.1 GHz i5 beat out the 2.5 GHz i7 on single and double threaded runs, but starts to lose starting with 3 threads. And the 2.9 GHz i7 was faster than both of the other processors on every run.

![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/cpu-performance-graph.png)


### BVH Acceleration Structure

I implemented a BVH tree acceleration structure to decrease the number of object intersection tests needed for each ray casted. The basic process for this is that a ray traveling through a scene will only get close enough to a few objects and mostly miss all of the other objects. So in order to decrease the amount of these intersection tests for the missed objects. I create bounding boxes around some number of the objects. If the ray intersects a box then it is possible that it will intersect the objects inside the box. If it doesn’t intersect the bounding box then it can't intersect any of the objects within the box. I make a tree out of these bounding boxes where the root represents the whole scene’s bounding box, and each node represents some smaller bounding box with some subset of the objects. I stop building the tree when I get to a small number of objects within a bounding box, I said 3 or less, these represent the leaves of the tree. 

So now when a ray is traveling through a scene, I test whether or not the ray intersects any of the bounding boxes that represent the objects. If it intersects any of the leaf bounding boxes I can tests only those objects (3 or less), instead of the entire scene.

![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/nvidia-bvh.png)

This is an example of how the tree might look if you used one object per leaf (taken from nvidia).


![alt text](https://raw.githubusercontent.com/boonemiller/Ray-Tracer/master/RayTracer/bvh-performance.png)

Graph of the speed ups achieved by adding a bvh acceleration structure into my raytracing code. This shows the percentage speed up of using a bvh acceleration structure over no acceleration structure at all, on a variety of scenes. 
