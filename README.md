raytra
======

Simple C++ raytracer


Requirements
-------------
* linux
* gcc
* *OpenEXR library and headers


How to compile
------------
Including is a makefile.  Just type "make" in the directory to compile.  Make sure you have the OpenEXR library and headers installed and change the library path in the Makefile if it's not the default (/usr/include/OpenEXR).

Scene files
------------
The scenefiles directory contains different scenes to render.  They contain all the objects, lighting and position information to pass to the renderer.  The demo_scene.txt shows all the features I implemented: glossy reflection, area lights/soft shadows, anti-aliasing and refraction.  The lexus.txt will demonstrate the speed of the BVH structure, since it contains about 40K triangles.


File descriptions
-----------------
raytra.cpp: main application control, just a few lines to call the parser, and then the renderer.
objects.cpp: Contains all the object classes which includes lights, surfaces, boundingboxes and the bvh implementation (as BvhNode).
readscene.cpp: reads the scene files.
scene.cpp: the main renderer class.
support.cpp: support code which includes points, vectors, rays and the Utility class which implements some useful functions.
camera.cpp: camera class.


Custom options
---------------------------------------
The following is an explanation of the custom options available to be used in a scene file.  The demo_scene.txt file uses all of these options as an example.
The first is the "o" option.  It requires four parameters.  If no "o" is specified, the built-in default values are used.

The format is:
o [max_depth] [max_samples] [max_shadow_samples] [max_glossy_samples]

max_depth: the depth of raycasting recursion in the main renderer.  DEFAULT: 3
max_samples: how many samples per pixel to take (anti-aliasing).  It should be a perfect square (i.e. a value of 16 will sample 4 in the x direction, 4 in the y direction). DEFAULT: 1
max_shadow_samples: How many shadow samples to take when there is an area light in the scene.  If there is no area lights in the scene, this option is ignored. DEFAULT: 1
max_glossy_samples: How many glossy reflection samples to take when there is an object with glossy material in the scene.  If there are no objects with glossy materials, this option is ignored. DEFAULT: 3


The second option is the "g" option.  This represents the glossiness factor of the specular reflection component of the material.  It requires one parameter:

g [glossy_factor]

glossy_factor: value (ideally between 0 and 1) that represents how "glossy" (blurry) the specular reflection component will be.


The third option is the "r" option.  This represents the refractive property of the object (can really only be implemented with the sphere right now).  It requires 4 parameteres:

r [refraction_index] [kr] [kg] [kb]

refraction_index: represents how large the angle of refraction is as the ray passes through the object. 
kr, kg, kb: The red, green and blue values for the refractive constant (we assume the constants have the ln operation folded into them).
 

