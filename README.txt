Subject 	: CSCI420 - Computer Graphics 
Assignment 3: Ray Tracing
Author		: Dhruv Tarsadiya

Core Credit Features: 
======================
All core credits and requirements satisfied. 
Still Images saved in ../hw3/stills
To run code: 
    cd ../../hw3
    make 
    ./hw3 <FILE_NAME>.scene <Optional Screenshot File>


Extra Credit Features: 
======================
1. Recursive Reflection: 
    Allowed recursive reflection with  I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)  and  (1 - ks/rec_ref_coef) * localPhongColor + ks/rec_ref_coef * colorOfReflectedRay
    where rec_ref_coef controls the degree/ extent of the recursive effect wanted. It and also MAX_RECURSIVE_DEPTH can be changed to vary recursive reflection effect. See implementation in
    computeColor()

2. Good Antialiasing:
    See draw.scene() for implementation. Adjust samples_per_pixel for super sampling vs efficiency. Samples drawn from uniform distrubtion for each sample. Effect of each sample normalised
    per pixel. Introduces smoothness and reduces jaggedness. 

3. Soft Shadows:
    See computeColor() and Light struct for implementation. Modified Light Struct to have default area and num_samples per light source. 1x1 light sourcegrid divided into further minigrids for 
    drawing samples with offsets. These samples are mini light sources whose effect is normalised per area light source (i.e. per light source in scene file ) Used the center of the area light 
    source for intersection calculation. The proportion of unblocked shadow rays per light source determines the softness of the shadow. 

4. Included vec3.h helper class with overloaded operators to help with defining colors, points and directions. 
5. Included extra scene files like test3,test4 etc. See 007.jpg onwards 

Keyboard/Mouse controls: (Please document Keyboard/Mouse controls if any)
1. Exit screen with ESC

Comments : Used https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm as reference for triangle ray intersection i.e. hitTraingle() method in code.
           Change focal_length to change camera focal len. 
           Change background_color if needed 
           

