# CS580-RayTracer

Example input 

1. ./ray_tracer ./cornell-box/cornell-box/cornellScene.txt ./cornell-box/cornell-box.txt
use cornellScene setting and mesh data to render cornell box just display the result 

2. ./ray_tracer NAME_OF_OBJFILE.obj 
Use application as parser of obj file and will generate a Mesh.txt file that is in format of our input 

Parameters:
SS_ON: control soft shadowing true = on; false = off 
reflect_num: reflectivity 
SAMPLE_RATE: antialiasing
material: 0 flat color, 1 marble (passed in from mesh file)
---------------------------
scenefile format example: 

ambient 0.3 0.3 0.3
eye 0 0 0                               # camera setting would be were ray is orginated from 
look 0 0 -1
up 0 1 0
light
1                                        # number of light source 
pos 0 0 0
col 1 1 1 

---------------------------
meshfile format example: 

3                               # number of object/ Mesh
Mesh                                # a new object 
Material: 0                             # material 0 as flat color, 1 as marble, need to be set manually if wanna do procedural texturing
2                                           # number of triangle inside this triangle mesh
triangle
pos: -0.884011 2.3185 -6.56797 
nor: 0 -1 0.0008 
dif: 1 1 1 
spe: 0 0 0 
shi: 0
triangle
pos: -0.884011 2.3185 -6.56797 
nor: 0 -1 0.0008 
dif: 1 1 1 
spe: 0 0 0 
shi: 0
sphere                               # a new object
Material: 0                             # material 0 as flat color, 1 as marble, need to be set manually if wanna do procedural texturing
pos: 15.0 15.0 -2.0
rad: 10
dif: 0.3 0.3 0.3
spe: 0.0 0.0 0.0
shi: 1
sphere                               # a new object
Material: 0                             # material 0 as flat color, 1 as marble, need to be set manually if wanna do procedural texturing
pos: 10.0 10.0 -2.0
rad: 10
dif: 0.3 0.3 0.3
spe: 0.0 0.0 0.0
shi: 1