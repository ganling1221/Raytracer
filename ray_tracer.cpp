/*
CSCI 580 Final Project
RayTracer
*/

#include <string.h>
#include <random>
#include "OBJ_Loader.h"
#include "Vector.h"
#include "Marble.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <algorithm>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#define MAX_TRIANGLES 8000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define RANGE 0.0001

unsigned char buffer[HEIGHT][WIDTH][3];

struct Ray
{
	Vector origin;
	Vector direction;

	Ray() : origin(), direction()
	{}

	Ray(Vector ori, Vector dir) : origin(ori), direction(dir)
	{}
};
struct Texture 
{
	double x;
	double y;	
};
struct Vertex
{
	Vector position;
	Texture texture;
	Vector color_diffuse;
	Vector color_specular;
	Vector normal;
	double shininess;
};

typedef struct _Triangle
{
	struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
	Vector position;
	Vector color_diffuse;
	Vector color_specular;
	double shininess;
	double radius;
	int material;

} Sphere;

typedef struct _Light
{
	Vector position;
	Vector color;
} Light;

typedef struct _Object 
{
	Triangle triangles[MAX_TRIANGLES];
	int material;
	//material =0 means flat color
	//material =1 means prodcedural mapped with marble 
}Object;

Object objects[MAX_SPHERES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
Vector ambient_light;
Vector eye;
Vector look;
Vector up;

Vector** pixels;

int num_objects=0;
int num_triangles=0;
int num_spheres=0;
int num_lights=0;

/* for anti-alising */
const int SAMPLE_RATE = 4;

/* for recursive ray */
const int reflect_num = 3;
const double ratio = 0.1; 	// determine ratio of reflected color
const double ERR = 1e-12; 	// add to reflection origin to avoid intersection with reflector

/* for soft shadow */
const bool SS_ON = true; 
const double AREA_LIGHT_WIDTH = 2.0;
const double AREA_LIGHT_LENGTH = 2.0;

const double LIGHT_CELL_SIZE = 4.0;

/* view parameters */
double aspect_ratio = (double) WIDTH / (double) HEIGHT;
double view_min_y = -tan(fov * M_PI / 2.0 / 180);
double view_min_x = aspect_ratio * view_min_y;
double view_z = -1.0;
double view_width = -2 * view_min_x;
double view_height = -2 * view_min_y;

void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

void ray_tracer(int x, int y, Vector &res_color);

Vector get_color(Ray ray, int cnt);

bool get_object_intersection(Ray ray, int &obj_type, int &obj_id, int &tri_id, double &dist);
bool get_light_intersection(Ray ray, int li_id, double &dist);
bool get_sphere_intersection(Ray ray, int sp_id, double &dist);
bool get_triangle_intersection(Ray ray, int object_id,  int tr_id, double &dist);

// for reflective ray
Vector get_sphere_color(Ray ray, int sp_id, double inter_dist, Ray &ref_ray, Vector &ref_ks);

// for reflective ray
Vector get_triangle_color(Ray ray, int object_id, int tr_id, double inter_dist, Ray &ref_ray, Vector &ref_ks);

Vector color_PhongIlluminate(Vector point, Vector normal, Vector viewer, Vector kd, Vector ks, double shi);

/* ---------------------- ASSISTANT FUNCTIONS ---------------------- */

double get_dist(Vector v1, Vector v2)
{
	double sum = (v1.x-v2.x)*(v1.x-v2.x) + (v1.y-v2.y)*(v1.y-v2.y) + (v1.z-v2.z)*(v1.z-v2.z);
	return sqrt(sum);
}

double dot_product(Vector v1, Vector v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

Vector cross_product(Vector v1, Vector v2)
{
	Vector res;

	res.x = v1.y * v2.z - v1.z * v2.y;
	res.y = v1.z * v2.x - v1.x * v2.z;
	res.z = v1.x * v2.y - v1.y * v2.x;

	return res;
}

/* ----------------------------------------------------------------- */

void draw_scene()
{
	// init pixels to store colors
	pixels = new Vector* [WIDTH * SAMPLE_RATE];
	for (int i = 0; i < WIDTH * SAMPLE_RATE; i++)
	{
		pixels[i] = new Vector [HEIGHT * SAMPLE_RATE];
	}

	// get colors
	for (int i = 0; i < WIDTH * SAMPLE_RATE; i++) 
	{
		for (int j = 0; j < HEIGHT * SAMPLE_RATE; j++)
		{
			ray_tracer(i, j, pixels[i][j]);
		}
	}

	// print pixels
	for(int x = 0; x < WIDTH; x++)
	{
		glPointSize(2.0);  
		glBegin(GL_POINTS);
		for(int y = 0; y < HEIGHT; y++)
		{
			Vector color;
			// average pixel color for anti-alising
			if (SAMPLE_RATE > 1) 
			{
				for (int i = 0; i < SAMPLE_RATE; i++)
				{
					for (int j = 0; j < SAMPLE_RATE; j++)
					{
						color += pixels[x*SAMPLE_RATE + i][y*SAMPLE_RATE + j];
					}
				}
				color /= (SAMPLE_RATE*SAMPLE_RATE);
			}
			else
			{
				color = pixels[x][y];
			}

			plot_pixel(x, y, color.x, color.y, color.z);
		}
		glEnd();
		glFlush();
	}
	printf("Done!\n"); fflush(stdout);
}

void ray_tracer(int x, int y, Vector &res_color)
{
	double view_x = view_min_x + (x*1.0 / (WIDTH*SAMPLE_RATE)) * view_width;
	double view_y = view_min_y + (y*1.0 / (HEIGHT*SAMPLE_RATE)) * view_height;
	Vector wDir = Vector(-look.x,-look.y,-look.z);
	 wDir.normalize(); 
     Vector uDir = cross_product(up, wDir);
	 uDir.normalize(); 
     Vector vDir = cross_product(wDir, uDir);  

	Vector ray_dir=(uDir*view_x + vDir*view_y - wDir);
	ray_dir.normalize(); // NORMALIZE!

	Ray cur_ray(eye, ray_dir);

	// for recursion
	res_color = get_color(cur_ray, reflect_num);

	res_color += ambient_light * 255.0;

	res_color.x = std::min((std::max(res_color.x, 0.0)), 255.0); 
	res_color.y = std::min((std::max(res_color.y, 0.0)), 255.0); 
	res_color.z = std::min((std::max(res_color.z, 0.0)), 255.0); 

	return;
}

Vector get_color(Ray ray, int cnt)
{
	int obj_type = -1;
	int obj_id = -1;
	int tri_id = -1;

	double inter_dist = DBL_MAX;

	bool is_interseced = get_object_intersection(ray, obj_type, obj_id,tri_id, inter_dist);

	Vector res_color;

	Ray reflect_ray;
	Vector reflect_ks;

	if (!is_interseced)
	{
		// set background color to white
		res_color = Vector(255.0, 255.0, 255.0); 
		return res_color;
	}

	// light
	if (obj_type == 0)
	{
		res_color = lights[obj_id].color * 255.0;
		// if encounter light, no reflected rays
		return res_color;
	}
	// sphere
	else if (obj_type == 1)
	{
		// printf("get sphere color\n");
		res_color = get_sphere_color(ray, obj_id, inter_dist, reflect_ray, reflect_ks);
	}
	// triangle
	else if (obj_type == 2)
	{
		res_color = get_triangle_color(ray, obj_id,tri_id, inter_dist, reflect_ray, reflect_ks);
	}
	// this shouldn't happen
	else
	{
		printf("There's something wrong with get_object_intersecction!\n");
		exit(1);
	}

	// add recursion to compute color of reflected ray
	// only when cnt > 1, there are reflected rays
	if (cnt > 0) {
		Vector reflect_color = get_color(reflect_ray, cnt-1);

		Vector unit_vector =  Vector(1.0, 1.0, 1.0);
		res_color = res_color * (1.0 - ratio) + reflect_color * reflect_ks * ratio;
	}

	return res_color;
}

// done
bool get_object_intersection(Ray ray, int &obj_type, int &obj_id, int &tri_id, double &dist)
{
	bool is_intersected = false;
	double cur_dist = DBL_MAX;

	for (int i = 0; i < num_lights; i++)
	{
		if (get_light_intersection(ray, i, cur_dist))
		{
			if (cur_dist < dist)
			{
				is_intersected = true;
				dist = cur_dist;	
				obj_type = 0;
				obj_id = i;
			}
		}
	}

	for (int i = 0; i < num_spheres; i++)
	{
		if (get_sphere_intersection(ray, i, cur_dist))
		{
			if (cur_dist < dist)
			{
				is_intersected = true;
				dist = cur_dist;
				obj_type = 1;
				obj_id = i;
			}
		}
	}
	for(int j =0; j< num_objects;j++){
		for (int i = 0; i < num_triangles; i++)
		{
			if (get_triangle_intersection(ray, j, i, cur_dist))
			{
				if (cur_dist < dist)
				{
					is_intersected = true;
					dist = cur_dist;
					obj_type = 2;
					obj_id = j;
					tri_id = i;

				}
			}
		}
	}
	

	return is_intersected;
}

// done
bool get_light_intersection(Ray ray, int li_id, double &dist) 
{
	if (ray.origin == lights[li_id].position)
	{
		return false;
	}

	// compare light direction to ray direction
	Vector light_dir = lights[li_id].position - ray.origin;
	Vector dir_cmp = light_dir / ray.direction;

	if (abs(dir_cmp.x - dir_cmp.y) < RANGE)
	{
		if (abs(dir_cmp.y - dir_cmp.z) < RANGE)
		{
			if (dir_cmp.x > 0)
			{
				dist = dir_cmp.x;
				return true;
			}
		}
	}

	return false;
}

// done
bool get_sphere_intersection(Ray ray, int sp_id, double &dist)
{
	Vector tmp_b = ray.direction * (ray.origin - spheres[sp_id].position);
	double b = 2.0 * (tmp_b.x + tmp_b.y + tmp_b.z);

	Vector tmp_c = ray.origin - spheres[sp_id].position;
	double c = pow(tmp_c.x, 2) + pow(tmp_c.y, 2) + pow(tmp_c.z, 2) - pow(spheres[sp_id].radius, 2);

	double delta = pow(b, 2) - 4.0 * c;
	if (delta < 0.0)
	{
		return false;
	}

	double t0 = (-b + sqrt(delta)) / 2.0;
	double t1 = (-b - sqrt(delta)) / 2.0;

	if (t0 <= 0 && t1 <= 0) {
		return false;
	}
	else if (t0 > 0 && t1 >0)
	{
		dist = t0 < t1? t0 : t1; 
	}
	else
	{
		dist = t0 > t1? t0 : t1;
	}

	if (dist < RANGE)
	{
		return false;
	}

	return true;	
}

// done
bool get_triangle_intersection(Ray ray, int object_id, int tr_id, double &dist)
{
	Vector p1 = objects[object_id].triangles[tr_id].v[0].position;
	Vector p2 = objects[object_id].triangles[tr_id].v[1].position;
	Vector p3 = objects[object_id].triangles[tr_id].v[2].position;

	Vector normal = cross_product(p2 - p1, p3 - p1);
	normal.normalize();

	double n_dot_d = dot_product(normal, ray.direction);

	if (n_dot_d == 0.0)
	{
		return false;
	}

	double d = -1.0 * dot_product(normal, p1);
	dist = - (dot_product(normal, ray.origin) + d) / n_dot_d;

	if (dist < RANGE)
	{
		return false;
	}

	Vector p = ray.origin + ray.direction * dist;
	double area_p1p2p3 = dot_product(cross_product(p2-p1, p3-p1), normal);
	double area_pp2p3 = dot_product(cross_product(p2-p, p3-p), normal);
	double area_p1pp3 = dot_product(cross_product(p-p1, p3-p1), normal);
	double area_p1p2p = dot_product(cross_product(p2-p1, p-p1), normal);

	if (area_pp2p3 >= 0 && area_p1pp3 >= 0 && area_p1p2p >= 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

// done
Vector get_sphere_color(Ray ray, int sp_id, double inter_dist, Ray &ref_ray, Vector &ref_ks)
{
	Marble marble =  Marble();
	Vector inter_point = ray.origin + ray.direction * inter_dist;

	Vector normal = inter_point - spheres[sp_id].position;
	normal.normalize();

	Vector viewer = ray.direction * -1.0;
	// viewer.normalize();
	
	Vector kd = spheres[sp_id].color_diffuse;
	Vector ks = spheres[sp_id].color_specular;
	double shi = spheres[sp_id].shininess;
	if(spheres[sp_id].material == 1){
		 kd = marble.getColor(inter_point);
		 shi = marble.getShininess();
	}
	Vector phong_color = color_PhongIlluminate(inter_point, normal, viewer, kd, ks, shi);

	// compute reflected ray
	double l_dot_n = dot_product(viewer, normal);
	l_dot_n = l_dot_n < 0? 0 : l_dot_n;

	ref_ray.direction = normal * (2 * l_dot_n) - viewer;
	ref_ray.direction.normalize();

	ref_ray.origin = inter_point + normal * ERR;

	ref_ks = ks;

	return phong_color * 255.0;
}

// done
Vector get_triangle_color(Ray ray, int object_id, int tr_id, double inter_dist, Ray &ref_ray, Vector &ref_ks)
{
	Marble marble =  Marble();
	Vertex v1 = objects[object_id].triangles[tr_id].v[0];
	Vertex v2 = objects[object_id].triangles[tr_id].v[1];
	Vertex v3 = objects[object_id].triangles[tr_id].v[2];

	Vector p1 = v1.position;
	Vector p2 = v2.position;
	Vector p3 = v3.position;

	Vector tri_normal = cross_product(p2 - p1, p3 - p2);
	tri_normal.normalize();

	Vector p = ray.origin + ray.direction * inter_dist;
	double area_p1p2p3 = dot_product(cross_product(p2-p1, p3-p1), tri_normal);
	double area_pp2p3 = dot_product(cross_product(p2-p, p3-p), tri_normal);
	double area_p1pp3 = dot_product(cross_product(p-p1, p3-p1), tri_normal);
	double area_p1p2p = dot_product(cross_product(p2-p1, p-p1), tri_normal);

	double alpha = area_pp2p3 / area_p1p2p3;
	double beta = area_p1pp3 / area_p1p2p3;
	double gamma = 1.0 - alpha - beta;

	Vector inter_point = ray.origin + ray.direction * inter_dist;

	Vector normal = v1.normal * alpha + v2.normal * beta + v3.normal * gamma;
	normal.normalize();

	Vector viewer = ray.direction * -1.0;
	// viewer.normalize();

	Vector kd = v1.color_diffuse * alpha + v2.color_diffuse * beta + v3.color_diffuse * gamma;

	Vector ks = v1.color_specular * alpha + v2.color_specular * beta + v3.color_specular * gamma;
	double shi = v1.shininess * alpha + v2.shininess * beta + v3.shininess * gamma;
	
	if(objects[object_id].material==1){
		kd = marble.getColor(inter_point); 
		shi = marble.getShininess();
	}

	Vector phong_color = color_PhongIlluminate(inter_point, normal, viewer, kd, ks, shi);

	// compute reflected ray
	double l_dot_n = dot_product(viewer, normal);
	l_dot_n = l_dot_n < 0? 0 : l_dot_n;

	ref_ray.direction = normal * (2 * l_dot_n) - viewer;
	ref_ray.direction.normalize();

	ref_ray.origin = inter_point + normal * ERR;

	ref_ks = ks;

	return phong_color * 255.0;
}

// done
Vector color_PhongIlluminate(Vector point, Vector normal, Vector viewer, Vector kd, Vector ks, double shi)
{
	Vector color;
	for (int i = 0; i < num_lights; i++)
	{
		double intensity = 0;
		if(SS_ON){
			double cellMinX = lights[i].position.x - AREA_LIGHT_WIDTH /2;
			double cellMaxX = lights[i].position.x + AREA_LIGHT_WIDTH /2;

			double cellMinZ = lights[i].position.z - AREA_LIGHT_LENGTH /2;
			double cellMaxZ = lights[i].position.z + AREA_LIGHT_LENGTH /2;

			double cellXStep = (cellMaxX - cellMinX) / LIGHT_CELL_SIZE;
			double cellZStep = (cellMaxZ - cellMinZ) / LIGHT_CELL_SIZE;

			double sample_count = 0;		

			for(double cellX = cellMinX; cellX<=cellMaxX-cellXStep; cellX+= cellXStep){
				for(double cellZ = cellMinZ; cellZ<=cellMaxZ-cellZStep; cellZ+= cellZStep){
					double cellXRand = cellX + (double)(rand()) / ((double)(RAND_MAX/(cellXStep)));
					double cellZRand = cellZ + (double)(rand()) / ((double)(RAND_MAX/(cellZStep)));

					Vector ray_dir(cellXRand,lights[i].position.y,cellZRand);
					// Vector ray_dir(cellX,lights[i].position.y,cellZ);
					ray_dir = ray_dir - point;
					ray_dir.normalize();

					Ray shadow_ray(point, ray_dir);
					double dist_to_light = get_dist(lights[i].position, point);

					bool is_blocked = false;
					double obj_dist = DBL_MAX;

					
					for (int j = 0; j < num_spheres; j++)
					{
						if (get_sphere_intersection(shadow_ray, j, obj_dist))
						{
							if (obj_dist < dist_to_light)
							{
								is_blocked = true;
								break;
							}
						}
					}
					for(int i =0; i < num_objects;i++)
					{
						for (int j = 0; j < num_triangles; j++)
						{
							if (get_triangle_intersection(shadow_ray, i, j, obj_dist))
							{
								if (obj_dist < dist_to_light)
								{
									is_blocked = true;
									break;
								}
							}
						}
					}				
				

					if (!is_blocked) intensity++;

					sample_count++;
				}
			}
			intensity /= sample_count;
		}
		else{
			Vector ray_dir = lights[i].position - point;
			ray_dir.normalize();

			Ray shadow_ray(point, ray_dir);
			double dist_to_light = get_dist(lights[i].position, point);

			bool is_blocked = false;
			double obj_dist = DBL_MAX;

			
			for (int j = 0; j < num_spheres; j++)
			{
				if (get_sphere_intersection(shadow_ray, j, obj_dist))
				{
					if (obj_dist < dist_to_light)
					{
						is_blocked = true;
						break;
					}
				}
			}
			for(int i =0; i < num_objects;i++)
			{
				for (int j = 0; j < num_triangles; j++)
				{
					if (get_triangle_intersection(shadow_ray, i, j, obj_dist))
					{
						if (obj_dist < dist_to_light)
						{
							is_blocked = true;
							break;
						}
					}
				}
			}				

			if (!is_blocked) intensity = 1.0;
		}
		
		Vector ray_dir = lights[i].position - point;
		ray_dir.normalize();
		
		double l_dot_n = dot_product(ray_dir, normal);
		l_dot_n = l_dot_n < 0? 0 : l_dot_n;

		Vector r = normal * (2 * (dot_product(ray_dir, normal))) - ray_dir;
		r.normalize();

		double r_dot_v = dot_product(r, viewer);
		r_dot_v = r_dot_v < 0? 0 : r_dot_v;

		color += (lights[i].color * (((kd * l_dot_n) + ks * pow(r_dot_v, shi) ))) * intensity;
	}

	return color;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
	glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
	glVertex2i(x,y);

}

void parse_check(char *expected,char *found)
{
	if(strcasecmp(expected,found))
	{
		char error[100];
		printf("Expected '%s ' found '%s '\n",expected,found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}
}

void parse_doubles(FILE*file, char *check, double p[3])
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
	printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

// overload to accept Vector
void parse_doubles(FILE*file, char *check, Vector &p)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf %lf", &p.x, &p.y, &p.z);
	printf("%s %lf %lf %lf\n", check, p.x, p.y, p.z);
}
// overload to accept Vector
void parse_doubles2D(FILE*file, char *check, Texture &p)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf", &p.x, &p.y);
	printf("%s %lf %lf\n", check, p.x, p.y);
}

void parse_rad(FILE*file,double *r)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check("rad:",str);
	fscanf(file,"%lf",r);
	printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
	char s[100];
	fscanf(file,"%s",s);
	parse_check("shi:",s);
	fscanf(file,"%lf",shi);
	printf("shi: %f\n",*shi);
}

void parse_mat(FILE*file,int *mat)
{
	char s[100];
	fscanf(file,"%s",s);
	parse_check("Material:",s);
	fscanf(file,"%d",mat);
	printf("Material: %d\n",*mat);
}

int loadScene(char *argv)
{
	FILE *file = fopen(argv,"r");
	int number_of_Mesh;

	char type[50];
	int m;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file,"%i",&number_of_Mesh);

	printf("number_of_Mesh: %i\n",number_of_Mesh);
	

	char str[200];

	for(m=0; m < number_of_Mesh; m++)
	{
		fscanf(file,"%s\n",type);
		printf("%s\n",type);
		if(strcasecmp(type,"Mesh")==0)
		{
			Object o;
			printf("found triangle Mesh\n");
			int material =0;
			parse_mat(file,&material);
			o.material = material;
			int number_of_objects;
			fscanf(file,"%i",&number_of_objects);
			printf("number_of_Mesh: %i\n",number_of_objects);
			for(int i=0; i < number_of_objects; i++)
			{
				fscanf(file,"%s\n",type);
				printf("%s\n",type);
				int j;
				for(j=0;j < 3;j++)
				{
					parse_doubles(file,"pos:",t.v[j].position);
					parse_doubles(file,"nor:",t.v[j].normal);
					parse_doubles(file,"dif:",t.v[j].color_diffuse);
					parse_doubles(file,"spe:",t.v[j].color_specular);
					parse_shi(file,&t.v[j].shininess);
				}

				if(num_triangles == MAX_TRIANGLES)
				{
					printf("too many triangles, you should increase MAX_TRIANGLES!\n");
					exit(0);
				}
				o.triangles[num_triangles++] = t;
			}
			objects[num_objects++] = o;

		}else if(strcasecmp(type,"sphere")==0)
		{
			printf("found sphere\n");
			int material =0;
			parse_mat(file,&material);
			s.material = material;
			parse_doubles(file,"pos:",s.position);
			parse_rad(file,&s.radius);
			parse_doubles(file,"dif:",s.color_diffuse);
			parse_doubles(file,"spe:",s.color_specular);
			parse_shi(file,&s.shininess);

			if(num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
			
		else
		{
			printf("unknown type in scene description:\n%s\n",type);
			exit(0);
		}	
	}
	
	return 0;
}

void setScene(char *sceneName){
	//load scenefile first 
	FILE *scenefile = fopen(sceneName,"r");
	int number_of_item;

	parse_doubles(scenefile,"ambient",ambient_light);
	parse_doubles(scenefile,"eye",eye);
	parse_doubles(scenefile,"look",look);
	parse_doubles(scenefile,"up",up);
	char type[50];
	fscanf(scenefile,"%s\n",type);
	fscanf(scenefile,"%i",&number_of_item);
	int i;
	for(i=0; i < number_of_item; i++)
	{
		Light l;
		parse_doubles(scenefile,"pos",l.position);
		parse_doubles(scenefile,"col",l.color);
		if(num_lights == MAX_LIGHTS)
		{
			printf("too many lights, you should increase MAX_LIGHTS!\n");
			exit(0);
		}
		lights[num_lights++] = l;
	}
}
//--------------------------/obj loader/ ---------------------//
/**
 * @brief parse .obj to our traingle list file 
 * run once and save file in scene
 * 
 * @param objFileName 
 */
void parseObj(char *objFileName ){
// Initialize Loader
	objl::Loader Loader;
	// Load .obj File
	bool loadout = Loader.LoadFile(objFileName);

	// Check to see if it loaded
  
	// If so continue
	if (loadout)
	{
		//output into .txt file 
		std::ofstream file("mesh.scn");
	
		file <<  Loader.LoadedMeshes.size()<<"\n";
		// Go through each loaded mesh and out its contents
		for (int i =0; i < Loader.LoadedMeshes.size(); i++)
		{

			objl::Mesh curMesh = Loader.LoadedMeshes[i];
			file << "Mesh "<<"\n";
			file << "Material: 0" <<"\n";
			file <<  curMesh.Indices.size()/3<<"\n";

			// Copy one of the loaded meshes to be our current mesh
			// Print Indices
			// Go through every 3rd index and print the
			//	triangle that these indices represent
			for (int j = 0; j < curMesh.Indices.size(); j += 3)
			{
				file << "triangle" << "\n";
				for(int i = 0; i<3;i++){
					objl::Vertex v0 = curMesh.Vertices[curMesh.Indices[j+i]];
					file << "pos: " << v0.Position.X<<" "<< v0.Position.Y<<" "<< v0.Position.Z-3<<" "<< "\n";
					file << "nor: " << v0.Normal.X<<" "<< v0.Normal.Y<<" "<< v0.Normal.Z<<" "<< "\n";
					file << "dif: " << curMesh.MeshMaterial.Kd.X<<" "<< curMesh.MeshMaterial.Kd.Y<<" "<< curMesh.MeshMaterial.Kd.Z<<" "<< "\n";
					file << "spe: " << curMesh.MeshMaterial.Ks.X<<" "<< curMesh.MeshMaterial.Ks.Y<<" "<< curMesh.MeshMaterial.Ks.Z<<" "<< "\n";
					file << "shi: " << curMesh.MeshMaterial.Ns<< "\n";
				}
				
			}
			// Leave a space to separate from the next mesh
		}

		// Close File
		file.close();
	}
	// If not output an error
	else
	{
		std::ofstream file("mesh.txt");

		// Output Error
		file << "Failed to Load File. May have failed to find it or it was not an .obj file.\n";

		// Close File
		file.close();
	}
}
void display()
{

}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
        case 'q': case 'Q':
            exit(0);
            break;
        default:
            break;
    }
}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0,WIDTH,0,HEIGHT,1,-1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
	//hack to make it only draw once
	static int once=0;
	if(!once)
	{
		draw_scene();
	}
	once=1;
}

int main (int argc, char ** argv)
{
	if (argc < 2 || argc > 3)
	{  
		printf ("usage: %s <scenefile> <meshfile>\n", argv[0]);
		printf ("usage: %s <objfile>\n", argv[0]);
		
		exit(0);
	}
	if (argc == 2){
		//  (from .obj file to meshfile used as our input)
		printf ("using this application as a obj parser to meshfile");

		parseObj(argv[1]);
		return 0;
	}

	glutInit(&argc,argv);

	setScene(argv[1]);
	loadScene(argv[2]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(WIDTH,HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutKeyboardFunc(keyboard);

	init();
	glutMainLoop();
}
