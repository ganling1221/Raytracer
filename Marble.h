#ifndef __MARBLE_H__
#define __MARBLE_H__

#include <math.h>
#include "Vector.h"
#include "PerlinNoise.h"

struct Marble{
   PerlinNoise perlin;
   Vector color1 = Vector(0.4,0.3,0.2);
   Vector color2 = Vector(0.6,0.7,0.8);
   double shininess = 50;
   double reflectivity = 0.15;

   	Marble() {
	
	}
	Vector getColor(Vector point) {
	double x = point.x;
	double y = point.y;
	double z = point.z;
	double noiseCoef = 0;
	
	for (int level = 1; level < 5; level ++) {
		noiseCoef +=  pow(2 ,-level) * std::abs(perlin.noise(
			pow(2 ,level) * x,
			pow(2 ,level) * y,
			pow(2 ,level) * z
		));
	}
	noiseCoef = 0.5f * sinf((x + y) * 0.05f + noiseCoef) + 0.5f;

	return color2 * noiseCoef + color1 * (1.0f - noiseCoef);
	}

	double getShininess() {
		return shininess;
	}
};

#endif
