#ifndef __MARBLE_H__
#define __MARBLE_H__

#include <math.h>
#include "Vector.h"
#include "PerlinNoise.h"

struct Marble{
   PerlinNoise perlin;
   Vector color1 = Vector(0.1,0.1,0.1);
   Vector color2 = Vector(0.9,0.9,0.89);
   double shininess = 30;

   	Marble() {
	
	}
	Vector getColor(Vector point) {
	double x = point.x;
	double y = point.y;
	double z = point.z;
	double turbulence = 0;
	
	for (int level = 1; level <10; level ++) {
		turbulence +=  std::abs(perlin.noise(
			level * 0.1 * x,
         level * 0.25 * y,
         level * 0.1 * z
		));
	}
	double marbleFrequency =  0.4f;
	double noiseAmplitude = 0.4f;
	double noise = 0.4f * sin((x + y+turbulence * noiseAmplitude) * marbleFrequency ) + 0.8f;

	return color1 * noise + color2 * (1.0f - noise);
	}

	double getShininess() {
		return shininess;
	}
};

#endif
