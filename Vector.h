#ifndef __VECTOR_H__
#define __VECTOR_H__
#include <math.h>

struct Vector
{
	/* data */
	double x;
	double y;
	double z;

	/* functions */
	Vector() : x(0.0), y(0.0), z(0.0)
	{}

	Vector(double i, double j, double k) : x(i), y(j), z(k)
	{}

	Vector operator+(const double c)
	{
		return Vector(x + c, y + c, z + c);
	}

	Vector operator+(const Vector v)
	{
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	Vector operator-(const Vector v)
	{
		return Vector(x - v.x, y - v.y, z - v.z);
	}

	Vector operator*(const double c)
	{
		return Vector(c * x, c * y, c * z);
	}

	Vector operator*(const Vector &v)
	{
		return Vector(x * v.x, y * v.y, z * v.z);
	}

	Vector operator/(const Vector &v)
	{
		return Vector(x / v.x, y / v.y, z / v.z);
	}

	bool operator==(const Vector v)
	{
		return x == v.x && y == v.y && z == v.z;
	}

	void operator=(const Vector &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	void operator+=(const Vector &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	void operator/=(const double c)
	{
		x /= c;
		y /= c;
		z /= c;
	}

	void normalize()
	{
		double len = sqrt(x*x + y*y + z*z);
		x /= len;
		y /= len;
		z /= len;
	}

	double sum()
	{
		return (x + y + z);
	}

	double magnitude()
	{
		return sqrt(x*x + y*y + z*z);
	}
};


#endif
