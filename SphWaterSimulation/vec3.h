#ifndef VEC3_H_
#define VEC3_H_

#include <iostream>
#include <cassert>
#include <cmath>

class Vec3{
public:
	double x, y, z;
	Vec3(double x = 0., double y = 0., double z = 0.): x(x), y(y), z(z){}
	Vec3(const Vec3& a) = default;
	Vec3 operator+(const Vec3& a) const{ return Vec3(a.x + x, a.y + y, a.z + z); }
	Vec3 operator-(const Vec3& a) const{ return Vec3(x - a.x, y - a.y, z - a.z); }
	Vec3& operator+=(const Vec3& a){ return *this = *this + a; }
	Vec3& operator-=(const Vec3& a){ return *this = *this - a; }
	Vec3& operator*=(double p){ return *this = *this * p; }
	Vec3& operator/=(double p){ return *this = *this / p; }

	Vec3 operator+(double a) const{ return Vec3(a + x, a + y, a + z); }
	Vec3 operator-(double a) const{ return Vec3(x - a, y - a, z - a); }
	Vec3 operator*(double a) const{ return Vec3(a * x, a * y, a * z); }
	Vec3 operator/(double a) const{
		assert(a != 0. && "Divided by zero!");
		return Vec3(x / a, y / a, z / a);
	}
	Vec3 operator-() const{ return Vec3(-x, -y, -z); }

	friend std::ostream& operator<<(std::ostream& os, const Vec3& a){
		os << "Vec3(" << a.x << ", " << a.y << ", " << a.z << ")";
		return os;
	}

	bool operator==(const Vec3& a) const{ return x == a.x && y == a.y && z == a.z; }
	bool operator!=(const Vec3& a) const{ return x != a.x || y != a.y || z != a.z; }

	double maxCoord() const{
		return x > y && x > z ? x : y > z ? y : z;
	}
	double maxAbsCoord() const{
		return abs(x) > abs(y) && abs(x) > abs(z) ? abs(x) : abs(y) > abs(z) ? abs(y) : abs(z);
	}

	double dot(const Vec3& a) const{
		return (a.x * x + a.y * y + a.z * z);
	}

	// Perform element-wise multiplication.
	Vec3 multiply(const Vec3& a) const{ return Vec3(x * a.x, y * a.y, z * a.z); }
	double norm() const{
		return sqrt(x * x + y * y + z * z);
	}

	// Return the square norm of vector || x || ^ 2.
	double sqrNorm()const{
		return x * x + y * y + z * z;
	}
	Vec3 normalized(){
		double n_norm = this->norm();
		if(n_norm == 0.) return Vec3(0., 0., 0.);
		return Vec3(x / n_norm, y / n_norm, z / n_norm);
	}
	Vec3 crossProduct(const Vec3& a)const{
		return Vec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
	}

	// Return a new Vec3 with each coordinate clipped in range [low, high].
	Vec3 clip(double low = 0, double high = 1) const{
		return Vec3(x > high ? high : x < low ? low : x,
			y > high ? high : y < low ? low : y,
			z > high ? high : z < low ? low : z);
	}

	double getx() const{ return x; }
	double gety() const{ return y; }
	double getz() const{ return z; }

	void setx(double x_){ x = x_; }
	void sety(double y_){ y = y_; }
	void setz(double z_){ z = z_; }
};

#endif