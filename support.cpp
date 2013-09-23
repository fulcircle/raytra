/*
 * support.c
 *
 *  Created on: Mar 3, 2010
 *      Author: Vikram Sunkavalli
 */

#include <iomanip>
#include <cmath>
#include <iostream>
#include "support.h"
#include "objects.h"


//POINT
Point::Point() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
}
Point::Point(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

void Point::setPoint(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

double Point::distance(const Point& p1, const Point& p2) {
		double dx = p1.x - p2.x;
		double dy = p1.y - p2.y;
		double dz = p1.z - p2.z;
		return sqrt(dx*dx + dy*dy + dz*dz);

}

double Point::operator[](const int nIndex) {
	if(nIndex==0) return x;
	if(nIndex==1) return y;
	if(nIndex==2) return z;
	return NULL;
}

std::ostream& operator<<(std::ostream& outstream, const Point& p) {
	outstream.precision(2);
	outstream << "[" << std::fixed << p.x << " " << p.y << " " << p.z << "]";
	return outstream;
}

//OVERLOADED POINT OPERATIONS
Vector operator-(const Point& p1, const Point& p2) {
	Vector n(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z);
	return n;
}

Point& Point::operator=(const Point& p2) {
	this->x = p2.x;
	this->y = p2.y;
	this->z = p2.z;
	return *this;
}


//VECTOR
Vector::Vector():Point() {}

Vector::Vector(double x, double y, double z):Point(x, y, z) {}

void Vector::setVector(double x, double y, double z) {
		setPoint(x,y,z);
}

double Vector::length() {
	return sqrt(x*x+y*y+z*z);
}
void Vector::normalize() {
	double len = length();
	if(len==0 || len==1) return;
	x/=len;
	y/=len;
	z/=len;
}

//OVERLOADED VECTOR OPERATIONS
Vector& Vector::operator=(const Vector& p2) {
	this->x = p2.x;
	this->y = p2.y;
	this->z = p2.z;
	return *this;
}

Vector operator+(const Vector& v1, const Vector& v2) {
	Vector n(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
	return n;
}

Vector operator-(const Vector& v1, const Vector& v2) {
	Vector n(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
	return n;
}

Vector operator*(const Vector& v, const double c) {
	Vector n(c*v.x, c*v.y, c*v.z);
	return n;
}

Vector operator*(const double c, const Vector& v) {
	Vector n(c*v.x, c*v.y, c*v.z);
	return n;
}

Vector operator/(const Vector& v, const double c) {
	Vector n(v.x/c, v.y/c, v.z/c);
	return n;
}

Vector operator/(const double c, const Vector& v) {
	Vector n(v.x/c, v.y/c, v.z/c);
	return n;
}

// Dot product
double operator*(const Vector& v1, const Vector& v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

// Cross product
Vector operator^(const Vector& v1, const Vector& v2) {
	Vector n(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
	return n;
}

Vector& operator+=(Vector& v1, const Vector& v2) {
	v1.x+=v2.x;
	v1.y+=v2.y;
	v1.z+=v2.z;
	return v1;
}
Vector& operator-=(Vector& v1, const Vector& v2) {
	v1.x-=v2.x;
	v1.y-=v2.y;
	v1.z-=v2.z;
	return v1;
}

Vector& operator*=(Vector& v, const double c) {
	v.x*=c;
	v.y*=c;
	v.z*=c;
	return v;
}

Vector& operator/=(Vector& v, const double c) {
	v.x/=c;
	v.y/=c;
	v.z/=c;
	return v;
}

Vector operator-(const Vector& v) {
	Vector n(-v.x, -v.y, -v.z);
	return n;
}


//Treat a point as a vector for vector addition/multiplications between point and vector
Vector operator+(const Vector& v, const Point& p) {
	Vector n(v.x+p.x, v.y+p.y, v.z+p.z);
	return n;
}

Vector operator+(const Point& p, const Vector& v) {
	Vector n(v.x+p.x, v.y+p.y, v.z+p.z);
	return n;
}

double operator*(const Vector& v, const Point &p) {
	return v.x*p.x+v.y*p.y+v.z*p.z;
}

double operator*(const Point &p, const Vector &v) {
	return v.x*p.x+v.y*p.y+v.z*p.z;
}


//RAY
Ray::Ray() {}
Ray::Ray(Point origin, Vector dir) {
	this->origin = origin;
	this->dir = dir;
	this->dir.normalize();
}
Ray::Ray(Point origin, Point p) {
	this->origin = origin;
	this->dir = p - origin;
	this->dir.normalize();
}

void Ray::normalize() {
	this->dir.normalize();
}

Point Ray::getPoint(double t) {
	return this->origin+t*this->dir;

}

Render::Render(int height, int width) {
	resize(height, width);

}

void Render::resize(int height, int width) {
	this->width = width;
	this->height = height;
	pixels.resizeErase(height, width);
}

void Render::setPixel(int i, int j, const Imf::Rgba *color) {
	pixels[i][j] = *color;

}

Imf::Rgba *Render::getPixel(int i, int j) {
	return &pixels[i][j];
}

void Utility::renderToFile (const char fileName[], const Render *render) {
	writeRgba(fileName, &render->pixels[0][0], render->width, render->height);
}

void Utility::writeRgba (const char fileName[], const Imf::Rgba *pixels, int width, int height) {

	Imf::RgbaOutputFile file (fileName, width, height, Imf::WRITE_RGBA);
	file.setFrameBuffer (pixels, 1, width);
	file.writePixels (height);
}

void Utility::printArray2D(const Imf::Array2D<int> &array2D, int height, int width) {
	std::cout.precision(1);
	for(int y=0; y < height; ++y) {
		std::cout <<"| ";
		for (int x=0; x < width; ++x) {
			std::cout << std::fixed << array2D[y][x] << " ";
		}
		std::cout <<"|" << std::endl;
	}
}

std::vector<Surface*> Utility::quickSort(std::vector<Surface*> surfs, int axis) {
	std::vector<Surface*> less, greater, combined;
	if(surfs.size()<=1) {
		return surfs;
	}

	int pivot_index = surfs.size()/2;
	Surface *pivot = surfs[pivot_index];
	surfs.erase(surfs.begin()+pivot_index);


	for(uint i=0; i < surfs.size(); ++i) {
		if(surfs[i]->origin[axis]<=pivot->origin[axis])
			less.push_back(surfs[i]);
		else greater.push_back(surfs[i]);
	}
	std::vector<Surface*> sortedL = Utility::quickSort(less, axis);
	std::vector<Surface*> sortedG = Utility::quickSort(greater, axis);
	combined.insert(combined.begin(), sortedL.begin(), sortedL.end());
	combined.push_back(pivot);
	combined.insert(combined.end(), sortedG.begin(), sortedG.end());
	return combined;

}
