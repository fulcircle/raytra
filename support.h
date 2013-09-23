/*
 * support.h
 *
 *  Created on: Mar 3, 2010
 *      Author: Vikram Sunkavalli
 */


#ifndef SUPPORT_H_
#define SUPPORT_H_

#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <vector>

class Vector;
class Surface;

template <class T>
class array2D {

private:
	std::vector<std::vector<T> > vect2D;

public:
	void resizeErase(int height, int width) {
		vect2D.clear();
		vect2D.resize(height);
		for(int y=0; y<height; ++y) {
			vect2D[y].clear();
			vect2D[y].resize(width);
		}
	}
	std::vector<T>& operator[](const int nIndex) {
		return vect2D[nIndex];
	}
};

class Point {
public:
	double x, y, z;
	Point();
	Point(double x, double y, double z);

	friend Vector operator-(const Point& p1, const Point& p2);
	friend Vector operator+(const Vector& v, const Point& p);
	friend Vector operator+(const Point& p, const Vector& v);
	friend double operator*(const Vector& v, const Point &p);
	friend double operator*(const Point &p, const Vector& v);
	Point& operator=(const Point& p2);
	double operator[](const int nIndex);
	static double distance(const Point& p1, const Point& p2);
	void setPoint(double x, double y, double z);
	friend std::ostream& operator<<(std::ostream& output, const Point& p);
};

class Vector:public Point {
public:
	Vector();
	Vector(double x, double y, double z);

	friend Vector operator+(const Vector& v1, const Vector& v2);
	friend Vector operator-(const Vector& v1, const Vector& v2);
	friend double operator*(const Vector& v1, const Vector& v2);
	friend Vector operator^(const Vector& v1, const Vector& v2);
	friend Vector operator*(const Vector& v, const double c);
	friend Vector operator*(const double c, const Vector& v);
	friend Vector operator/(const double c, const Vector& v);
	friend Vector operator/(const Vector& v, const double c);
	friend Vector& operator+=(Vector& v1, const Vector& v2);
	friend Vector& operator-=(Vector& v1, const Vector& v2);
	friend Vector& operator*=(Vector& v, const double c);
	friend Vector& operator/=(Vector& v, const double c);
	friend Vector operator-(const Vector& v);
	friend Vector operator+(const Vector& v, const Point& p);
	friend Vector operator+(const Point& p, const Vector& v);
	friend double operator*(const Vector& v, const Point &p);
	friend double operator*(const Point &p, const Vector& v);
	Vector& operator=(const Vector& p2);
	void setVector(double x, double y, double z);
	double length();
	void normalize();
};

class Ray {
public:
	Point origin;
	Vector dir;
	Ray();
	Ray(Point origin, Vector dir);
	Ray(Point origin, Point p);
	void normalize();
	Point getPoint(double t);
};

class Render {
friend class Utility;
public:
	Render(int height, int width);
	void resize(int height, int width);
	Imf::Rgba *getPixel(int i, int j);
	void setPixel(int i, int j, const Imf::Rgba *color);

private:
	Imf::Array2D<Imf::Rgba> pixels;
	int width, height;
};

class Utility {
public:
	static void renderToFile(const char fileName[], const Render* const render);
	static void printArray2D(const Imf::Array2D<int> &array2D, int height, int width);
	static std::vector<Surface*> quickSort(std::vector<Surface*> surfs, int axis);
private:
	static void writeRgba (const char fileName[], const Imf::Rgba *pixels, int width, int height);
};

#endif /* SUPPORT_H_ */
