/*
 * objects.h
 *
 *  Created on: Mar 5, 2010
 *      Author: Vikram Sunkavalli
 * Description: Contains all objects classes: bounding boxes, surfaces (incl. bvh) and lights.
 */

#ifndef OBJECTS_H_
#define OBJECTS_H_

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>
#include <iostream>
#include "support.h"
#include "globals.h"

typedef struct _hit_record {
	int num_hits;
	std::vector<double> hits;
	Vector N;
	Imf::Rgba ka; // Ambient component
	Imf::Rgba kd; // Diffuse component
	Imf::Rgba ks; // Specular component
	Imf::Rgba km; // Ideal component
	Imf::Rgba kc; // Absorption coefficient
	double phong_exp;
	double refract_index;
	double glossy;
	std::string typeString;
	int id;
} HitRecord;

class PointLight;
class AmbientLight;
class AreaLight;

typedef struct _lights {
	std::vector<PointLight*> pointLights;
	std::vector<AmbientLight*> ambientLights;
	std::vector<AreaLight*> areaLights;
} Lights;

class BoundingBox {
public:
	Point minExtent;
	Point maxExtent;

	BoundingBox();
	BoundingBox(Point minExtent, Point maxExtent);
	static BoundingBox combine(BoundingBox &a, BoundingBox &b);
	static BoundingBox getBoundingBox(std::vector<Surface*> *surfs);
	bool hitBox(Ray *r, double t0, double t1);
};

class Surface {
public:
	std::string typeString;
	Point origin;
	Imf::Rgba ideal;
	Imf::Rgba diffuse;
	Imf::Rgba specular;
	double phong_exp;
	double refract_index;
	Imf::Rgba absorption_coefficient;
	double glossy;
	BoundingBox boundingBox;
	bool isNode;
	Surface();
	Surface(Point origin);
	void setMaterial(Imf::Rgba diffuse, Imf::Rgba specular, double phong_exp, Imf::Rgba ideal, double refract_index=0, Imf::Rgba absorption_coefficient=Imf::Rgba(0,0,0), double glossy=0);
	void setRefract(double refract_index, Imf::Rgba absorption_coefficient);
	void setGlossy(double glossy);
	void copyHitRecordInfo(HitRecord *rec, std::vector<double> *hits, Ray *r);
	virtual bool hit(Ray *ray, double t0, double t1, HitRecord *rec) = 0;
	virtual Vector getNormal(Point p)=0;

};

class BvhNode:public Surface {
public:
	Surface *left;
	Surface *right;
	BvhNode(std::vector<Surface*> surfs, int axis);
	void expand();
	virtual bool hit(Ray *r, double t0, double t1, HitRecord *rec);
	virtual Vector getNormal(Point p);
	~BvhNode();
private:
	std::vector<Surface*> surfs;
	std::vector<Surface*> leftSurfs;
	std::vector<Surface*> rightSurfs;
	int axis;
	void split();
};

class Sphere:public Surface {
public:
	Sphere(Point center, double radius);
	bool hit(Ray *ray, double t0, double t1, HitRecord *rec);
	Vector getNormal(Point p);
	double getRadius();
	double getSqrdRadius();
private:
	double radius;
	double sqrdRadius;
};

class Triangle:public Surface {
public:
	Point a;
	Point b;
	Point c;
	std::vector<Point> vertices;
	Vector N; // Normal

	Triangle(Point a, Point b, Point c);
	bool hit(Ray *ray, double t0, double t1, HitRecord *rec);
	Vector getNormal(Point p);
};

class Plane:public Surface {
public:
	Vector N; // Normal
	double d; // Scalar value

	Plane(Vector N, double d);
	bool hit(Ray *ray, double t0, double t1, HitRecord *rec);
	Vector getNormal(Point p);
};

class Light {
public:
	Point pos; // Light position
	Imf::Rgba I; //Intensity
};

class PointLight:public Light {
public:
	PointLight(Point pos, Imf::Rgba intensity);
};

class AmbientLight:public Light {
public:
	AmbientLight(Imf::Rgba intensity);
};

class AreaLight:public PointLight {
public:
	Vector edgeVectorA;
	Vector edgeVectorB;
	Vector N;

	AreaLight(Vector edgeVectorA, Vector edgeVectorB, Point pos, Imf::Rgba intensity);
};

#endif /* OBJECTS_H_ */
