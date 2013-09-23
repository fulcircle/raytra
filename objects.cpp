/*
 * objects.cpp
 *
 *  Created on: Mar 5, 2010
 *      Author: Vikram Sunkavalli
 * Description: Contains all objects classes: bounding boxes, surfaces (incl. bvh) and lights.
 */

#include "objects.h"


PointLight::PointLight(Point pos, Imf::Rgba intensity) {
	this->pos = pos;
	this->I = intensity;
}

AmbientLight::AmbientLight(Imf::Rgba intensity) {
	this->I = intensity;
	this->pos = Point(0,0,0);
}

AreaLight::AreaLight(Vector edgeVectorA, Vector edgeVectorB, Point pos, Imf::Rgba intensity):PointLight(pos, intensity) {
	this->edgeVectorA = edgeVectorA;
	this->edgeVectorB = edgeVectorB;
	this->N = edgeVectorA^edgeVectorB;
	this->N.normalize();
}

BoundingBox::BoundingBox() {
	minExtent = Point(0,0,0);
	maxExtent = Point(0,0,0);
}

BoundingBox::BoundingBox(Point minExtent, Point maxExtent) {
	this->minExtent = minExtent;
	this->maxExtent = maxExtent;
}

BoundingBox BoundingBox::combine(BoundingBox &a, BoundingBox &b) {
	Point newMinExtent, newMaxExtent;

	if(a.minExtent.x < b.minExtent.x) newMinExtent.x = a.minExtent.x;
	else newMinExtent.x = b.minExtent.x;

	if(a.minExtent.y < b.minExtent.y) newMinExtent.y = a.minExtent.y;
	else newMinExtent.y = b.minExtent.y;

	if(a.minExtent.z < b.minExtent.z) newMinExtent.z = a.minExtent.z;
	else newMinExtent.z = b.minExtent.z;

	if(a.maxExtent.x > b.maxExtent.x) newMaxExtent.x = a.maxExtent.x;
	else newMaxExtent.x = b.maxExtent.x;

	if(a.maxExtent.y > b.maxExtent.y) newMaxExtent.y = a.maxExtent.y;
	else newMaxExtent.y = b.maxExtent.y;

	if(a.maxExtent.z > b.maxExtent.z) newMaxExtent.z = a.maxExtent.z;
	else newMaxExtent.z = b.maxExtent.z;

	BoundingBox newBox(newMinExtent, newMaxExtent);
	return newBox;
}

BoundingBox BoundingBox::getBoundingBox(std::vector<Surface*> *surfs) {
	Point minExtent, maxExtent;
	if(surfs!=NULL && surfs->size() > 0) {
		minExtent = surfs->at(0)->boundingBox.minExtent;
		maxExtent = surfs->at(0)->boundingBox.maxExtent;
		for (std::vector<Surface*>::iterator it_surface = surfs->begin();
				it_surface != surfs->end(); ++it_surface) {
			if((*it_surface)->boundingBox.minExtent.x < minExtent.x) minExtent.x = (*it_surface)->boundingBox.minExtent.x;
			if((*it_surface)->boundingBox.minExtent.y < minExtent.y) minExtent.y = (*it_surface)->boundingBox.minExtent.y;
			if((*it_surface)->boundingBox.minExtent.z < minExtent.z) minExtent.z = (*it_surface)->boundingBox.minExtent.z;

			if((*it_surface)->boundingBox.maxExtent.x > maxExtent.x) maxExtent.x = (*it_surface)->boundingBox.maxExtent.x;
			if((*it_surface)->boundingBox.maxExtent.y > maxExtent.y) maxExtent.y = (*it_surface)->boundingBox.maxExtent.y;
			if((*it_surface)->boundingBox.maxExtent.z > maxExtent.z) maxExtent.z = (*it_surface)->boundingBox.maxExtent.z;
		}

	}
	return BoundingBox(minExtent, maxExtent);
}

bool BoundingBox::hitBox(Ray *r, double t0, double t1) {
	double a, txmin, txmax, tymin, tymax, tzmin, tzmax, tnear, tfar;

	a = 1/r->dir.x;
	if(a>=0) {
		txmin = a*(this->minExtent.x-r->origin.x);
		txmax = a*(this->maxExtent.x-r->origin.x);
	}
	else {
		txmin = a*(this->maxExtent.x-r->origin.x);
		txmax = a*(this->minExtent.x-r->origin.x);
	}

	a = 1/r->dir.y;
	if(a>=0) {
		tymin = a*(this->minExtent.y-r->origin.y);
		tymax = a*(this->maxExtent.y-r->origin.y);
	}
	else {
		tymin = a*(this->maxExtent.y-r->origin.y);
		tymax = a*(this->minExtent.y-r->origin.y);
	}

	a = 1/r->dir.z;
	if(a>=0) {
		tzmin = a*(this->minExtent.z-r->origin.z);
		tzmax = a*(this->maxExtent.z-r->origin.z);
	}
	else {
		tzmin = a*(this->maxExtent.z-r->origin.z);
		tzmax = a*(this->minExtent.z-r->origin.z);
	}

	tnear = std::max(tzmin, std::max(txmin, tymin));
	tfar = std::min(tzmax, std::min(txmax, tymax));

	if(tnear > tfar) return false;

	if ((tnear >= t0 && tnear <= t1) || (tfar >= t0 && tfar <= t1)) { // Ray origin inside volume, interval may span space outside volume, so we have a hit
		return true;
	}
	else if(tnear < t0 && tfar > t1) return true; // Ray and interval are inside bounding volume
	else return false;

}


BvhNode::BvhNode(std::vector<Surface*> surfs, int axis) {
	this->surfs = surfs;
	this->axis = axis;
	this->isNode = true;
	this->typeString = "BvhNode";
	left = NULL;
	right = NULL;
}

void BvhNode::expand() {

	int size = surfs.size();
	if(size==1) {
		left = surfs[0];
		right = NULL;
		boundingBox = surfs[0]->boundingBox;
	}
	else if(size==2) {
		left = surfs[0];
		right = surfs[1];
		boundingBox = BoundingBox::combine(surfs[0]->boundingBox, surfs[1]->boundingBox);
	}
	else {
		BvhNode *nodeToExpand;
		this->split();
		nodeToExpand = new BvhNode(leftSurfs, (axis+1)%3);
		nodeToExpand->expand();
		left = nodeToExpand;

		nodeToExpand = new BvhNode(rightSurfs, (axis+1)%3);
		nodeToExpand->expand();
		right = nodeToExpand;

		boundingBox = BoundingBox::combine(left->boundingBox, right->boundingBox);
	}
}

void BvhNode::split() {
	BoundingBox box = BoundingBox::getBoundingBox(&surfs);
	surfs = Utility::quickSort(surfs, axis);
	int split_index = surfs.size()/2;
	leftSurfs.clear();
	leftSurfs.assign(surfs.begin(), surfs.begin()+split_index);
	rightSurfs.clear();
	rightSurfs.assign(surfs.begin()+split_index, surfs.end());

}

bool BvhNode::hit(Ray *r, double t0, double t1, HitRecord *rec) {
	if(this->boundingBox.hitBox(r, t0, t1)) {
		HitRecord lrec, rrec;
		lrec.num_hits = 0;
		rrec.num_hits = 0;
		bool leftHit=false, rightHit=false;
		if(left!=NULL) leftHit = left->hit(r, t0, t1, &lrec);
		if(right!=NULL) rightHit = right->hit(r, t0, t1, &rrec);

		if(leftHit && rightHit) {
			if (lrec.hits[0] < rrec.hits[0]) *rec = lrec;
			else *rec = rrec;
			return true;
		}
		else if(leftHit) {
			*rec=lrec;
			return true;
		}
		else if(rightHit) {
			*rec = rrec;
			return true;
		}
		else return false;

	}
	else return false;
}

Vector BvhNode::getNormal(Point p) {
	return Vector();
}

BvhNode::~BvhNode() {
	delete left;
	delete right;
}

Surface::Surface() {
	this->typeString = "surface";
	this->setMaterial(Imf::Rgba(0,0,0,0), Imf::Rgba(0,0,0,0), 0, Imf::Rgba(0,0,0,0), 0, Imf::Rgba(0,0,0,0),0);
	this->isNode = false;
}

Surface::Surface(Point origin) {
	this->typeString = "surface";
	this->setMaterial(Imf::Rgba(0,0,0,0), Imf::Rgba(0,0,0,0), 0, Imf::Rgba(0,0,0,0), 0, Imf::Rgba(0,0,0,0), 0);
	this->origin = origin;
	this->isNode = false;
}

void Surface::setMaterial(Imf::Rgba diffuse, Imf::Rgba specular, double phong_exp, Imf::Rgba ideal, double refract_index, Imf::Rgba absorption_coefficient, double glossy) {
	this->diffuse = diffuse;
	this->specular = specular;
	this->phong_exp = phong_exp;
	this->ideal = ideal;
	this->refract_index = refract_index;
	this->absorption_coefficient = absorption_coefficient;
	this->glossy = glossy;
}

void Surface::setRefract(double refract_index, Imf::Rgba absorption_coefficient) {
	this->refract_index = refract_index;
	this->absorption_coefficient = absorption_coefficient;
}

void Surface::setGlossy(double glossy) {
	this->glossy = glossy;
}

void Surface::copyHitRecordInfo(HitRecord *rec, std::vector<double> *hits, Ray *r) {
	if(hits->size() < 1) {
		rec->num_hits = 0;
		return;
	}
	for(uint i=0; i<hits->size(); ++i) {
		rec->hits.push_back(hits->at(i));
		++rec->num_hits;
	}
	rec->ka = this->diffuse;
	rec->kd = this->diffuse;
	rec->ks = this->specular;
	rec->km = this->ideal;
	rec->N = this->getNormal(r->getPoint(hits->at(0)));
	rec->phong_exp = this->phong_exp;
	rec->refract_index = this->refract_index;
	rec->kc = this->absorption_coefficient;
	rec->glossy = this->glossy;
	rec->typeString = this->typeString;
}


Sphere::Sphere(Point center, double radius):Surface(center) {
	this->radius = radius;
	sqrdRadius = radius*radius;
	typeString = "sphere";
	boundingBox.minExtent = Point(center.x-radius, center.y-radius, center.z-radius);
	boundingBox.maxExtent = Point(center.x+radius, center.y+radius, center.z+radius);
}

Vector Sphere::getNormal(Point p) {
	Vector N = (p-origin)/radius;
	return N;
}

double Sphere::getRadius() {
	return radius;
}

double Sphere::getSqrdRadius() {
	return sqrdRadius;
}


//RAY-SPHERE INTERSECTION
// Very fast ray/sphere intersection testing
// Derived from: http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
// Further optimization done as explained here: http://www.devmaster.net/wiki/Talk:Ray-sphere_intersection
bool Sphere::hit(Ray *ray, double t0, double t1, HitRecord *rec) {
	ray->normalize();
	std::vector<double> hits;

	// First check if ray origin is inside sphere (generally needed for refraction)
	Vector dist = this->origin - ray->origin;
	double distDotdist = dist*dist;
	double diff = this->sqrdRadius-distDotdist;

	if (diff >= EPSILON) { // ray origin inside of the sphere
		double distDotdir = dist*ray->dir;
		double hcSquared = (diff) + distDotdir*distDotdir;
		double h0 =  distDotdir + sqrt(hcSquared);
		hits.push_back(h0);
		this->copyHitRecordInfo(rec, &hits, ray);
		return 1;
	}

	// Assume ray is normalized so A = 1;
	dist = -dist;
	double B = dist*ray->dir;
	double C = dist*dist-this->sqrdRadius;
	double D = B*B-C;

	if(fabs(D)<EPSILON) {
		double h0 = -B;
		if(h0 >= t0 && h0 <= t1) {
			hits.push_back(h0);
			this->copyHitRecordInfo(rec, &hits, ray);
			return 1;
		}
		else return 0;
	}

	if (D < 0) return 0;

	double sqrtD = sqrt(D);
	double h0 = -B - sqrtD;
	double h1 = -B + sqrtD;
	if(h0 >= t0 && h0 <= t1) hits.push_back(h0);
	if(h1 >= t0 && h1 <= t1) hits.push_back(h1);
	if(hits.size() > 0) {
		this->copyHitRecordInfo(rec, &hits, ray);
		return 1;
	}
	else return 0;
}

Triangle::Triangle(Point a, Point b, Point c):Surface() {
	this->a = a;
	this->b = b;
	this->c = c;
	vertices.clear();
	vertices.push_back(a); vertices.push_back(b); vertices.push_back(c);
	double oneThird = 1.0/3.0;
	double xcenter = oneThird*(a.x+b.x+c.x);
	double ycenter = oneThird*(a.y+b.y+c.y);
	double zcenter = oneThird*(a.z+b.z+c.z);
	this->origin = Point(xcenter, ycenter, zcenter); // Centroid
	typeString = "triangle";
	N = (b-a)^(c-a);
	N.normalize();

	boundingBox.minExtent = this->a;
	boundingBox.maxExtent = this->a;

	for(int i=1; i<3; ++i) {
		if (vertices[i].x < boundingBox.minExtent.x) boundingBox.minExtent.x = vertices[i].x;
		if (vertices[i].y < boundingBox.minExtent.y) boundingBox.minExtent.y = vertices[i].y;
		if (vertices[i].z < boundingBox.minExtent.z) boundingBox.minExtent.z = vertices[i].z;

		if (vertices[i].x > boundingBox.maxExtent.x) boundingBox.maxExtent.x = vertices[i].x;
		if (vertices[i].y > boundingBox.maxExtent.y) boundingBox.maxExtent.y = vertices[i].y;
		if (vertices[i].z > boundingBox.maxExtent.z) boundingBox.maxExtent.z = vertices[i].z;
	}
}

// Fast ray/triangle intersection
// From: http://www.devmaster.net/wiki/Ray-triangle_intersection
bool Triangle::hit(Ray *ray, double t0, double t1, HitRecord *rec) {
	std::vector<double> hits;
	Vector N = this->getNormal(Point());
	ray->normalize();
	Vector D = ray->dir;
	Point o = ray->origin;
	Point a = this->a;
	double t = -((o-a)*N)/(D*N); // First get distance from ray origin to plane of triangle

	if(t < 0) return false; // Plane is behind ray, so no intersection

	// Now find point P where ray intersects plane
	Point p = o+(t*D);

	//Check if distance is within interval
	if(t < t0 || t > t1) return false;

	//Now check if p lies in the triangle
	//First we will project the triangle to its dominant axis in order to reduce computation (R^3->R^2)

	//Find the dominant axis and set the indices to the corresponding planes
	int index0, index1;
	if(fabs(N.z) > fabs(N.x) && fabs(N.z) > fabs(N.y)) { index0 = 0; index1 = 1; }
	else if(fabs(N.y) > fabs(N.x)) { index0 = 0; index1 = 2; }
	else { index0 = 1; index1 = 2; }

	//Now compute u,v in R^2
	Vector B = this->b-this->a;
	Vector C = this->c-this->a;
	Vector P = p-a;

	double u = (P[index1]*C[index0]-P[index0]*C[index1])/(B[index1]*C[index0]-B[index0]*C[index1]);
	if(u < 0 ) return false;
	double v = (P[index1]*B[index0]-P[index0]*B[index1])/(C[index1]*B[index0]-C[index0]*B[index1]);
	if(v < 0) return false;

	if(u+v > 1) return false;

	hits.push_back(t);
	this->copyHitRecordInfo(rec, &hits, ray);
	return true;
}

Vector Triangle::getNormal(Point p) {
	return N;
}

Plane::Plane(Vector N, double d):Surface() {
	N.normalize();
	this->N = N;
	this->d = d;
	double infin = std::numeric_limits<double>::infinity();
	this->boundingBox.minExtent = Point(-infin, -infin, -infin);
	this->boundingBox.maxExtent = Point(infin, infin, infin);
	typeString = "plane";
}

bool Plane::hit(Ray *ray, double t0, double t1, HitRecord *rec) {

	std::vector<double> hits;
	ray->normalize();

	double Vd = ray->dir*this->N;
	if (Vd>=0) return false;

	double V0 = -(this->N*ray->origin+this->d);

	double t = V0/Vd;

	if(t < t0 || t > t1) return false;

	hits.push_back(t);
	this->copyHitRecordInfo(rec, &hits, ray);

	return true;
}

Vector Plane::getNormal(Point p) {
	return N;
}
