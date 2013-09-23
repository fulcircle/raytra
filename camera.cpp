/*
 * camera.cpp
 *
 *  Created on: Mar 13, 2010
 *      Author: Vikram Sunkavalli
 */

#include "camera.h"

Vector Camera::UP(0,1,0);

Camera::Camera() {
	d = 0;
	width = 0;
	height = 0;
	l = r = t = b = 0;
}

Camera::Camera(Point e, Vector dir, double lr, double tb, double fl, int width, int height) {
	constructBasis(dir);
	this->e = e;
	this->width = width;
	this->height = height;
	d = fl;
	l = -lr/2.0; r = lr/2.0; b = -tb/2.0; t = tb/2.0;


}

void Camera::getBasis(Vector *ub, Vector *vb, Vector *wb) {
	*ub = u;
	*vb = v;
	*wb = w;
}

void Camera::constructBasis(Vector& dir) {
	dir.normalize();
	u = dir^UP;
	v = u^dir;
	w = -dir;
	u.normalize(); v.normalize(); w.normalize();
}
