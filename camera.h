/*
 * camera.h
 *
 *  Created on: Mar 13, 2010
 *      Author: vik
 */

#ifndef CAMERA_H_
#define CAMERA_H_

#include "globals.h"
#include "support.h"
#include "objects.h"

class Camera {
public:
	Point e;
	Vector u, v, w;

	double d; //focal length
	double width, height;
	double l, r, t, b;

	static Vector UP;

	Camera();
	Camera(Point e, Vector dir, double lr, double tb, double fl, int width, int height);
	void getBasis(Vector *ub, Vector *vb, Vector *wb);
private:
	array2D<Ray> rays;
	void constructBasis(Vector& dir);
};

#endif /* CAMERA_H_ */
