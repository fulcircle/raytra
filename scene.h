/*
 * scene.h
 *
 *  Created on: Mar 13, 2010
 *      Author: Vikram Sunkavalli
 */

#ifndef SCENE_H_
#define SCENE_H_

#include <iomanip>
#include "camera.h"
#include "objects.h"
#include "support.h"

class Scene {
public:
	Scene(Camera *cam, std::vector<Surface*> *surfaces, Lights *lights, std::vector<double> *opts);
	Render *render();
private:
	array2D<Ray> rays;
	Camera *cam;
	std::vector<Surface*> *surfaces;
	Lights *lights;
	AmbientLight *al;
	static Imf::Rgba black;
	BvhNode *bvhTree;
	int max_depth;
	int samples;
	int shadow_samples;
	int glossy_samples;
	Vector reflect(const Vector& d, const Vector& n);
	bool refract(const Vector &d, const Vector &n, double refract_index1, double refract_index2, Vector *t);
	Imf::Rgba rayColor(Ray *r, double t0, double t1, int depth=0);
};


#endif /* SCENE_H_ */
