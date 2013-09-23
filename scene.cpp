/*
 * scene.cpp
 *
 *  Created on: Mar 13, 2010
 *      Author: Vikram Sunkavalli
 * Description: Main renderer.
 */

#include "scene.h"

Imf::Rgba Scene::black = Imf::Rgba(0,0,0,1);

Scene::Scene(Camera *cam, std::vector<Surface*> *surfaces, Lights *lights, std::vector<double> *opts) {
	this->cam = cam;
	this->surfaces = surfaces;
	this->lights = lights;
	if(opts->size()==4) {
		this->max_depth = opts->at(0);
		this->samples = opts->at(1);
		if(this->samples < 1) this->samples = 1;
		this->shadow_samples = opts->at(2);
		if(this->shadow_samples < 1) this->shadow_samples = 1;
		this->glossy_samples = opts->at(3);
		if(this->glossy_samples < 1) this->glossy_samples = 1;
	}
	else {
		this->max_depth = DEFAULT_MAX_DEPTH;
		this->samples = DEFAULT_SAMPLES;
		this->shadow_samples = DEFAULT_SHADOW_SAMPLES;
		this->glossy_samples = DEFAULT_GLOSSY_SAMPLES;
	}

	this->bvhTree = new BvhNode(*surfaces, X_AXIS);
	this->bvhTree->expand();
}

Render *Scene::render() {


	int height = cam->height;
	int width = cam->width;

	double t = cam->t;
	double b = cam->b;

	double l = cam->l;
	double r = cam->r;
	int rminusl = r-l;
	int bminust = b-t;


	Render *aRender = new Render(cam->height, cam->width);

	if(lights->ambientLights.size()>0)
		al = lights->ambientLights[0];
	else {
		al = new AmbientLight(black);
	}

	Ray view_ray;
	double infinity = std::numeric_limits<double>::infinity();
	int pq_samples = sqrt(this->samples);
	int total_pixels = height*width;

	srand(time(NULL));
	double rand_n;

	int rendered_pixels = 0;
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			double percent = ((double)rendered_pixels/(double)total_pixels)*100;
			std::cout << "Rendered: " <<  std::setprecision(5) <<  percent << std::endl;
			Imf::Rgba c=black;
			for(int p=0; p < pq_samples; ++p) {
				for(int q=0; q < pq_samples; ++q) {
					if(pq_samples > 1)
						rand_n = ((double) rand()) / RAND_MAX;
					else
						rand_n = 0.5;
					double u = l + rminusl*(x+(p+rand_n)/pq_samples)/width;
					double v = t + bminust*(y+(q+rand_n)/pq_samples)/height;
					view_ray.origin = cam->e;
					view_ray.dir = u*cam->u+v*cam->v-cam->d*cam->w;
					view_ray.normalize();
					Imf::Rgba rayc = rayColor(&view_ray, 0, infinity);
					c.r += rayc.r; c.g += rayc.g; c.b += rayc.b;
				}
			}
			c.r = c.r/(this->samples);
			c.g = c.g/(this->samples);
			c.b = c.b/(this->samples);
			aRender->setPixel(y,x,&c);
			++rendered_pixels;
		}
	}

	return aRender;
}

Imf::Rgba Scene::rayColor(Ray *r, double t0, double t1, int depth) {

	++depth;
	HitRecord rec, srec;
	Imf::Rgba ray_color=black;

	rec.hits.resize(0);
	rec.num_hits = 0;
	bool result = bvhTree->hit(r, t0, t1, &rec);

	if(result) {

		double infinity = std::numeric_limits<double>::infinity();
		double thit = rec.hits[0];
		Point p = r->getPoint(thit);

		//Lighting
		double ambr = rec.ka.r*al->I.r;
		double ambg = rec.ka.g*al->I.g;
		double ambb = rec.ka.b*al->I.b;

		double lamr=0, lamg=0, lamb=0;
		double bpr=0, bpg=0, bpb=0;
		double reflectr=0, reflectg=0, reflectb=0;
		double refractr=0, refractg=0, refractb=0;
		double kr=1, kg=1, kb=1;

		Vector v = p-r->origin;
		v.normalize();

		// Point Light shadows
		for(std::vector<PointLight*>::iterator it_pointLight = lights->pointLights.begin();
				it_pointLight != lights->pointLights.end(); ++it_pointLight) {
			PointLight *pl = *it_pointLight;
			Vector l = pl->pos-p;
			double tmax = l.length();
			l.normalize();
			Ray *sr = new Ray(p, l); // Shadow ray
			if(!bvhTree->hit(sr, EPSILON, tmax, &srec)) {
				Vector h = -v+l;
				h.normalize();
				double ndotl = rec.N*l;
				double ndoth = rec.N*h;
				double maxndotl = std::max(0.0,ndotl);
				double maxndothp = pow(std::max(0.0,ndoth),rec.phong_exp);

				lamr += rec.kd.r*pl->I.r*maxndotl;
				lamg += rec.kd.g*pl->I.g*maxndotl;
				lamb += rec.kd.b*pl->I.b*maxndotl;

				bpr += rec.ks.r*pl->I.r*maxndothp;
				bpg += rec.ks.g*pl->I.g*maxndothp;
				bpb += rec.ks.b*pl->I.b*maxndothp;

			}
			delete sr;
		}

		// Area light shadows
		int edge_samples = sqrt(this->shadow_samples);
		for(std::vector<AreaLight*>::iterator it_areaLight = lights->areaLights.begin();
				it_areaLight != lights->areaLights.end(); ++it_areaLight) {
			AreaLight *al = (*it_areaLight);
			double intensityr = al->I.r/this->shadow_samples;
			double intensityg = al->I.g/this->shadow_samples;
			double intensityb = al->I.b/this->shadow_samples;
			double la = al->edgeVectorA.length();
			double lb = al->edgeVectorB.length();
			double step_a = la/edge_samples;
			double step_b = lb/edge_samples;
			for(double sx=0; sx < la; sx+=step_a) {
				for(double sy=0; sy < lb; sy+=step_b) {
					double rand_a = ((double) rand()) / RAND_MAX;
					double rand_b = ((double) rand()) / RAND_MAX;

					double a = sx+(step_a*rand_a);
					double b = sy+(step_b*rand_b);

					Point pos = al->pos+a*al->edgeVectorA+b*al->edgeVectorB;
					Vector l = pos-p;
					double tmax = l.length();
					l.normalize();

					Ray *sr = new Ray(p, l); // Shadow ray
					if(!bvhTree->hit(sr, EPSILON, tmax, &srec)) {
						Vector h = -v+l;
						h.normalize();
						double ndotl = rec.N*l;
						double ndoth = rec.N*h;
						double maxndotl = std::max(0.0,ndotl);
						double maxndothp = pow(std::max(0.0,ndoth),rec.phong_exp);

						lamr += rec.kd.r*intensityr*maxndotl;
						lamg += rec.kd.g*intensityg*maxndotl;
						lamb += rec.kd.b*intensityb*maxndotl;

						bpr += rec.ks.r*intensityr*maxndothp;
						bpg += rec.ks.g*intensityg*maxndothp;
						bpb += rec.ks.b*intensityb*maxndothp;
					}
					delete sr;
				}
			}

		}

		// Specular reflection
		if((rec.km.r!=0 || rec.km.g!=0 || rec.km.b!=0) && depth<=this->max_depth) {
			Vector spec_reflect_dir = reflect(v, rec.N);
			if(rec.glossy==0) {
				Ray *ref_ray = new Ray(p, spec_reflect_dir);
				Imf::Rgba reflectColor = this->rayColor(ref_ray, EPSILON, infinity, depth);
				reflectr = rec.km.r*reflectColor.r;
				reflectg = rec.km.g*reflectColor.g;
				reflectb = rec.km.b*reflectColor.b;
				delete ref_ray;
			}
			else if(rec.glossy>0 && depth < 2) {
				//Make a square orthogonal to reflect direction and randomly sample it
				Vector u = spec_reflect_dir^Camera::UP;
				Vector v = u^spec_reflect_dir;
				double term1 = -0.5*rec.glossy;
				u.normalize(); v.normalize(); spec_reflect_dir.normalize();
				for(int s=0; s < this->glossy_samples; ++s) {
					double rand_u = ((double) rand()) / RAND_MAX;
					double rand_v = ((double) rand()) / RAND_MAX;
					double ur = term1+rand_u*rec.glossy;
					double vr = term1+rand_v*rec.glossy;
					Vector random_reflect_dir = spec_reflect_dir+ur*u+vr*v;
					Ray *randomReflectRay = new Ray(p, random_reflect_dir);
					Imf::Rgba reflectColor = this->rayColor(randomReflectRay, EPSILON, infinity, depth);
					reflectr += rec.km.r*reflectColor.r/this->glossy_samples;
					reflectg += rec.km.g*reflectColor.g/this->glossy_samples;
					reflectb += rec.km.b*reflectColor.b/this->glossy_samples;
					delete randomReflectRay;
				}
			}
		}

		// Refraction
		double R = 1;
		if(rec.refract_index!=0 && depth<=this->max_depth) {
			Vector t;
			double c;
			bool refract = true;
			if(v*rec.N < 0) {
				this->refract(v, rec.N, 1, rec.refract_index, &t);
				c = -v*rec.N;
			}
			else {
				kr = exp(-rec.kc.r*thit);
				kg = exp(-rec.kc.g*thit);
				kb = exp(-rec.kc.b*thit);
				if(this->refract(v, rec.N, rec.refract_index, 1, &t))
					c = t*rec.N;
				else refract = false;
			}
			if(refract) {
				Ray *refractRay = new Ray(p, t);
				refractRay->origin = refractRay->getPoint(EPSILON);
				refractRay->normalize();
				Imf::Rgba refractColor = this->rayColor(refractRay, 0, infinity, depth);
				refractr = refractColor.r;
				refractg = refractColor.g;
				refractb = refractColor.b;
				delete refractRay;

				double refractMinus1 = rec.refract_index-1;
				double refractPlus1 = rec.refract_index+1;
				double R0 = (refractMinus1*refractMinus1)/(refractPlus1*refractPlus1);
				R = R0 + (1-R0)*pow((1-c),5);
			}
		}

		ray_color.r += kr*(ambr+lamr+bpr+R*reflectr+(1-R)*refractr);
		ray_color.g += kg*(ambg+lamg+bpg+R*reflectg+(1-R)*refractg);
		ray_color.b += kb*(ambb+lamb+bpb+R*reflectb+(1-R)*refractb);
		ray_color.a = 1.0;

	}
	return ray_color;
}

Vector Scene::reflect(const Vector& d, const Vector& n) {
	Vector r = d-2*(d*n)*n;
	return r;
}

bool Scene::refract(const Vector &d, const Vector &n, double refract_index1, double refract_index2, Vector *t) {
	double index_ratio = refract_index1/refract_index2;
	double dDotN = -d*n;
	double root_term_enum = 1.0-index_ratio*index_ratio*(1.0-dDotN*dDotN);
	if (root_term_enum < 0.0) return false;
	*t = index_ratio*(d-n*(dDotN))-n*sqrt(root_term_enum);
	return true;
}
