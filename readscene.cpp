#include "readscene.h"


using namespace std;


// this is called from the parseSceneFile function, which uses
// it to get the float from the correspoding position on the line.
//
// return the corresponding token in the inString. Errors out
// if you've asked for more than are in the line.
//
// you really don't need to know what is going on in here, I think.
//
float getTokenAsFloat (string inString, int whichToken) {

	float thisFloatVal = 0.;    // the return value

	if (whichToken == 0) {
		cerr << "error: the first token on a line is a character!" << endl;
		exit (-1);
	}

	// c++ string class has no super-easy way to tokenize, let's use c's:
	char *cstr = new char [inString.size () + 1];

	strcpy (cstr, inString.c_str());

	char *p = strtok (cstr, " ");
	if (p == 0) {
		cerr << "error: the line has nothing on it!" << endl;
		exit (-1);
	}

	for (int i = 0; i < whichToken; i++) {
		p = strtok (0, " ");
		if (p == 0) {
			cerr << "error: the line is not long enough for your token request!" << endl;
			exit (-1);
		}
	}

	thisFloatVal = atof (p);

	delete[] cstr;

	return thisFloatVal;
}


//
// read the scene file.
//
// You'll need a few globals (or add arguments to this function): for the
// list of surfaces (geometric objects like spheres, triangles, planes) and
// another for the lights. These can be of type std::vector. You'll also
// need a global (or other parameter) for the camera.
//
// This should be pretty self-explanatory: it reads through the lines, and
// using the first character figures out what kind of command it is. It
// then calls the "getTokenAsFloat" routine to pull out the needed values.
// NOTE: since different commands take different number of arguments, you
// must take care not to call getTokenAsFloat for an index that is beyond the
// end of the line!
//
// One tricky bit: when an surface (sphere/triangle/plane) is read in, we want
// to keep a pointer to it so that if the VERY NEXT LINE is a material definition,
// we can add that material to the object. In the code that follows, I use the
// variable "lastSurfaceLoaded" to do that, but the code is commented out since
// I don't know the class names you will be using.
//
// Very minimal error check here. You might improve it slightly, but we'll
// only use "correct" scene files.
//
//
void parseSceneFile (char *filnam)
{

	ifstream inFile(filnam);    // open the file
	string line;

	if (! inFile.is_open ()) {  // if it's not open, error out.
		cerr << "can't open scene file" << endl;
		exit (-1);
	}


	// Note: you'll have to keep track of whatever the last surface (geometry) type
	// you loaded in was, so you can apply a material to it if "m..." is on the
	// next line. So here, you'll have something like:
	//
	Surface *lastSurfaceLoaded = 0;
	//
	// and each time you load in a new piece of geometry (sphere, triangle, plane)
	// you will set lastSurfaceLoaded to it.


	while (! inFile.eof ()) {   // go through every line in the file until finished

		getline (inFile, line); // get the line

		switch (line[0])  {     // we'll decide which command based on the first character

		//
		// geometry types:
		//
		// NOTE: whichever type of geo you load in, set the "lastSurfaceLoaded"
		// pointer to it, so we can add a material spec to it if a material is
		// specified on the next line
		//
		case 's': {
			// it's a sphere, load in the parameters

			float x, y, z, r;
			x = getTokenAsFloat (line, 1);
			y = getTokenAsFloat (line, 2);
			z = getTokenAsFloat (line, 3);
			r = getTokenAsFloat (line, 4);

			// build your sphere here from the parameters
			// i.e. you must call your sphere constructor and set its position
			// and radius from the above values. You must also put your new
			// sphere into the objects list (which can be global)
			// So something like;
			// mySphereClass *ms = new mySphereClass (x, y, z, r);   // make a new instance of your sphere class
			// objectsList->push_back (ms);  // objectsList is a global std:vector for example.
			// lastSurfaceLoaded = ms;       // this sphere was the last surface loaded (for materials)
			Sphere *s = new Sphere(Point(x,y,z), r);
			surfaces->push_back(s);
			lastSurfaceLoaded = s;

#ifdef IM_DEBUGGING
			// if we're debugging, show what we got:
			cout << "[Sphere] ";
			cout << "parameters: x: " << x << " y: " << y << " z: " << z << " r: " << r << endl;
#endif
		}

		break;

		case 't':  { // triangle
			double ax, ay, az, bx, by, bz, cx, cy, cz;
			ax = getTokenAsFloat (line, 1);
			ay = getTokenAsFloat (line, 2);
			az = getTokenAsFloat (line, 3);

			bx = getTokenAsFloat (line, 4);
			by = getTokenAsFloat (line, 5);
			bz = getTokenAsFloat (line, 6);

			cx = getTokenAsFloat (line, 7);
			cy = getTokenAsFloat (line, 8);
			cz = getTokenAsFloat (line, 9);

			Triangle *t = new Triangle(Point(ax, ay, az), Point(bx, by, bz), Point(cx, cy, cz));
			surfaces->push_back(t);
			lastSurfaceLoaded = t;
#ifdef IM_DEBUGGING
			cout << "[Triangle] vertices: a: " << t->a << " b: " << t->b << " c: " << t->c << endl;
#endif
		}

		break;

		case 'p': {  // plane
			double nx, ny, nz, d;
			nx = getTokenAsFloat (line, 1);
			ny = getTokenAsFloat (line, 2);
			nz = getTokenAsFloat (line, 3);
			d = getTokenAsFloat (line, 4);
			Vector N(nx, ny, nz);
			N.normalize();
			Plane *p = new Plane(N, d);
			surfaces->push_back(p);
			lastSurfaceLoaded = p;
#ifdef IM_DEBUGGING
			cout << "[Plane] normal: " << N << " d: " << d << endl;
#endif
		}
		break;

		//
		// camera:
		//
		case 'c':  { // camera
			double posx = getTokenAsFloat (line, 1);
			double posy = getTokenAsFloat (line, 2);
			double posz = getTokenAsFloat (line, 3);

			double dirx = getTokenAsFloat (line, 4);
			double diry = getTokenAsFloat (line, 5);
			double dirz = getTokenAsFloat (line, 6);

			double fl = getTokenAsFloat (line, 7);

			double lr = getTokenAsFloat (line, 8);
			double tb = getTokenAsFloat (line, 9);

			// one trick here: the cameras pixel count (width, height) are integers,
			// so cast them.

			int width = getTokenAsFloat (line, 10);
			int height = getTokenAsFloat (line, 11);

			Point e(posx, posy, posz);
			Vector v(dirx,diry,dirz);

			cam = new Camera(e, v, lr, tb, fl, width, height);

#ifdef IM_DEBUGGING
			cout << "[Camera] position: x: " << posx << " y: " << posy << " z: " << posz << endl;
			cout << "[Camera] orientation: dirx: " << dirx << " diry: " << diry << " dirz: " << dirz << endl;
			cout << "[Camera] f:" << fl << "mm, " << lr << "mmx" << tb << "mm, " << width << "x" << height << endl;
#endif
		}
		break;


		//
		// lights:
		//
		case 'l':   // light

			// slightly different from the rest, we need to examine the second param,
			// which is at the third position on the line:
			switch (line[2]) {
			case 'p': {  // point light
				double posx = getTokenAsFloat (line, 2);
				double posy = getTokenAsFloat (line, 3);
				double posz = getTokenAsFloat (line, 4);

				double ir = getTokenAsFloat (line, 5);
				double ig = getTokenAsFloat (line, 6);
				double ib = getTokenAsFloat (line, 7);

				Imf::Rgba c;
				c.r = ir;  c.g = ig; c.b = ib; c.a = 1.0;
				PointLight *l = new PointLight(Point(posx, posy, posz), c);
				lights->pointLights.push_back(l);
#ifdef IM_DEBUGGING
				cout << "[PointLight] position: x: " << posx << " y: " << posy << " z: " << posz << endl;
				cout << "[PointLight] color: r: " << ir << " g: " << ig << " b: " << ib << endl;
#endif
			}
			break;
			case 'r': {  // area light
				Vector edgeVectorA, edgeVectorB;
				Point position;
				Imf::Rgba intensity;

				edgeVectorA.x = getTokenAsFloat(line, 2);
				edgeVectorA.y = getTokenAsFloat(line, 3);
				edgeVectorA.z = getTokenAsFloat(line, 4);

				edgeVectorB.x = getTokenAsFloat(line, 5);
				edgeVectorB.y = getTokenAsFloat(line, 6);
				edgeVectorB.z = getTokenAsFloat(line, 7);

				position.x = getTokenAsFloat(line, 8);
				position.y = getTokenAsFloat(line, 9);
				position.z = getTokenAsFloat(line, 10);

				intensity.r = getTokenAsFloat(line, 11);
				intensity.g = getTokenAsFloat(line, 12);
				intensity.b = getTokenAsFloat(line, 13);

				AreaLight *areaLight = new AreaLight(edgeVectorA, edgeVectorB, position, intensity);
				lights->areaLights.push_back(areaLight);

#ifdef IM_DEBUGGING
				cout << "[DirLight] Position: " << position << " Vector A: " << edgeVectorA << " Vector B: " << edgeVectorB << endl;
				cout << "[DirLight] intensity: r: " << intensity.r << " g: " << intensity.g << " b: " << intensity.b << endl;
#endif
			}
			break;
			case 'a':   { // ambient light
				double ar = getTokenAsFloat (line, 2);
				double ag = getTokenAsFloat (line, 3);
				double ab = getTokenAsFloat (line, 4);

				Imf::Rgba c;
				c.r = ar; c.g = ag; c.b = ab;
				AmbientLight *l = new AmbientLight(c);
				lights->ambientLights.push_back(l);
#ifdef IM_DEBUGGING
				cout << "[AmbientLight] r: " << ar << " g: " << ag << " b: " << ab << endl;
#endif
			}
			break;

			}

			break;

			//
			// materials:
			//
			case 'm':   {// material
				// the trick here: we keep a pointer to the last surface we read in,
				// so we can apply this material to it. Say it's called "lastSurfaceLoaded"
				// we migh then do something like this:
				//
				//  1. read in the 10 material parameters: dr, dg, db, sr, sg, sb, r, ir, ig, ib
				//  2. call lastSurfaceLoaded->setMaterial(dr, dg, db,...);
				//
				Imf::Rgba diffuse, specular, ideal;
				diffuse.r = getTokenAsFloat (line, 1);
				diffuse.g = getTokenAsFloat (line, 2);
				diffuse.b = getTokenAsFloat (line, 3);
				diffuse.a = 1;

				specular.r = getTokenAsFloat (line, 4);
				specular.g = getTokenAsFloat (line, 5);
				specular.b = getTokenAsFloat (line, 6);
				specular.a = 1;

				double phong_exp = getTokenAsFloat (line, 7);

				ideal.r = getTokenAsFloat (line, 8);
				ideal.g = getTokenAsFloat (line, 9);
				ideal.b = getTokenAsFloat (line, 10);
				ideal.a = 1;

				if(lastSurfaceLoaded!=NULL)
					lastSurfaceLoaded->setMaterial(diffuse, specular, phong_exp, ideal);
			}
			break;

			case 'r': {
				Imf::Rgba absorption_coefficient;
				double refract_index = getTokenAsFloat (line, 1);

				absorption_coefficient.r = getTokenAsFloat (line, 2);
				absorption_coefficient.g = getTokenAsFloat (line, 3);
				absorption_coefficient.b = getTokenAsFloat (line, 4);


				if(lastSurfaceLoaded!=NULL)
					lastSurfaceLoaded->setRefract(refract_index, absorption_coefficient);
			}
			break;

			case 'g': {

				double glossy = getTokenAsFloat (line, 1);

				if(lastSurfaceLoaded!=NULL)
					lastSurfaceLoaded->setGlossy(glossy);
			}
			break;


			case '/':
				// don't do anything, it's a comment
				break;


				//
				// options
				//
			case 'o':  { // make your own options if you wish
				double max_depth = getTokenAsFloat (line, 1);
				double samples = getTokenAsFloat(line, 2);
				double shadow_samples = getTokenAsFloat(line, 3);
				double glossy_samples = getTokenAsFloat(line, 4);
				opts->push_back(max_depth);
				opts->push_back(samples);
				opts->push_back(shadow_samples);
				opts->push_back(glossy_samples);
#ifdef IM_DEBUGGING
				cout << "[Options] Max Depth: " << max_depth << " Samples: " << samples << " Shadow Samples: " << shadow_samples << " Glossy Samples: " << glossy_samples << endl;
#endif
			}
				break;
		}

	}
}
