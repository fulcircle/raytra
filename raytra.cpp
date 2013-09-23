//============================================================================
// Name        : raytra.cpp
// Author      : Vikram Sunkavalli
// Main function.
//============================================================================

#include "camera.h"
#include "globals.h"
#include "scene.h"
#include "readscene.h"

using namespace std;

Camera *cam;
vector<Surface*> *surfaces = new vector<Surface*>;
Lights *lights = new Lights;
vector<double> *opts = new vector<double>;

int main(int argc, char *argv[]) {

	if ( argc < 3) {
	    cout<<"Usage: "<< argv[0] <<" <scene file> <output file>\n";
	}
	 else {
		 parseSceneFile(argv[1]);
		 Scene scene(cam, surfaces, lights, opts);
		 Utility::renderToFile(argv[2], scene.render());
	   }

	return 0;

}
