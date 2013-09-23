
#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <vector>
class Camera;
typedef struct _lights Lights;
class Surface;

extern Camera *cam;
extern std::vector<Surface*> *surfaces;
extern Lights *lights;
extern std::vector<double> *opts;

//#define IM_DEBUGGING

#define EPSILON 1e-6

#define DEFAULT_MAX_DEPTH 3
#define DEFAULT_SAMPLES 1
#define DEFAULT_SHADOW_SAMPLES 1
#define DEFAULT_GLOSSY_SAMPLES 1

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#endif
